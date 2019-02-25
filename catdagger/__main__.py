#!/usr/bin/env python

import copy
import numpy as np
from astropy.io import fits
from astropy import wcs
import argparse
import scipy.stats as sstats
import scipy.signal as ssig
import scipy.spatial as spat
import Tigger
from catdagger import logger
log = logger.getLogger("CATDagger")

'''
The following definition can be found in Table 28,
Definition of the Flexible Image Transport System (FITS),   version 3.0
W. D.  Pence, L.  Chiappetti, C. G.  Page, R. A.  Shaw, E.  Stobie
A&A 524 A42 (2010)
DOI: 10.1051/0004-6361/201015362
'''
FitsStokesTypes = {
    "I" : 1, #Standard Stokes unpolarized
    "Q" : 2, #Standard Stokes linear
    "U" : 3, #Standard Stokes linear
    "V" : 4, #Standard Stokes circular
    "RR": -1, #Right-right circular
    "LL": -2, #Left-left circular
    "RL": -3, #Right-left cross-circular
    "LR": -4, #Left-right cross-circular
    "XX": -5, #X parallel linear
    "YY": -6, #Y parallel linear
    "XY": -7, #XY cross linear
    "YX": -8  #YX cross linear
}

class BoundingConvexHull():
    def __init__(self, list_hulls, sigma, name, wcs=None):
        self._wcs = wcs
        self._name = name
        self._vertices = points = np.vstack([b.corners
            if hasattr(b, "corners") else [b[0], b[1]] for b in list_hulls])
        self._hull = spat.ConvexHull(points)
        self._sigma = sigma

    def __str__(self):
        return "{0:.2f}x within region ".format(self._sigma) + \
               ",".join(["({0:d},{1:d})".format(x,y) for (x,y) in self.corners])
    @property
    def area(self):
        lines = np.hstack([self.corners, np.roll(self.corners, -1, axis=0)])
        return 0.5 * np.abs(np.sum([x1*y2-x2*y1 for x1,y1,x2,y2 in lines]))
    @property
    def name(self):
        return self._name

    @property
    def area_sigma(self):
        return self._sigma

    @property
    def wcs(self):
        return self._wcs

    @property
    def corners(self):
        """ Returns vertices and guarentees clockwise winding """
        return self._vertices[self._hull.vertices]

    def normals(self, left = True):
        """ return a list of left normals to the hull """
        normals = []
        for i in xrange(self.corners.shape[0]):
            # assuming clockwise winding
            j = (i + 1) % self.corners.shape[0]
            edge = self.corners[j, :] - self.corners[i, :]
            if left:
                normals.append((-edge[1], edge[0]))
            else:
                normals.append((edge[1], -edge[0]))
        return np.asarray(normals, dtype=np.double)

    @property
    def lnormals(self):
        return self.normals(left = True)

    @property
    def rnormals(self):
        return self.normals(left=False)
    
    def is_neighbour(self, other, min_sep_dist=1.0e-4):
        """ 
            Implements the separating lines collision detection theorem 
        """
        if not isinstance(other, BoundingConvexHull):
            raise TypeError("rhs must be a BoundingConvexHull")

        # get the projection axes
        normals = np.vstack([self.lnormals, other.lnormals])
        norms = np.linalg.norm(normals, axis=1)
        normals = normals / norms[None, 2]

        # compute vectors to corners from origin
        vecs_reg1 = self.corners
        vecs_reg2 = other.corners

        # compute projections onto normals
        for ni, n in enumerate(normals):
            projs = np.dot(vecs_reg1, n.T)
            minproj_reg1 = np.min(projs)
            maxproj_reg1 = np.max(projs)
            projs = np.dot(vecs_reg2, n.T)
            minproj_reg2 = np.min(projs)
            maxproj_reg2 = np.max(projs)
            if minproj_reg2 - maxproj_reg1 > 1.0e-4 or minproj_reg1 - maxproj_reg2 > 1.0e-4:
                return False
        return True

    @property
    def centre(self):
        # Barycentre of polygon
        return np.mean(self._vertices, axis=1)

    def __contains__(self, s):
        if not isinstance(s, Tigger.Models.SkyModel.Source):
            raise TypeError("Source must be a Tigger lsm source")
        ra = np.rad2deg(s.pos.ra)
        dec = np.rad2deg(s.pos.dec)
        x, y, _, _ = self._wcs.all_world2pix([[ra, dec, 0, 0]], 1)[0]

        dot = 0
        for i in range(len(self.corners)):
            j = (i + 1) % len(self.corners)
            v1 = self.corners[i] - np.array([x, y])
            v2 = self.corners[j] - np.array([x, y])
            dot += np.arccos(np.clip(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)), -1, +1))
        return np.abs(360 - np.rad2deg(dot)) < 1.0e-6

class BoundingBox(BoundingConvexHull):
    def __init__(self, xl, xu, yl, yu, sigma, name, wcs=None):
        BoundingConvexHull.__init__(self,
                                    [[xl,yl],[xl,yu],[xu,yu],[xu,yl]],
                                    sigma,
                                    name,
                                    wcs)
def merge_regions(regions, min_sep_distance=1.0e-4, min_area=0):
    """ Merge neigbouring regions into convex hulls """
    for reg in regions:
        if not isinstance(reg, BoundingConvexHull):
           raise TypeError("Expected BoundingConvexHull as argument")

    merged = True
    orig_regs = len(regions)
    while merged:
        merged = False
        new_regions = []
        exclude_list = []
        for me_i in range(len(regions)):
            me = regions[me_i]
            if me in exclude_list: continue # already merged with another region
            nreg = [me]
            for other_i in range(me_i + 1, len(regions)):
                other = regions[other_i]
                if me.is_neighbour(other, min_sep_distance):
                    merged = True
                    exclude_list.append(other)
                    nreg.append(other)
                    print>>log, "\t - Merged regions {0:s} and {1:s}".format(me.name, other.name)
            new_regions.append(BoundingConvexHull(nreg,
                                                  sigma=np.mean([reg.area_sigma for reg in nreg]),
                                                  name="&".join([reg.name for reg in nreg]),
                                                  wcs=me.wcs))
        regions = new_regions
    discard = []
    for reg in regions:
        if reg.area < min_area:
            discard.append(reg)
            print>>log, "\t - Discarding region {0:s} because of its small size".format(reg.name)
    for reg in discard:
        regions.remove(reg)
    return regions

def tag_regions(stokes_cube,  
                regionsfn = "dE.reg", 
                sigma = 2.3, 
                block_size=80, 
                hdu_id = 0, 
                use_stokes="I", 
                global_stat_percentile=30.0,
                min_blocks_in_region = 3):
    """
        Method to tag regions with higher than sigma * percentile noise
    """
    with fits.open(stokes_cube) as img:
        cube = img[hdu_id].data
        hdr = img[hdu_id].header
        w = wcs.WCS(hdr)
    types = {hdr["CTYPE{0:d}".format(ax + 1)]: (ax + 1) for ax in range(hdr["NAXIS"])}
    if set(types.keys()) != set(["FREQ", "STOKES", "RA---SIN", "DEC--SIN"]):
        raise TypeError("FITS must have FREQ, STOKES and RA and DEC ---SIN axes")
    stokes_axis = np.arange(hdr["CRVAL{0:d}".format(types["STOKES"])] - hdr["CRPIX{0:d}".format(types["STOKES"])] * (hdr["CDELT{0:d}".format(types["STOKES"])] - 1),
                            (hdr["NAXIS{0:d}".format(types["STOKES"])] + 1) * hdr["CDELT{0:d}".format(types["STOKES"])],
                            hdr["CDELT{0:d}".format(types["STOKES"])])
    reverse_stokes_map = {FitsStokesTypes[k]: k for k in FitsStokesTypes.keys()}
    print>>log, "Stokes in the cube:", [reverse_stokes_map[s] for s in stokes_axis]
    sel_stokes = [reverse_stokes_map[s] for s in stokes_axis].index(use_stokes)
    print>>log, "Stokes slice selected:", sel_stokes
    sel_stokes = np.take(cube, sel_stokes, axis=(hdr["NAXIS"] - types["STOKES"]))
    chan_axis = hdr["NAXIS"] - types["FREQ"] if types["FREQ"] > types["STOKES"] else hdr["NAXIS"] - types["FREQ"] - 1
    print>>log, "Collapsing axis:", types["FREQ"]
    band_avg = np.mean(sel_stokes, axis=chan_axis)
    bin_lower = np.arange(0, hdr["NAXIS{0:d}".format(types["RA---SIN"])], block_size)
    bin_upper = np.clip(bin_lower + block_size, 0, hdr["NAXIS{0:d}".format(types["RA---SIN"])])
    assert bin_lower.shape == bin_upper.shape
    if band_avg.shape[0] != band_avg.shape[1]:
        raise TypeError("Image must be square!")
    print>>log, "Creating regions of {0:d} px".format(block_size)
    binned_stats = np.zeros((bin_lower.shape[0],
                             bin_lower.shape[0]))
    for y, (ly, uy) in enumerate(zip(bin_lower, bin_upper)):
        for x, (lx, ux) in enumerate(zip(bin_lower, bin_upper)):
            wnd = band_avg[ly:uy, lx:ux].flatten()
            binned_stats[y, x] = np.std(wnd)
    percentile_stat = np.nanpercentile(binned_stats, global_stat_percentile)
    print>>log, "Computed regional statistics (global std of {0:.2f} mJy)".format(percentile_stat * 1.0e3)
    tagged_regions = []
    for (y, x) in np.argwhere(binned_stats > percentile_stat * sigma):
        det = binned_stats[y, x] / float(percentile_stat)
        reg_name = "reg[{0:d},{1:d}]".format(x, y)
        tagged_regions.append(BoundingBox(bin_lower[x], bin_upper[x], bin_lower[y], bin_upper[y], det, reg_name, wcs=w))

    print>>log, "Merging regions" 
    tagged_regions = [i for i in merge_regions(tagged_regions, min_area=min_blocks_in_region * block_size**2)] 
    with open(regionsfn, "w+") as f:
        f.write("# Region file format: DS9 version 4.0\n")
        f.write("global color=red font=\"helvetica 6 normal roman\" edit=1 move=1 delete=1 highlite=1 include=1 wcs=wcs\n")
        for reg in tagged_regions:
            f.write("physical; polygon({0:s}) #select=1 text={1:s}\n".format(",".join(map(str, reg.corners.flatten())),
                                                                             "{mean area deviation %.2fx}" % reg._sigma))
        print>>log, "Writing dE regions to DS9 regions file {0:s}".format(regionsfn)
    print>>log, "The following regions must be tagged for dEs ({0:.2f}x{1:.2f} mJy)".format(sigma, percentile_stat * 1.0e3)
    if len(tagged_regions) > 0:
        for r in tagged_regions:
            print>>log, "\t - {0:s}".format(str(r))
    else:
        print>>log, "No regions met cutoff criterion. No dE tags shall be raised."
    return tagged_regions

def tag_lsm(lsm,
            stokes_cube,
            tagged_regions,
            hdu_id=0,
            regionsfn = "dE.srcs.reg",
            taggedlsm_fn="tagged.catalog.lsm.html",
            de_tag="dE"):
    with fits.open(stokes_cube) as img:
        cube = img[hdu_id].data
        hdr = img[hdu_id].header
        w = wcs.WCS(hdr)

    with open(regionsfn, "w+") as f:
        f.write("# Region file format: DS9 version 4.0\n")
        f.write("global color=green font=\"helvetica 6 normal roman\" edit=1 move=1 delete=1 highlite=1 include=1 wcs=wcs\n")

        mod = Tigger.load(lsm)
        for ireg, reg in enumerate(tagged_regions):
            print>>log, "Tagged sources in Region {0:d}:".format(ireg), str(reg)
            encircled_sources = filter(lambda s: s in reg, mod.sources)
            encircled_fluxes = [s.flux.I for s in encircled_sources]
            for s in encircled_sources:
                s.setTag("dE", True)
                s.setTag("cluster", reg.name) #recluster sources
            if len(encircled_fluxes) > 0:
                argmax = np.argmax(encircled_fluxes)
                s = encircled_sources[argmax]
                s.setTag("cluster_lead", True)
                ra = np.rad2deg(s.pos.ra)
                dec = np.rad2deg(s.pos.dec)
                x, y, _, _ = w.all_world2pix([[ra, dec, 0, 0]], 1)[0]
                x = int(x)
                y = int(y)
                f.write("physical;circle({0:d}, {1:d}, 20) # select=1 text={2:s}\n".format(x, y,
                        "{%.2f mJy}" % (s.flux.I * 1.0e3)))
                print>>log, "\t - {0:s} tagged as '{1:s}' cluster lead".format(s.name, de_tag)
        print>>log, "Writing tagged leads to DS9 regions file {0:s}".format(regionsfn)
    print>>log, "Writing tagged LSM to {0:s}".format(taggedlsm_fn)
    mod.save(taggedlsm_fn)

def main():
    parser = argparse.ArgumentParser("CATDagger - an automatic differential gain tagger (C) SARAO, Benjamin Hugo 2019")
    parser.add_argument("noise_map",
                        type=str,
                        help="Residual / noise FITS map to use for estimating local RMS")
    parser.add_argument("--stokes",
                        type=str,
                        default="I",
                        help="Stokes to consider when computing global noise estimates. Ideally this should be 'V', if available")
    parser.add_argument("--min-tiles-region",
                        type=int,
                        default=3,
                        help="Minimum number of tiles per region. Regions with fewer tiles will not be tagged as dE")
    parser.add_argument("--input-lsm",
                        type=str,
                        default=None,
                        help="Tigger LSM to recluster and tag. If this is not specified only DS9 regions will be written out")
    parser.add_argument("--ds9-reg-file",
                        type=str,
                        default="dE.reg",
                        help="SAODS9 regions filename to write out")
    parser.add_argument("--ds9-tag-reg-file",
                        type=str,
                        default="dE.tags.reg",
                        help="SAODS9 regions filename to contain tagged cluster leads as circles")
    parser.add_argument("-s",
                        "--sigma",
                        type=float,
                        default=2.3,
                        help="Threshold to use in detecting outlier regions")
    parser.add_argument("--tile-size",
                        type=int,
                        default=80,
                        help="Number of pixels per region tile axis")
    parser.add_argument("--global-rms-percentile",
                        type=float,
                        default=30,
                        help="Percentile tiles to consider for global rms calculations")
    parser.add_argument("--de-tag-name",
                        type=str,
                        default="dE",
                        help="Tag name to use for tagged sources in tigger LSM")
    args = parser.parse_args()
    tagged_regions = tag_regions(args.noise_map,
                                 regionsfn = args.ds9_reg_file,
                                 sigma = args.sigma,
                                 block_size = args.tile_size,
                                 hdu_id = 0,
                                 use_stokes = args.stokes,
                                 global_stat_percentile = args.global_rms_percentile,
                                 min_blocks_in_region = args.min_tiles_region)
    if args.input_lsm is not None:
        tag_lsm(args.input_lsm,
                args.noise_map,
                tagged_regions,
                hdu_id=0,
                regionsfn = args.ds9_tag_reg_file,
                taggedlsm_fn=args.input_lsm + ".de_tagged.lsm.html",
                de_tag=args.de_tag_name)

if __name__ == "__main__":
    main()

