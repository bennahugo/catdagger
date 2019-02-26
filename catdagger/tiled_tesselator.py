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
from catdagger.filters import within_radius_from, \
    notin, arealess, skewness_more, pos2neg_more
from catdagger.geometry import BoundingBox, BoundingConvexHull, merge_regions
from catdagger.fits_tools import FitsStokesTypes
log = logger.getLogger("tiled_tesselator")

def tag_regions(stokes_cube,  
                regionsfn = "dE.reg", 
                sigma = 2.3, 
                block_size=80, 
                hdu_id = 0, 
                use_stokes="I", 
                global_stat_percentile=30.0,
                min_blocks_in_region = 3,
                min_distance_from_centre = 0,
                exclusion_zones=[],
                max_right_skewness=np.inf,
                max_abs_skewness=np.inf,
                max_positive_to_negative_flux=np.inf):
    """
        Tiled tesselator

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
    print>>log, "Stokes in the cube: {0:s}".format(",".join([reverse_stokes_map[s] for s in stokes_axis]))
    sel_stokes = [reverse_stokes_map[s] for s in stokes_axis].index(use_stokes)
    print>>log, "Stokes slice selected: {0:d} (Stokes {1:s})".format(sel_stokes, use_stokes)
    sel_stokes = np.take(cube, sel_stokes, axis=(hdr["NAXIS"] - types["STOKES"]))
    chan_axis = hdr["NAXIS"] - types["FREQ"] if types["FREQ"] > types["STOKES"] else hdr["NAXIS"] - types["FREQ"] - 1
    print>>log, "Collapsing axis: {0:d}".format(types["FREQ"])
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
    segment_cutoff = percentile_stat * sigma
    print>>log, "Computed regional statistics (global std of {0:.2f} mJy)".format(percentile_stat * 1.0e3)
    tagged_regions = []
    for (y, x) in np.argwhere(binned_stats > segment_cutoff):
        det = binned_stats[y, x] / float(percentile_stat)
        reg_name = "reg[{0:d},{1:d}]".format(x, y)
        tagged_regions.append(BoundingBox(bin_lower[x], bin_upper[x], 
                                          bin_lower[y], bin_upper[y], 
                                          det, reg_name, w, band_avg))
    
    if min_distance_from_centre > 0:
        print>>log, "Enforsing radial exclusion zone of {0:.2f} px form " \
                    "phase tracking centre".format(min_distance_from_centre)
        exclusion_zones.append((hdr["CRPIX{0:d}".format(types["RA---SIN"])],
                                hdr["CRPIX{0:d}".format(types["DEC--SIN"])],
                                float(min_distance_from_centre)))

    # enforce all exclusion zones
    print>>log, "Enforsing exclusion zones:"
    for (cx, cy, exclrad) in exclusion_zones:
        tagged_regions = filter(notin(filter(within_radius_from(exclrad, cx, cy), 
                                             tagged_regions)), 
                                tagged_regions)

    print>>log, "Merging regions" 
    tagged_regions = [i for i in merge_regions(tagged_regions,  
                                               exclusion_zones=exclusion_zones)] 
    # apply regional filters
    print>>log, "Culling regions based on filtering criteria:"
    min_area=min_blocks_in_region * block_size**2
    tagged_regions = filter(notin(filter(arealess(min_area=min_area), 
                                         tagged_regions)), 
                            tagged_regions)
    tagged_regions = filter(notin(filter(skewness_more(max_skewness=max_right_skewness,
                                                      absskew=False), 
                                         tagged_regions)),
                            tagged_regions)
    tagged_regions = filter(notin(filter(skewness_more(max_skewness=max_abs_skewness,
                                                       absskew=True), 
                                         tagged_regions)),
                            tagged_regions)
    tagged_regions = filter(notin(filter(pos2neg_more(max_positive_to_negative_flux), 
                                         tagged_regions)),
                            tagged_regions)

    # finally we're done
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