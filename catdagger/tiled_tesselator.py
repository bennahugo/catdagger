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
from catdagger.fits_tools import FitsStokesTypes, read_stokes_slice, getcrpix
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
    fn = stokes_cube
    w, hdr, band_avg = read_stokes_slice(stokes_cube, hdu_id, use_stokes, average_channels=True)
    bin_lower = np.arange(0, band_avg.shape[0], block_size)
    bin_upper = np.clip(bin_lower + block_size, 0, band_avg.shape[0])
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
        crra, crdec = getcrpix(fn, hdu_id, use_stokes)
        exclusion_zones.append((crra,
                                crdec,
                                float(min_distance_from_centre)))

    # enforce all exclusion zones
    print>>log, "Enforsing exclusion zones:"
    for (cx, cy, exclrad) in exclusion_zones:
        tagged_regions = filter(notin(filter(within_radius_from(exclrad, cx, cy), 
                                             tagged_regions)), 
                                tagged_regions)
    if len(exclusion_zones) == 0: 
        print>>log, "\t - No exclusion zones"
    print>>log, "Merging regions:" 
    prev_tagged_regions = copy.deepcopy(tagged_regions)
    tagged_regions = [i for i in merge_regions(tagged_regions,  
                                               exclusion_zones=exclusion_zones)]
    if prev_tagged_regions == tagged_regions: 
        print>>log, "\t - No mergers" 
    # apply regional filters
    print>>log, "Culling regions based on filtering criteria:"
    prev_tagged_regions = copy.deepcopy(tagged_regions)
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
    if prev_tagged_regions == tagged_regions: 
        print>>log, "\t - No cullings"
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
        print>>log, "\t - No regions met cutoff criteria. No dE tags shall be raised."
    return tagged_regions