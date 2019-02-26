#!/usr/bin/env python

# CATDagger: an automatic differential gain catalog tagger
# (c) 2019 South African Radio Astronomy Observatory, B. Hugo
# This code is distributed under the terms of GPLv2, see LICENSE.md for details
# -*- coding: utf-8 -*-
#
# Copyright (c) 2019 SARAO
#
# This file is part of CATDagger.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import argparse

from catdagger import logger
from catdagger.tiled_tesselator import tag_regions
from catdagger.lsm_tools import tag_lsm
import logging
logging.getLogger("matplotlib").disabled=True
logger.setGlobalVerbosity(["matplotlib=40"])
log = logger.getLogger("__main__")

def exclusion_zone(val):
    if not isinstance(val, str):
        raise argparse.ArgumentTypeError("Exclusion zone must be of type string")
    valtup = val.split(",") if isinstance(val, str) else tuple(val) if isinstance(val, list) else []

    if len(valtup) != 3:
        raise argparse.ArgumentTypeError("Exclusion zone must be a tripple like (cx, cy, radius)")
    try:
        cx = int(valtup[0]); cy = int(valtup[1]); exclrad = float(valtup[2])
    except:
        raise argparse.ArgumentTypeError("Exclusion zone must be a tripple like (int, int, float)")
    return (cx, cy, exclrad)

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
    parser.add_argument("--min-distance-from-tracking-centre",
                        type=int,
                        default=0,
                        help="Cutoff distance from phase centre in which no tags be raised or statistics be computed. "
                             "This can be used to effectively exclude the FWHM of an  parabolic reflector-based interferometer.")
    parser.add_argument("--add-custom-exclusion-zone",
                        type=exclusion_zone,
                        default=None,
                        nargs="+")
    args = parser.parse_args()
    import time
    tic = int(time.time())
    exclusion_zones = [exclz for exclz in args.add_custom_exclusion_zone] \
        if args.add_custom_exclusion_zone is not None else []
    tagged_regions = tag_regions(args.noise_map,
                                 regionsfn = args.ds9_reg_file,
                                 sigma = args.sigma,
                                 block_size = args.tile_size,
                                 hdu_id = 0,
                                 use_stokes = args.stokes,
                                 global_stat_percentile = args.global_rms_percentile,
                                 min_blocks_in_region = args.min_tiles_region,
                                 min_distance_from_centre = args.min_distance_from_tracking_centre,
                                 exclusion_zones=exclusion_zones)
    if args.input_lsm is not None:
        tag_lsm(args.input_lsm,
                args.noise_map,
                tagged_regions,
                hdu_id=0,
                regionsfn = args.ds9_tag_reg_file,
                taggedlsm_fn=args.input_lsm + ".de_tagged.lsm.html",
                de_tag=args.de_tag_name)
    toc = int(time.time())
    print>>log, "CATDagger ran successfully in {0:.0f}:{1:02.0f} minutes".format((toc - tic) // 60,
                                                                                 (toc - tic) % 60)
if __name__ == "__main__":
    main()
