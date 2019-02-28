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

import copy
import numpy as np
from astropy.io import fits
from astropy import wcs
import Tigger
from catdagger import logger
from catdagger.geometry import BoundingBox, BoundingConvexHull
log = logger.getLogger("lsm_tools")

def tag_lsm(lsm,
            stokes_cube,
            tagged_regions,
            hdu_id=0,
            regionsfn = "dE.srcs.reg",
            taggedlsm_fn="tagged.catalog.lsm.html",
            de_tag="dE",
            store_only_dEs=False):
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
                s.setTag(de_tag, True)
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
    if store_only_dEs:
        print>>log, "Removing direction independent components from catalog before writing LSM"
        ncomp_di_dies = len(mod.sources)
        mod.sources = filter(lambda s: de_tag in s.getTagNames(), mod.sources)
        ncomp_des = len(mod.sources)
        print>>log, "\t - Removed {0:d} direction independent sources from catalog".format(ncomp_di_dies - ncomp_des)
    print>>log, "Writing tagged LSM to {0:s}".format(taggedlsm_fn)
    mod.save(taggedlsm_fn)
    return mod.sources
