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

import numpy as np
import scipy.stats as sstats
import scipy.signal as ssig
import scipy.spatial as spat
import Tigger
from astropy import wcs
from catdagger import logger
log = logger.getLogger("geometry")

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