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
from astropy import wcs
from catdagger import logger
log = logger.getLogger("filters")

class arealess():
    def __init__(self, min_area=0):
        self._min_area = min_area

    def __call__(self, reg):
        if reg.area < self._min_area:
            print>>log, "\t - Discarding region {0:s} because of its small size".format(reg.name)
            return True
        return False

class notin():
    def __init__(self, discard_list):
        self._dl = set(discard_list)

    def __call__(self, reg):
        return reg not in self._dl

class within_radius_from():
    """ centre within radius of (x, y) px"""
    def __init__(self, min_radius=0, cx=0, cy=0):
        self._min_radius = min_radius
        self._x = cx
        self._y = cy

    def __call__(self, reg):
        if np.sum((reg.centre - np.array([self._x, self._y]))**2) < self._min_radius**2:
            print>>log, "\t - Discarding region {0:s} because of its " \
                        "radial proximity to exclusion zone " \
                        "({1:.2f}, {2:.2f}, {3:.2f})".format(reg.name,
                                                            self._x,
                                                            self._y,
                                                            self._min_radius)
            return True
        return False