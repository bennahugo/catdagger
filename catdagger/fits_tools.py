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
from astropy.io import fits
from astropy import wcs
from catdagger import logger
log = logger.getLogger("FITS_tools")

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

def getcrpix(fn, hdu_id, use_stokes="I"):
    stokes_cube = fn
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
    return hdr["CRPIX{0:d}".format(types["RA---SIN"])], \
           hdr["CRPIX{0:d}".format(types["DEC--SIN"])]

def read_stokes_slice(fn,
                      hdu_id = 0, 
                      use_stokes="I",
                      average_channels=True):
    stokes_cube = fn
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
    if average_channels:
        print>>log, "Collapsing axis: {0:d} (FREQ)".format(types["FREQ"])
        band_avg = np.mean(sel_stokes, axis=chan_axis)
        return w, band_avg 
    else:
        return w, sel_stokes

def blank_components(fn, list_ra_dec_rad, hdu_id = 0, use_stokes="I"):
    w, data = read_stokes_slice(fn, hdu_id, use_stokes)
    pass
