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
import scipy.signal as ssig
from catdagger.gauss2 import twodgaussian
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

def get_fitted_beam(fn, hdu_id):
    print>>log, "Finding fitted CLEAN beam parameters in {0:s}".format(fn)
    stokes_cube = fn
    with fits.open(stokes_cube) as img:
        hdr = img[hdu_id].header
        w = wcs.WCS(hdr)
    if not all([k in hdr for k in ["BPA", "BMAJ", "BMIN"]]):
        raise KeyError("Fitted clean beam parameters is not in FITS file")
    return hdr["BMIN"], hdr["BMAJ"], hdr["BPA"]

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
        return w, hdr, band_avg 
    else:
        return w, hdr, sel_stokes

def save_stokes_slice(fn,
                      cube_slice,
                      hdu_id = 0, 
                      use_stokes="I",
                      backup=True):
    stokes_cube = fn
    with fits.open(stokes_cube) as img:
        cube = img[hdu_id].data
        hdr = img[hdu_id].header
        w = wcs.WCS(hdr)
        backup and img.writeto(stokes_cube + ".orig.fits", overwrite=True)
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
    stokes_slice_indx = tuple([slice(None) if k != "STOKES" else sel_stokes for k in sorted(types.keys(), key=lambda k: types[k], reverse=True)])
    cube[stokes_slice_indx] = cube_slice
    chan_axis = hdr["NAXIS"] - types["FREQ"] if types["FREQ"] > types["STOKES"] else hdr["NAXIS"] - types["FREQ"] - 1
    print>>log, "Saving model FITS back to disk: {0:s}".format(stokes_cube)
    with fits.open(stokes_cube) as img:
        img[hdu_id].data = cube
        img.writeto(stokes_cube, overwrite=True)


def blank_components(fn, rmsmap, psf_image, list_src, hdu_id = 0, use_stokes="I"):

    w, hdr, data = read_stokes_slice(fn, hdu_id, average_channels=False)
    types = {hdr["CTYPE{0:d}".format(ax + 1)]: (ax + 1) for ax in range(hdr["NAXIS"])}
    cdelt = float(max(np.abs(hdr["CDELT{0:d}".format(types["RA---SIN"])]), \
                      np.abs(hdr["CDELT{0:d}".format(types["DEC--SIN"])])))
    BMIN, BMAJ, BPA = get_fitted_beam(psf_image, hdu_id) 
    BMAJ=int(BMAJ / cdelt) # in pixels
    BMIN=int(BMIN / cdelt) # in pixels
    print>>log, "Blanking the following positions with fitted resolution:"
    for s in list_src:
        ra = s.pos.ra
        dec = s.pos.dec
        ex = int(s.shape.ex / np.rad2deg(cdelt)) if s.shape is not None and \
                        hasattr(s.shape, "typecode") and \
                        s.shape.typecode == "Gau" \
            else 0
        ey = int(s.shape.ex / np.rad2deg(cdelt)) if s.shape is not None and \
                        hasattr(s.shape, "typecode") and \
                        s.shape.typecode == "Gau" \
            else 0
        epa = s.shape.pa if s.shape is not None and \
                         hasattr(s.shape, "typecode") and \
                         s.shape.typecode == "Gau" \
            else 0
        emaj = max(ex, ey)
        emin = min(ex, ey)

        # Build a mask of source convolved with resolution
        if emaj != 0 and emin != 0:
            wnd_size = (2 * 10 * (max(emaj, BMAJ) + 1))
        else:
            wnd_size = (2 * 10 * BMAJ + 1)
        x, y = np.meshgrid(np.arange(-wnd_size//2, wnd_size//2),
                           np.arange(-wnd_size//2, wnd_size//2))
        wnd = twodgaussian([1.0, 0, 0, BMAJ, BMIN, BPA],
                           circle=0, rotate=1, vheight=0)(x, y)
        if emaj != 0 and emin != 0:
            src = twodgaussian([1.0, 0, 0, emaj, emin, np.rad2deg(epa)],
                               circle=0, rotate=1, vheight=0)(x, y)
        else:
            src = np.zeros_like(wnd)
            src[wnd_size//2, wnd_size//2] = 1
        conv_wnd = ssig.convolve(src, wnd, mode="same")
        # place the mask over the image portion and set everything in it to 0
        wnd_mask = np.zeros_like(wnd, dtype=np.bool)
        wnd_mask[conv_wnd >= 0.5 * np.max(wnd)] = True
        wnd_mask_cube = np.empty((data.shape[0], wnd.shape[0], wnd.shape[1]), dtype=np.bool)
        wnd_mask_cube[:, :, :] = wnd_mask[None, :, :]

        y, x, _, _ = w.all_world2pix([[np.rad2deg(ra), np.rad2deg(dec), 0, 0]], 1)[0]
        x = int(x)
        y = int(y)
        if np.isnan(x) or np.isnan(y): continue
        in_image = x >= 0 and x < data.shape[1] and \
                   y >= 0 and y < data.shape[2]
        wnd_size_half = wnd_size // 2
        xmin = np.clip(x - wnd_size_half - 1, 0, data.shape[1] - 1)
        xmax = np.clip(x + wnd_size_half, 0, data.shape[1] - 1)
        ymin = np.clip(y - wnd_size_half - 1, 0, data.shape[2] - 1)
        ymax = np.clip(y + wnd_size_half, 0, data.shape[2] - 1)
        wnd_cut_min_x = np.clip(0 - (x - wnd_size_half - 1), 0, x + wnd_size_half)
        wnd_cut_max_x = wnd_size - np.clip(x + wnd_size_half - (data.shape[1] - 1), 0, x + wnd_size_half)
        wnd_cut_min_y = np.clip(0 - (y - wnd_size_half - 1), 0, y + wnd_size_half)
        wnd_cut_max_y = wnd_size - np.clip(y + wnd_size_half - (data.shape[1] - 1), 0, y + wnd_size_half)
        assert xmax - xmin == wnd_cut_max_x - wnd_cut_min_x
        assert ymax - ymin == wnd_cut_max_y - wnd_cut_min_y
        if in_image:
            mod_flux = np.sum(data[:, xmin:xmax, ymin:ymax][wnd_mask_cube[:,
                                                                          wnd_cut_min_x:wnd_cut_max_x, 
                                                                          wnd_cut_min_y:wnd_cut_max_y]])

            print>>log, "\t - {0:d}, {1:d} with {2:0.2f} integrated flux (mJy) within resolution".format(
                x, y, mod_flux * 1.0e3 / data.shape[0])
            data[:, xmin:xmax, ymin:ymax][wnd_mask_cube[:,
                                                        wnd_cut_min_x:wnd_cut_max_x, 
                                                        wnd_cut_min_y:wnd_cut_max_y]] = 0.0
    save_stokes_slice(fn, data, hdu_id, use_stokes)
