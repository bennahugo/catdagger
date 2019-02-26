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