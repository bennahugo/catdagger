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

import os
import sys
import logging
from setuptools import setup, find_packages

requirements = ['numpy', 
                'astropy',
                'scipy',
                'astro-tigger-lsm']

setup(name='catdagger',
      version='0.1.0',
      description='An automatic differential gain catalog tagger',
      url='https://github.com/bennahugo/catdagger',
      classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Astronomy"],
      author='Benjamin Hugo',
      author_email='bhugo@ska.ac.za',
      license='GNU GPL v3',
      packages=['catdagger'],
      scripts=['bin/dagger'],
      install_requires=requirements,
      include_package_data=True,
      zip_safe=False,
)


