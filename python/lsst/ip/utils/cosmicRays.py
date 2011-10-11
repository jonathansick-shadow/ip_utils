#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

"""
Utilities to handle CRs
"""
import math
import os
import unittest

import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg

def findCosmicRays(exposure, crRejectPolicy, defaultFwhm, keepCRs):
    """defaultFwhm is in arcsec"""

    mi = exposure.getMaskedImage()
    wcs = exposure.getWcs()

    scale = wcs.pixelScale().asArcseconds()
    defaultFwhm /= scale            # convert to pixels
    ksize = 4*int(defaultFwhm) + 1

    psf = afwDetection.createPsf('DoubleGaussian', ksize, ksize, defaultFwhm/(2*math.sqrt(2*math.log(2))))

    bg = afwMath.makeStatistics(mi, afwMath.MEDIAN).getValue()
    crs = measAlg.findCosmicRays(mi, psf, bg, crRejectPolicy, keepCRs)

    return crs
