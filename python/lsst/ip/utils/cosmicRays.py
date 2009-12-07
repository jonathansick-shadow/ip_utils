#!/usr/bin/env python
"""
Utilities to handle CRs
"""
import math
import os
import pdb                          # we may want to say pdb.set_trace()
import unittest

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg

def findCosmicRays(exposure, crRejectPolicy, defaultFwhm, keepCRs):
    """defaultFwhm is in arcsec"""

    mi = exposure.getMaskedImage()
    wcs = exposure.getWcs()

    scale = math.sqrt(wcs.pixArea(afwImage.PointD(mi.getWidth()/2, mi.getHeight()/2)))*3600 # arcsec/pixel
    defaultFwhm /= scale            # convert to pixels

    psf = measAlg.createPSF('DoubleGaussian', 0, 0, defaultFwhm/(2*math.sqrt(2*math.log(2))))

    bg = afwMath.makeStatistics(mi, afwMath.MEDIAN).getValue()
    crs = measAlg.findCosmicRays(mi, psf, bg, crRejectPolicy, keepCRs)

    return crs
