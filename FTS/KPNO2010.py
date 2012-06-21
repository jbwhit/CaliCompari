#!/bin/python
# Copyright Jonathan Whitmore
# Distributed under the Boost Software License, Version 1.0.
# 
# Permission is hereby granted, free of charge, to any person or organization 
# obtaining a copy of the software and accompanying documentation covered by 
# this license (the "Software") to use, reproduce, display, distribute, 
# execute, and transmit the Software, and to prepare derivative works of the 
# Software, and to permit third-parties to whom the Software is furnished to 
# do so, all subject to the following:
# 
# The copyright notices in the Software and this entire statement, including 
# the above license grant, this restriction and the following disclaimer, 
# must be included in all copies of the Software, in whole or in part, and 
# all derivative works of the Software, unless such copies or derivative 
# works are solely in the form of machine-executable object code generated 
# by a source language processor.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND 
# NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE 
# DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, 
# WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# 
# Fit a continuum to the FTS spectrum KPNO2010.
# http://kurucz.harvard.edu/sun/irradiance2005/irradthu.dat
# Solar irradiance (Kurucz, R.L. 2005)   revised 02jan2010
# High resolution irradiance spectrum from 300 to 1000 nm.
# Thuillier et al (2004) calibration

import pylab
import numpy
import numpy as np
import scipy.interpolate as si

wav, flx = np.loadtxt('KPNO2010.txt', unpack=True, comments="#")
wav = 10. * wav # convert from nm to Angstroms)

binSize = 3000 # number of pixels of iodine per bin
smoothFactor = 0.0 # 0 forces spline through all points.
argMaxValues = [np.argmax(flx[i:i+binSize]) + i for i in \
                xrange(0, len(wav), binSize)]
spline = si.UnivariateSpline(wav[argMaxValues], flx[argMaxValues], \
                              s=smoothFactor)

# Plot the resulting continuum normalized spectra.
def plotSplineSmooth(binSize=binSize, smoothFactor=smoothFactor, \
                      wav=wav, flx=flx):
  """docstring for splinesmooth"""
  binSize = binSize # number of pixels of iodine
  argMaxValues = [np.argmax(flx[i:i+binSize]) + i for i in \
                  xrange(0, len(wav), binSize)] 
  spline = si.UnivariateSpline(wav[argMaxValues], flx[argMaxValues], \
                                s=smoothFactor)
  plot(wav[argMaxValues], flx[argMaxValues])
  plot(wav, flx/spline(wav))
  pass

np.savetxt('KPNO2010.contin.ascii', (wav, flx/spline(wav))) # Faster to load.
# np.savetxt('KPNO2010.contin-old.ascii', zip(wav, flx/spline(wav)))