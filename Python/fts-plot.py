#!/usr/bin/env python
# encoding: utf-8
"""
fts-plot.py

Created by Jonathan Whitmore on 2010-05-27.
Modified: 
"""

version = 1.0

import csv
import numpy as np
import scipy as sp
from scipy import arange, optimize, special, interpolate
import scipy.interpolate as si
import scipy.signal as ss 
from scipy.ndimage import *
import minuit as mi 
import matplotlib.pylab as pl
import sys
import os
import glob
import datetime
from optparse import OptionParser
from ConfigParser import RawConfigParser
from scipy.stats import linregress



wav, flx, err, b1, b2, b3 = np.loadtxt('Test.FTS.01.l.01.ascii',unpack=True)


print "Last data point weird -- excluding it"
wav = wav[:-2]
flx = flx[:-2]
err = err[:-2]

print "Maximum deviation: ", np.max(flx) - np.min(flx), "m/s\n\n"

pl.plot(wav,flx)
pl.errorbar(wav,flx,yerr=err)

print "Regression: ", linregress(wav,flx)[0]



# Run from uves/HR1996 