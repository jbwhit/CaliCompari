#!/usr/bin/env python
# encoding: utf-8
"""
helper.py

Created by Jonathan Whitmore on 2011-09-23.
Copyright (c) 2011 All rights reserved.
Email: JBWHIT@gmail.com
"""

import sys
import getopt
# import argparse
# parser.add_argument('foo', nargs='?', default=42) 
import os
import glob
import numpy as np

import scipy as sp
from scipy import arange, optimize, special, interpolate
from scipy.ndimage import *
import scipy.interpolate as si
import scipy.signal as ss 
from scipy.stats import linregress

import minuit as mi 
import matplotlib
import matplotlib.pylab as pl

import time
import datetime
from optparse import OptionParser
from ConfigParser import RawConfigParser
import random as ra
from scipy.optimize import leastsq
import scipy.constants
import pyfits as pf

c_light = sp.constants.c

# Directory structure:
# how to test for directory and how to create directory if not there.
# if not os.path.exists(dir):
    # os.makedirs(dir)
# Find out how to automatically load pylab for ipython and python startups.

help_message = '''
The help message goes here.
'''
# TODO 
# Integrate previous python code into this one which can take commandline 
# arguments
# fit a continuum 
# find the difference in wavelength calibration between two reference spectra
# allow for a slope to the fit
# write out results in fits files? automatically create QA (quality assurance) plots?
# Calculate total "information" in spectra
# Get the pixel numbers from the cleanup stage and store as an array -- if UVES

# ====================
# = Helper Functions =
# ====================

def weighted_std(measurements,errors):
  """Takes measurement array and errors and 
  turns it into the weighted standard deviation."""
  omega = (1.0 / errors) ** 2.0
  return np.sqrt(np.sum(omega * (measurements - np.average(measurements)) ** 2.0) / np.sum(omega))
  
def weighted_av(measurements,errors):
  """takes measurement array and errors and 
  turns it into the weighted average."""
  omega = (1.0 / errors) ** 2.0
  return np.sum(omega * measurements) / np.sum(omega)
  
def weighted_av_and_std(measurements,errors):
  """Prints both the weighted average and 
  the weighted standard deviation. """
  return weighted_av(measurements,errors), weighted_std(measurements,errors)

def normal_gaussian(elements, sigma):
  """docstring for norm_gaussian"""
  return ss.gaussian(elements, sigma) / sum(ss.gaussian(elements, sigma))

def normal_skew_gaussian(elements, sigma, skew):
  return skew_gaussian(elements, sigma, skew) / sum(skew_gaussian(elements, sigma, skew))

def skew_gaussian(elements, fsigma, skew):
  """docstring for skew_gaussian"""
  return whit_skew(arange(int(elements)),fsigma,int(elements/2),skew)

def whit_skew(x,o,c,a):
  """Returns an array 
  the gaussian of array(x), 
  scaling parameter o
  shift parameter c
  and skew parameter a ( a = 0 returns standard normal)
  """
  return (2.0 / o) * ( 1.0 / np.sqrt(2.0) ) * \
         np.exp( - ( ( (x - c) / o ) ** 2.) / 2.0 ) * \
         0.5 * (1.0 + sp.special.erf( a * ((x - c) / o) ) )        

  
if __name__ == "__main__":
  sys.exit(main())

