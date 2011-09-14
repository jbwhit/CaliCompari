#!/usr/bin/env python
# encoding: utf-8
"""
2011-03-01-load.py

Modified by Jonathan Whitmore on 2011-03-01
"""
import sys
import os
import glob
import csv
import numpy as np
import scipy as sp
from scipy import arange, optimize, special, interpolate
from scipy.signal import *
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
from numpy import ma

# ============================
# = Make directory variables =
# ============================

# TODO Fix this crap

AsciiDir = "Ascii_Output/"
PdfDir = "PDF_Output/"
PsDir = "PS_Output/"

RawCopyDir = "Ascii_Output/Raw_Copy/"
RawMaskedDir = "Ascii_Output/Raw_Masked/"

OutputContinuumDir = "Ascii_Output/Continuum/"
OutputCalibrationDir = "Ascii_Output/Calibration/"
OutputResolutionDir = "Ascii_Output/Resolution/"

CalibrationChipDir = "Ascii_Output/Chip/"

PdfContinuumDir = "PDF_Output/Continuum/"
PsContinuumDir = "PS_Output/Continuum/"

PdfCalibrationDir = "PDF_Output/Calibration/"
PsCalibrationDir = "PS_Output/Calibration/"

PdfQADir = "PDF_Output/QA_Bins/"
PsQADir = "PS_Output/QA_Bins/"

PdfQAOrderResidualsDir = "QA_PDF/Order-Residuals/"
PdfQAOrderChiSquareDir = "QA_PDF/Order-Chi-Square/"
PdfQABinResidualsDir = "QA_PDF/Bins-Residuals/"
PdfQABinChiSquareDir = "QA_PDF/Bins-Chi-Square/"
PdfQABinOverlayDir = "QA_PDF/Bins-Overlay/"

PsQAOrderResidualsDir = "QA_PS/Order-Residuals/"
PsQAOrderChiSquareDir = "QA_PS/Order-Chi-Square/"
PsQABinResidualsDir = "QA_PS/Bins-Residuals/"
PsQABinChiSquareDir = "QA_PS/Bins-Chi-Square/"
PsQABinOverlayDir = "QA_PS/Bins-Overlay/"

OutputMonteCarloDir = "MonteCarlo/"

# ========================
# = Function Definitions =
# ========================


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

#
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
