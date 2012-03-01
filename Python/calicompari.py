#!/usr/bin/env python
# encoding: utf-8
"""
calicomparison.py

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

# TODO uncomment and provide useful stuff
class Usage(Exception):
  def __init__(self, msg):
    self.msg = msg

# ===========
# = Cleanup =
# ===========
def cleanup(inputarray):
  """Cleanup the input spectrum"""
  if 'verbose' in globals():
    print "Beginning cleanup of data...", datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
  try:
    wav, flx, err, trash = inputarray
  except:
    print "Don't have 4 values, huh?"
    wav, flx, err = inputarray
  # wav, flx, err, trash = inputarray.data
  print "Max Flux", np.max(flx)
  print "Min Flux", np.min(flx)
  positiveerror = np.where(err > error_lower)
  highsn = np.where(flx[positiveerror]/err[positiveerror] > sn_lower)
  # Remove odd discrepancies near the edge
  cleanindex = highsn[0][edgebuffer:-edgebuffer]

  cleanwav = wav[cleanindex]
  cleanflx = flx[cleanindex]
  cleanerr = err[cleanindex]

  # Sometimes this makes things too small, so punt early 
  if len(cleanwav) < 200:
    print "Too small -- something's wrong..."
  
  dummywav = cleanwav
  for x in zip(begin_kill_array, end_kill_array):
    dummywav = [y for y in dummywav if not (y > x[0] and y < x[1])]
  # finalindex = [np.argwhere(cleanwav == x)[0][0] for x in dummywav]
  finalindex = [np.argwhere(wav == x)[0][0] for x in dummywav]
  
  if 'verbose' in globals():
    print "Finished cleaning up. ", datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
  print "Max Flux", np.max(flx[finalindex])
  print "Min Flux", np.min(flx[finalindex])  
  return wav[finalindex], flx[finalindex], err[finalindex]

# =====================
# = Continuum Fitting =
# =====================
def continuumfit(wav, flx, err):
  """docstring for continuum"""
  if 'verbose' in globals():
    print "Beginning continuum fitting of data...", datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
  x0 = wav
  y0 = flx
  
  beginxcontin = wav[0] - 5.5
  endxcontin = wav[-1] + 5.5
  numberofknots = 10
  xeval = np.arange(beginxcontin,endxcontin, \
                    float((endxcontin - beginxcontin)/numberofknots))
  
  padding = np.array([0,0,0,0])
  
  part1 = np.insert(xeval,0,beginxcontin)
  part1 = np.insert(part1,0,beginxcontin)
  part1 = np.insert(part1,0,beginxcontin)
  part1 = np.append(part1,endxcontin)
  part1 = np.append(part1,endxcontin)
  part1 = np.append(part1,endxcontin)
  part1 = np.append(part1,endxcontin)
  
  part2 = np.append(-1*np.ones(len(part1) - len(padding)),padding)
  
  splineorder = 3
  
  fitfunc = lambda part2: interpolate.splev(x0,(part1,part2,splineorder))
  errfunc = lambda part2: fitfunc(part2) - y0
  p2, success = optimize.leastsq(errfunc, part2[:])
  
  continuum_fit_y = flx/interpolate.splev(x0,(part1,p2,splineorder))
  continuum_fit_err = err/interpolate.splev(x0,(part1,p2,splineorder))
  continuum_points = interpolate.splev(x0,(part1,p2,splineorder))
  
  if 'verbose' in globals():
    print "Finished continuum fitting.", datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
  return continuum_points
  

# =======================
# = Calibration Fitting =
# =======================
minPixels = 100.
minReferenceOverlap = 80.
maxReferenceOverlap = 100.
# ==============================
# = Configuration file parsing =
# ==============================
filename='config'
config = RawConfigParser()
config.read(filename)
error_lower = float(config.get('Parameters','error_lower'))
# error_lower = 0.0
sn_lower = float(config.get('Parameters','sn_lower'))
# TODO Make this s_upper optional
s_upper = float(config.get('Parameters','s_upper'))
edgebuffer = int(config.get('Parameters','edgebuffer'))
telescope = config.get('Information','telescope')
astro_object = config.get('Information','astro_object')
spectrograph = config.get('Information','spectrograph')

bin_size = float(config.get('Parameters','bin_size'))  # 300.0 # in km/s
step = float(config.get('Parameters','step_size')) #  50.0 # km/s
gauss_elements = float(config.get('Parameters','gauss_elements'))
sigma = float(config.get('Parameters','sigma'))
i2exposures = config.get('Information','Exposures').strip().splitlines()
chips = config.get('Information','Multi_Chips').strip().splitlines()
destroy = config.get('Parameters','Wavelength Remove').strip().splitlines()
begin_kill_array = []
end_kill_array = []
for i in range(len(destroy)):
  t1, t2 = destroy[i].split()
  begin_kill_array.append(float(t1))
  end_kill_array.append(float(t2))

begin_kill_array = np.array(begin_kill_array)
end_kill_array = np.array(end_kill_array)

if spectrograph=="UVES": 
  print "Instrument: ", spectrograph
  FTSFile = config.get('Information','path_to_UVES_FTS')
if spectrograph=="HIRES":
  print "Instrument: ", spectrograph
  FTSFile = config.get('Information','path_to_HIRES_FTS')


def calibration(wav, flx, err, con):
  """Does the hard work of calibrating things."""
  def initial_shift(fmultiple,fshift,fsigma):
    better_flx = starting_flx * fmultiple
    better_kernel = normal_gaussian(gauss_elements,fsigma)
    better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                           slice_iof, mode='same'))
    better_y2 = si.splev(wav+fshift,better_tck)
    return sum((better_y2 - better_flx) ** 2 / \
               (fmultiple * err / con) ** 2 )

  def fullordershift(fmultiple,fshift,fsigma):
    better_flx = starting_flx * fmultiple
    better_kernel = normal_gaussian(gauss_elements,fsigma)
    better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                           slice_iof, mode='same'))
    better_y2 = si.splev(wav+fshift,better_tck)
    return (better_y2 - better_flx) ** 2 / \
               (fmultiple * err / con) ** 2 

  # def continuum_shift(fmultiple,fshift,fsigma,cmultiple):
  #   better_flx = flx / (con + con * cmultiple) * fmultiple
  #   better_kernel = normal_gaussian(gauss_elements,fsigma)
  #   better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
  #                          slice_iof, mode='same'))
  #   better_y2 = si.splev(wav+fshift,better_tck)
  #   return sum((better_y2 - better_flx) ** 2 / \
  #              (fmultiple * err / (con + con * cmultiple)) ** 2 )
  # 
  # def continuumshiftarray(fmultiple,fshift,fsigma,cmultiple):
  #   better_flx = flx / (con + con * cmultiple) * fmultiple
  #   better_kernel = normal_gaussian(gauss_elements,fsigma)
  #   better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
  #                          slice_iof, mode='same'))
  #   better_y2 = si.splev(wav+fshift,better_tck)
  #   return (better_y2 - better_flx) ** 2 / \
  #              (fmultiple * err / (con + con * cmultiple)) ** 2 
  # 
  # def second_shift(fmultiple,fshift,fsigma,fskew):
  #   better_flx = starting_flx * fmultiple
  #   better_kernel = normal_skew_gaussian(gauss_elements,fsigma,fskew)
  #   better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
  #                          slice_iof, mode='same'))
  #   better_y2 = si.splev(wav+fshift,better_tck)
  #   return sum((better_y2 - better_flx) ** 2 / \
  #              (fmultiple * err / con) ** 2 )
  # 
  # def skewshiftarray(fmultiple,fshift,fsigma,fskew):
  #   better_flx = starting_flx * fmultiple
  #   better_kernel = normal_skew_gaussian(gauss_elements,fsigma,fskew)
  #   better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
  #                          slice_iof, mode='same'))
  #   better_y2 = si.splev(wav+fshift,better_tck)
  #   return (better_y2 - better_flx) ** 2 / \
  #              (fmultiple * err / con) ** 2 
  # 
  def shiftperbin(i,fmultiple,fshift,fsigma):
    better_flx = starting_flx[whit_index[int(i)]] * fmultiple
    better_kernel = normal_gaussian(gauss_elements,fsigma)
    better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                           slice_iof, mode='same'))
    better_y2 = si.splev(wav[whit_index[int(i)]]+fshift, better_tck)
    return sum((better_y2 - better_flx) ** 2 / \
                (fmultiple * starting_err[whit_index[int(i)]] ) ** 2)

  # def shiftperbinarray(i,fmultiple,fshift,fsigma):
  #   better_flx = starting_flx[whit_index[int(i)]] * fmultiple
  #   better_kernel = normal_gaussian(gauss_elements,fsigma)
  #   better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
  #                          slice_iof, mode='same'))
  #   better_y2 = si.splev(wav[whit_index[int(i)]]+fshift, better_tck)
  #   return (better_y2 - better_flx) ** 2 / \
  #               (fmultiple * err[whit_index[int(i)]] / con[whit_index[int(i)]]) ** 2
  # 
  # def plotresidualshiftperbin(fshift,i,fmultiple,fsigma):
  #   better_flx = starting_flx[whit_index[int(i)]] * fmultiple
  #   better_kernel = normal_gaussian(gauss_elements,fsigma)
  #   better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
  #                          slice_iof, mode='same'))
  #   better_y2 = si.splev(wav[whit_index[int(i)]]+fshift, better_tck)
  #   pl.plot((better_y2 - better_flx) )
  #   return
  # 
  # def plotsidebyside(i):
  #   pl.plot(wav[fitbinindex],better_y2,wav[fitbinindex],flx[fitbinindex] * fitbinmultiple)
  # 
  # def plotcontinuum_shift(fmultiple,fshift,fsigma,cmultiple):
  #   better_flx = flx / (con + con * cmultiple) * fmultiple
  #   better_kernel = normal_gaussian(gauss_elements,fsigma)
  #   better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
  #                          slice_iof, mode='same'))
  #   better_y2 = si.splev(wav+fshift,better_tck)
  #   return sum((better_y2 - better_flx) ** 2 / \
  #              (fmultiple * err / (con + con * cmultiple)) ** 2 )
  # 
  # def skewshiftperbin(fshift,i,fmultiple,fsigma,fskew):
  #   better_flx = starting_flx[whit_index[int(i)]] * fmultiple
  #   better_kernel = normal_skew_gaussian(gauss_elements,fsigma,fskew)
  #   better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
  #                          slice_iof, mode='same'))
  #   better_y2 = si.splev(wav[whit_index[int(i)]]+fshift, better_tck)
  #   return sum((better_y2 - better_flx) ** 2 / \
  #               (fmultiple * starting_err[whit_index[int(i)]]) ** 2)
  
  if 'verbose' in globals():
    print "Beginning calibration analysis of data.", datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
  if len(wav) < minPixels:
    print "Less than ", minPixels, " pixels, skipping..."
    # continue
  # ================================================
  # = Logic Test to see if file overlaps FTS data =
  # ================================================
  iow, iof = np.loadtxt(FTSFile, unpack=True) # Load in the reference spectrum
  
  if(wav[0] < iow[0] + minReferenceOverlap or wav[-1] > iow[-1] - maxReferenceOverlap):
    print "Outside of overlap"
    return

  # print today
  # =================================
  # = Slice FTS to manageable size  =
  # =================================
  slice_iow = iow[np.where(iow > wav[0] - 10)]
  slice_iow = slice_iow[np.where(slice_iow < wav[-1] + 10)]
  slice_iof = iof[np.where(iow > wav[0] - 10)]
  slice_iof = slice_iof[np.where(slice_iow < wav[-1] + 10)]

  # Divide continuum out from input file
  starting_flx = flx / con
  starting_err = err / con
  # pl.plot(wav, starting_flx)
  # pl.savefig("first.pdf")
  # pl.close()
  # ======================================
  # = Initial Shift; Kernel Size finding =
  # ======================================
  # Should I make this more robust? -- maybe have a flag for an uncertain offset?
  # print "Initial shift is: ", initialshift

  # TODO comment out using continuum; use slope on FTS instead.
  m = mi.Minuit(initial_shift,\
                fmultiple=0.82,\
                # fix_fmultiple=True,\
                fshift=0.05,\
                fsigma=sigma,\
                # fix_fsigma=True,\
                strategy=2)
  
  # m.printMode=1
  print "Scanning..."
  m.scan(("fshift",20,-1.5,1.5))
  print "done."
  try: 
    m.migrad()
  except:
    print "Serious problem with the entire order."
  if m.values["fsigma"] < 0.:
    m = mi.Minuit(initial_shift,\
                  fmultiple=m.values["fmultiple"],\
                  fshift=m.values["fshift"],\
                  fsigma=-m.values["fsigma"],\
                  strategy=2)
    m.migrad()
  
  #m.minos() 
  order_velocity = c_light / wav[len(wav)/2]
  # TODO add in the ability for a slope in the fit -- not just an overall multiple
  # slopeArray = 
  # Review finished to here... 
  fitordermultiple = m.values["fmultiple"]
  fitordermultiple_err = m.errors["fmultiple"]
  
  fitordershift = m.values["fshift"]
  fitordershift_err = m.errors["fshift"]
  
  fitordersigma = m.values["fsigma"]
  fitordersigma_err = m.errors["fsigma"]
  
  print "\nShift: ", round(fitordershift * order_velocity,2), "m/s", \
          round(fitordershift_err * order_velocity,4), "m/s"
  print "Sigma: ", round(fitordersigma,2), round(fitordersigma_err,4)
  
  orderchisquarearray = fullordershift(fitordermultiple,fitordershift,fitordersigma)
  orderchisquare = np.sum(orderchisquarearray)
  DOF = 3
  print "Chi-Square: ", round(orderchisquare,2), "Chi-Sq/DOF: ", round(orderchisquare / (len(orderchisquarearray) - DOF),4)
  # titlestr = nodotnaming + " Total Chi-Square: " + str(round(orderchisquare,2)) + \
  #     " Chi-Sq/DOF: " + str(round(orderchisquare / (len(orderchisquarearray) - DOF),2))
  # pl.title(titlestr)        
  # pl.plot(wav,fullordershift(fitordermultiple,fitordershift,fitordersigma))
  # pl.xlabel("Wavelength in Angstroms")
  # pl.ylabel("Chi-Square")
  # pl.savefig(PdfQAOrderChiSquareDir + astro_object + naming + ".pdf")
  # pl.savefig(PsQAOrderChiSquareDir + astro_object + naming + ".eps")
  # pl.close()
  # Warning if constraining fitting function
  # Does this make sense?
  if gauss_elements/fitordersigma < 7.5:
    print "Warning: You should use a larger gaussian window\nYour current array is: ", gauss_elements
    print "Your current sigma is: ", fitordersigma, \
        " and elements/sigma: ", gauss_elements/fitordersigma
  

  # ========================================
  # = Create the bins and structure needed =
  # ========================================
  # Find minimum spacing between data points: use this to figure out how to oversample
  # np.min([wav[i+1] - wav[i] for i in range(len(wav)-1)])
  
  
  whit_index = []
  whit_step = wav[0] * step * 1000.0 / c_light
  for i in np.arange(wav[0],wav[-1], whit_step ):
    temp = i + i * bin_size * 1000.0 / c_light
    if len(np.where(wav[np.where(wav < temp)] > i )[0]) > 100 :
      whit_index.append(np.where(wav[np.where(wav < temp)[0]] > i ))
  
  wav_bin_av = []
  calib_fit = []
  calib_err = []
  resolution_fit = []
  resolution_err = []
  bin_call = []
  pix_average = []

  for i in np.arange(len(whit_index)):
    if(len(wav[whit_index[i]]) < 100):
      continue
    bin_call = mi.Minuit(shiftperbin,\
      i=i,\
      fix_i=True,\
      fmultiple=fitordermultiple,\
      # fix_fmultiple=True,\
      fshift=fitordershift,\
      fsigma=fitordersigma,\
      # fix_fsigma=True,\ #test commented 
      strategy=2)
    try:
      bin_call.migrad()
    except:
      print "Bin Broke! Moving on. Wavelength region: ", wav[whit_index[i]][0], wav[whit_index[i]][-1]
      continue
      
    fitbinmultiple = bin_call.values["fmultiple"]
    fitbinshift =    bin_call.values["fshift"]
    fitbinsigma =    bin_call.values["fsigma"]
    fitbinindex =    whit_index[i]
    
    lamb = np.average(wav[fitbinindex])
    midpoint = np.argmin( np.abs( [ slice_iow - lamb ] ) )
    # if len(wav) == len(pix):
    #   temp_pix = np.average(pix[fitbinindex])
    FWHM_prefactor = 2.0 * np.sqrt( 2.0 * np.log( 2.0 ) )
    
    # Full difference across considered lambdas
    full_lambda = slice_iow[midpoint + gauss_elements /2 ] - slice_iow[midpoint - gauss_elements/2]
    
    sigma_lambda = full_lambda * (fitbinsigma / gauss_elements)
    sigma_err_plus = full_lambda * (fitbinsigma + bin_call.errors["fsigma"]) / gauss_elements
    sigma_err_minus = full_lambda * (fitbinsigma - bin_call.errors["fsigma"]) / gauss_elements
    FWHM = FWHM_prefactor * sigma_lambda 
    FWHM_err1 = FWHM_prefactor * sigma_err_plus
    FWHM_err2 = FWHM_prefactor * sigma_err_minus

    wav_bin_av.append(lamb)
    calib_fit.append(fitbinshift)    
    calib_err.append(bin_call.errors["fshift"])
    resolution_fit.append( lamb / FWHM )
    resolution_err.append( np.abs( lamb/FWHM_err1 - lamb/FWHM_err2) /2.0 ) 
    # if len(wav) == len(pix):
    #   pix_average.append(temp_pix)
    # if i == len(whit_index)/2:
    #   binchisquarearray = shiftperbinarray(i,fitbinmultiple,fitbinshift,fitbinsigma)
    #   binchisquare = np.sum(binchisquarearray)
    #   plot_better_g = ss.gaussian(gauss_elements,fitbinsigma)
    #   better_kernel = plot_better_g / sum(plot_better_g)
    #   better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
    #                          slice_iof, mode='same'))
    #   better_y2 = si.splev(wav[whit_index[int(i)]]+fitbinshift,better_tck)
    #   
    #   pl.plot(wav[fitbinindex],flx[fitbinindex] * fitbinmultiple / con[fitbinindex],label=spectrograph)
    #   tempx = pl.xlim()
    #   pl.plot(wav[fitbinindex],better_y2,label="FTS")
    #   pl.xlim(tempx)
    #   pl.legend()
    #   pl.ylabel("Normalized Flux")
    #   pl.xlabel("Wavelength in Angstroms")
    #   titlestr = nodotnaming + " Mid-order Bin Fit-check"
    #   pl.title(titlestr)
    #   PdfQABinOverlayFile = PdfQABinOverlayDir + astro_object + naming + ".pdf"
    #   PsQABinOverlayFile = PsQABinOverlayDir + astro_object + naming + ".eps"
    #   pl.savefig(PdfQABinOverlayFile)
    #   pl.savefig(PsQABinOverlayFile)
    #   pl.close()
      
  
  wav_bin_av = np.array(wav_bin_av) 
  calib_fit = np.array(calib_fit)
  calib_err = np.array(calib_err)
  resolution_fit = np.array(resolution_fit)
  resolution_err = np.array(resolution_err)
  pix_average = np.array(pix_average)
  
  
  # OutputCalibrationFile = OutputCalibrationDir + astro_object + naming + ".ascii"
  # OutputResolutionFile = OutputResolutionDir + astro_object + naming + ".ascii" 

  # print >>CalibrationLogFileHandle, CalibrationDataFile, wav[len(wav)/2], fitordershift * order_velocity, \
          # fitordershift_err * order_velocity
  
  velocity_fit = c_light * calib_fit / wav_bin_av
  velocity_err = c_light * calib_err / wav_bin_av
  

  # print >>CalibrationSummaryFileHandle, astro_object + naming, weighted_av(velocity_fit,velocity_err), \
  #     weighted_std(velocity_fit,velocity_err), \
  #     fitordermultiple,  fitordermultiple_err, \
  #     fitordershift, fitordershift_err, fitordersigma, fitordersigma_err, order_velocity, \
  #     FWHM, (FWHM_err1 - FWHM_err2) / 2.0
  
  # # moved to calibration-plot.py
  # if(len(calib_fit) > 0):
  #   PdfCalibrationFile = PdfCalibrationDir + astro_object + naming + ".pdf" 
  #   PsCalibrationFile = PsCalibrationDir + astro_object + naming + ".eps" 
  # 
  #   titlestr = "Exposure: " + secondsplit[1] + " Chip: " + secondsplit[2] + " Order: " + secondsplit[3]
  # 
  #   pl.subplot(311)
  #   pl.ylabel("Flux")
  #   pl.title(titlestr)
  #   pl.plot(wav,flx * fitordermultiple,wav,err * fitordermultiple)
  #   temp = pl.xlim()
  #   pl.subplot(312)
  #   pl.ylabel("v_cal (m/s)")
  #   pl.errorbar(wav_bin_av,velocity_fit,yerr=velocity_err)
  #   pl.plot(wav_bin_av,np.zeros(len(wav_bin_av)))
  #   pl.xlim(temp)
  #   pl.subplot(313)
  #   pl.ylabel("R value")
  #   pl.errorbar(wav_bin_av,resolution_fit,yerr=resolution_err)
  #   pl.xlim(temp)
  #   pl.savefig(PdfCalibrationFile)
  #   pl.savefig(PsCalibrationFile)
  #   pl.close()
  # 
  # OutputCalibrationFileHandle = open('output.ascii', 'a')
  print "Weighted Average: ", round(weighted_av(velocity_fit,velocity_err),2), "m/s"
  print "Weighted std: ", round(weighted_std(velocity_fit,velocity_err),2), "m/s"
  
  # for j in range(len(wav_bin_av)):
  #   # if len(wav) == len(pix):
  #   #   print >>OutputCalibrationFileHandle, wav_bin_av[j], velocity_fit[j], velocity_err[j], \
  #   #           resolution_fit[j], resolution_err[j], \
  #   #           pix_average[j]
  #   #   continue
  #   print >>OutputCalibrationFileHandle, wav_bin_av[j], velocity_fit[j], velocity_err[j], \
  #         resolution_fit[j], resolution_err[j]
  # 
  # OutputCalibrationFileHandle.close()
  # OutputResolutionFileHandle = open(OutputResolutionFile,'w')
  # for j in range(len(wav_bin_av)):
      # print >>OutputResolutionFileHandle, wav_bin_av[j], resolution_fit[j], resolution_err[j] 
  # 
  # OutputResolutionFileHandle.close()
  
  if 'verbose' in globals():
    print "Finished calibration run. ", datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
  
  return wav_bin_av, velocity_fit, velocity_err, resolution_fit, resolution_err

# =======================
# = Main Program begins =
# =======================
def main(argv=None):
  # initialshift = 0.01 # default shift
  if argv is None:
    argv = sys.argv
  try:
    try:
      opts, args = getopt.getopt(argv[1:], "hc:r:s:o:v", ["help", "compare=", "reference=", "shift=" "output="])
    except getopt.error, msg:
      raise Usage(msg)
    # TODO flag for reference file; 
    # TODO flag for plotting QA stuff
    # option processing
    for option, value in opts:
      if option == "-v":
        global verbose
        verbose = True
        print "Verbose turned on."
      if option in ("-c", "--compare"):
        inputfitsfile = value
      if option in ("-r", "--reference"):
        referencefile = value
      if option in ("-s", "--shift"):
        initialshift = value
      if option in ("-h", "--help"):
        raise Usage(help_message)
      if option in ("-o", "--output"):
        output = value
  
  except Usage, err:
    print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
    print >> sys.stderr, "\t for help use --help"
    return 2
  

  # ==================================
  # = Print input values to terminal =
  # ==================================
  if 'verbose' in globals():
    print "Chips: ", chips
    print "Telescope: ", telescope, " Object: ", astro_object
    print "Error Lower Bound: ", error_lower
    print "S/N Lower Bound: ", sn_lower
    print "Edgebuffer: ", edgebuffer
    print "FTS File: ", FTSFile
    print "Bin size: ", bin_size, "km/s", "Step size: ", step, "km/s"
    print "# of Elements in gaussian array: ", gauss_elements, "Initial sigma: ", sigma
  
  
  # core_test = 0
  # pf.writeto(file, clobber=True)
  today = datetime.date.today()
  calibrationfitsfile = str("cali_" + inputfitsfile)
  if inputfitsfile.split('.')[-1] == "dat":
    x = np.loadtxt(inputfitsfile,usecols=(0,1,2,3),unpack=True)
    wav, flx, err = cleanup(x)
    if not len(wav) == len(flx) == len(err):
      print "Something's wrong with the cleanup! Breaking!"
    con = continuumfit(wav,flx,err)
    avwav, cal, calerr, res, reserr = calibration(wav, flx, err, con)
    pf.append(calibrationfitsfile, array([avwav, cal, calerr, res, reserr]), verify=False)
  if inputfitsfile.split('.')[-1] == "fits":
    for x in pf.open(inputfitsfile):
      if type(x) == pf.core.ImageHDU:
        try:
          wav, flx, err = cleanup(x.data)
        except: 
          print "That entire order broke... Moving on to the next one!"
          continue
        if not len(wav) == len(flx) == len(err):
          print "Something's wrong with the cleanup! Breaking!"
          break
        # print type(wav)
        try: 
          con = continuumfit(wav,flx,err)
        except: print "whatever1"
        try:
          avwav, cal, calerr, res, reserr = calibration(wav, flx, err, con)
        except Exception, err: 
          print "whatever2", err
        try:
          pf.append(calibrationfitsfile, np.array([avwav, cal, calerr, res, reserr]), verify=False)
        except: print "whatever3"
    # Do I need residuals? 

    # printmode for minuit =1 in very verbose mode
  
  
  
if __name__ == "__main__":
  sys.exit(main())


# Fancy things to add in when all else is working: 
# # =================
# # = Parse Options =
# # =================
# 
# parser = OptionParser()
# parser.add_option("-e","--exposure", 
#                   action="store", type="string", dest="exponame", 
#                   help="name exposures", metavar="EXPO")
# parser.add_option("-c","--core", 
#                 action="store", type="string", dest="corename", 
#                 help="Dual Core Advantage", metavar="CORE")
# 
# (options, args) = parser.parse_args()
# 
# # Optional use of dual core stuff
# try: core_test += float(options.corename) 
# except: pass
# 
# # Run with flag -c 1    or --core 1
# core1 = i2exposures[:len(i2exposures)/2]
# core2 = i2exposures[len(i2exposures)/2:]
# 
# if core_test == 1:
#     exposures_to_analyze = core1
#     print "Exposures: ", exposures_to_analyze
# if core_test == 2:
#     exposures_to_analyze = core2
#     print "Exposures: ", exposures_to_analyze
# if core_test != 1 and core_test !=2: 
#     exposures_to_analyze = core1 + core2
#     print "Exposures: ", exposures_to_analyze
# 
# # hdulist = pf.open('fxb_ceres_sci_437_01_b_up_001.fits')
# #     for x in hdulist:
#        # ....:     if type(x) == pf.core.ImageHDU:
#        # ....:         wav, flx, err, con = x.data
#        # ....:         print wav[0],wav[-1]
# TODO print wavelength range of order, bin, etc.
# TODO print "Order _ of _" if fits