#!/home/ssi/jwhitmore/progs/epd-7.1-2-rh5-x86_64/bin/python
# 2011-11-16 helpreduce.py

import sys
import shutil
import os
import glob
import csv
import numpy as np
import scipy as sp
import subprocess as subprocess
import time
import datetime
from optparse import OptionParser
import argparse
from ConfigParser import RawConfigParser
import random as ra
import tempfile
import itertools
import fileinput
import shlex
import pylab as pl
import pyfits as pf

import scipy.interpolate as si
import scipy.signal as ss

import minuit as mi 
import scipy.constants as spc

c_light = spc.c

# TODO ceres.cleanup(): invalid value encountered in divide
# TODO overlap bins
# TODO use full order fit results as first-guesses for bins
# TODO figure out how to save results
# TODO check difference between slope/no-slope
# TODO wavelength cut-out regions
# TODO minimum bin-size requirement
# TODO integrate exactly the smooth iodine spectrum

help_message = '''
Various limitations: 
Must have an FTS spectrum w/o gaps
Must have a telescope spectrum w/ increasing wavelength (gaps are OK)
The spacing of the nearest two pixels in the telescope spectrum is used to subsample; so unevenly sub-grid your data at your peril.
'''
# TODO make a flag that says whether each method is actually run.
# TODO implement logging so everything is reproducible. 

class Exposure(object):
  """docstring for Exposure"""
  def __init__(self, arcFile='', reductionProgram='', calibrationType='', calibrationFile='', exposureFile=''):
    """docstring for __init__"""
    super(Exposure, self).__init__()
    self.arcFile = arcFile # a calibration Arc File
    self.exposureFile = exposureFile # a calibration Arc File
    self.reductionProgram = reductionProgram # reduction software used
    self.calibrationType = calibrationType # Calibration type: iodine, asteroid, none
    self.calibrationFile = calibrationFile # Calibration File
    self.fitGuess = { 'fshift':0.02, 'fix_fshift':False, 'limit_fshift':(-1.0,1.0) ,'err_fshift':0.005 }
    self.fitGuess.update({ 'fsigma':10.5, 'fix_fsigma':False, 'limit_fsigma':(2.0,2000) ,'err_fsigma':5 })
    self.fitGuess.update({ 'fmultiple':1.25, 'fix_fmultiple':False, 'limit_fmultiple':(0.1, 100.0) ,'err_fmultiple':0.2 })
    self.fitGuess.update({ 'fslope':0.0005, 'fix_fslope':False, 'limit_fslope':(-1.0,1.0) ,'err_fslope':0.05 })
    self.fitGuess.update({ 'elements':100, 'fix_elements':True })
    self.fitGuess.update({ 'fwidth':200, 'fix_fwidth':True, 'limit_fwidth':(2., self.fitGuess['elements']/2 - 1), 'err_fwidth':5 })
    self.fitGuess.update({ 'strategy':2 })
    self.fitResults = {}
    self.tiltfitResults = {}
    self.BinResults = {}
    self.Bins = {}
    if self.exposureFile.split('.')[-1] == 'fits':
      print "A fits exposure file."
      self.Orders = {}
      hdu = pf.open(self.exposureFile)
      self.header = hdu[0].header
      for i,x in enumerate(hdu):
        self.Orders[i] = {}
        self.Orders[i]['wav'] = x.data[0]
        self.Orders[i]['flx'] = x.data[1]
        self.Orders[i]['err'] = x.data[2]
    else:
      print "Not a fits file.", self.exposureFile
    pass

  def usage(self):
    """docstring for usage"""
    print "The general order goes: "
    print "loadReferenceSpectra, cleanup, continuumFit, chop, overSample, fullOrderShift, binShift"
    pass
  
  def loadReferenceSpectra(self):
    """docstring for loadReferenceSpectra"""
    iow, iof = np.loadtxt(self.calibrationFile, unpack='True')
    print iow[0], iow[-1]
    for x in self.Orders:
      if self.Orders[x]['wav'][0] > iow[0] + 50.0:
        try:
          ok = (self.Orders[x]['wav'][0] - 10 < iow) & (self.Orders[x]['wav'][-1] + 10 > iow)
          if len(iow[ok]) > 200:
            self.Orders[x]['iow'] = iow[ok]
            self.Orders[x]['iof'] = iof[ok]
          "Reference spectra worked"
        except:
          print "Outside overlap."
    pass
  
  def cleanup(self,verbose=False):
    """mask out bad regions of the spectra"""
    # Think about whether to overwrite the input files
    if verbose==True:
      print "Beginning cleanup of data...", datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    errorcutoff = 0.0
    sncutoff = 10.0
    edgebuffer = 50
    for x in self.Orders:
      ok = self.Orders[x]['err'] >= errorcutoff
      ok2 = self.Orders[x]['flx']/self.Orders[x]['err'] >= sncutoff
      for key in ['wav', 'flx', 'err']:
        self.Orders[x][key] = self.Orders[x][key][(ok & ok2)]
    # # TODO Deal with sky lines/known regions to exclude
    # dummywav = cleanwav
    # for x in zip(begin_kill_array, end_kill_array):
    #   dummywav = [y for y in dummywav if not (y > x[0] and y < x[1])]
    # finalindex = [np.argwhere(wav == x)[0][0] for x in dummywav]
    pass
  
  def continuumFit(self, knots=10, plot=False, verbose=False):
    """fits a continuum via a spline through the flux values."""
    knots = 10
    edgeTolerance = 0.1
    for x in self.Orders:
      s = si.LSQUnivariateSpline(self.Orders[x]['wav'],\
                                self.Orders[x]['flx'],\
                                np.linspace(self.Orders[x]['wav'][0]+edgeTolerance, self.Orders[x]['wav'][-1]-edgeTolerance, knots),\
                                w=self.Orders[x]['err'])
      self.Orders[x]['con'] = s(self.Orders[x]['wav']) # new array is made -- continuum
    pass
  
  def chop(self, edgebuffer=50):
    """program chops off the offending beginning and ending few pixels of each order"""
    # edgebuffer = edgebuffer
    for x in self.Orders:
      for key in ['wav', 'flx', 'err', 'con']:
        self.Orders[x][key] = self.Orders[x][key][edgebuffer:-edgebuffer]
    print "Chopped", edgebuffer, "pixels."
    pass
  
  # def overSample(self):
  #   """sets the minimum spacing in telescope spectra:mindel
  #   Also oversamples the FTS spectrum and sets to oiow, oiof."""
  #   for x in self.Orders:
  #     self.Orders[x]['mindel'] = self.Orders[x]['wav'][-1] - self.Orders[x]['wav'][0]
  #     overSampleFactor = 5.0
  #     ftsoverSampleFactor = 10.0
  #     for i in range(len(self.Orders[x]['wav'])-1):
  #       if self.Orders[x]['mindel'] > self.Orders[x]['wav'][i+1] - self.Orders[x]['wav'][i]: 
  #         self.Orders[x]['mindel'] = self.Orders[x]['wav'][i+1] - self.Orders[x]['wav'][i]
  #     # # oiow; oiof
  #     if 'iow' in self.Orders[x]:
  #       iow = self.Orders[x]['iow']
  #       iof = self.Orders[x]['iof']
  #       s = si.UnivariateSpline(iow, iof, s=0)
  #       self.Orders[x]['oiow'] = np.linspace(iow[0], iow[-1], len(iow) * ftsoverSampleFactor)
  #       self.Orders[x]['oiof'] = s(self.Orders[x]['oiow'])
  #     if overSampleFactor % 2 == 1:
  #       subPixelSize = self.Orders[x]['mindel'] / overSampleFactor
  #       self.Orders[x]['overwav'] = []
  #       for y in self.Orders[x]['wav']:
  #         self.Orders[x]['overwav'].append([y + (i - overSampleFactor/2 + 0.5) * self.Orders[x]['mindel'] / overSampleFactor for i in range(int(overSampleFactor))])
  #       self.Orders[x]['overwav'] = np.array(self.Orders[x]['overwav'])
  #     else:
  #       print "overSampleFactor must be an odd number"
  #     pass
  #   pass
  
  def newOverSample(self):
    """sets the minimum spacing in the telescope spectra (mindel) for each order over the whole exposure.
    Rename. """
    for x in self.Orders:
      self.Orders[x]['mindel'] = self.Orders[x]['wav'][-1] - self.Orders[x]['wav'][0]
      for i in range(len(self.Orders[x]['wav']) - 1):
        if self.Orders[x]['mindel'] > self.Orders[x]['wav'][i+1] - self.Orders[x]['wav'][i]: 
          self.Orders[x]['mindel'] = self.Orders[x]['wav'][i+1] - self.Orders[x]['wav'][i]
    pass
  
  def fullExposureShift(self, verbose=False, veryVerbose=False, robustSearch=False, elements=1000, sigma=50):
    """docstring for fullExposureShift"""
    for x in self.Orders:
      if 'oiow' in self.Orders[x]:
        print "Working on order: ", x
        self.dictionaryOrderShift(order=x, verbose=verbose, veryVerbose=veryVerbose, robustSearch=robustSearch)
    pass

  def dictionaryShift(self, order, fmultiple, fshift, fsigma, elements, **kwargs):
    """docstring for dictionaryShift"""
    kernel = self.gaussKernel(self.fitGuess['elements'], fsigma)
    tck = si.splrep(self.Orders[order]['oiow'], np.convolve(kernel, self.Orders[order]['oiof'] * fmultiple, mode='same'))
    overflx = np.average(si.splev(np.hstack(self.Orders[order]['overwav']) + fshift, tck).reshape(np.shape(self.Orders[order]['overwav'])), axis=1)
    return np.sum( ((overflx - self.Orders[order]['flx'] / self.Orders[order]['con']) / \
                    (self.Orders[order]['err'] / self.Orders[order]['con']) ) ** 2 )

  
  def dictionaryOrderShift(self, order=7, verbose=False, veryVerbose=False, robustSearch=False):
    """docstring for dictionaryShift"""
    m = mi.Minuit(self.dictionaryShift, order=order, fix_order=True, **self.fitGuess)
    if veryVerbose==True:
      m.printMode=1
    if robustSearch==True:
      print "Robust search. Beginning initial scan..."
      m.scan(("fshift",20,-0.5,0.5))
      print "done."
    try: 
      print "Finding initial shift/fit", '\n', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") 
      m.migrad()
      self.fitResults[str(order)] = {}
      self.fitResults[str(order)]['values'] = m.values
      self.fitResults[str(order)]['errors'] = m.errors
    except:
      print "Serious problem with order:", order
    pass

  # def OrderShiftandTilt(self, order=7, verbose=False, veryVerbose=False, robustSearch=False):
  #   """docstring for dictionaryShift"""
  #   m = mi.Minuit(self.shiftandtilt, order=order, fix_order=True, **self.fitGuess)
  #   if veryVerbose==True:
  #     m.printMode=1
  #   if robustSearch==True:
  #     print "Robust search. Beginning initial scan..."
  #     m.scan(("fshift",20,-0.5,0.5))
  #     print "done."
  #   try: 
  #     print "Finding initial shift/fit", '\n', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") 
  #     m.migrad()
  #     self.tiltfitResults[str(order)] = {}
  #     self.tiltfitResults[str(order)]['values'] = m.values
  #     self.tiltfitResults[str(order)]['errors'] = m.errors
  #   except:
  #     print "Serious problem with order:", order
  #   pass
  
  def newOrderShiftandTilt(self, order=7, verbose=False, veryVerbose=False, robustSearch=False):
    """docstring for dictionaryShift"""
    m = mi.Minuit(self.newshiftandtilt, order=order, fix_order=True, **self.fitGuess)
    if veryVerbose==True:
      m.printMode=1
    if robustSearch==True:
      print "Robust search. Beginning initial scan..."
      m.scan(("fshift",20,-0.5,0.5))
      print "done."
    # try: 
    print "Finding initial shift/fit", '\n', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") 
    m.migrad()
    self.tiltfitResults[str(order)] = {}
    self.tiltfitResults[str(order)]['values'] = m.values
    self.tiltfitResults[str(order)]['errors'] = m.errors
    # except:
      # print "Serious problem with order:", order
    pass


  # 
  # def createMinuit(minimizationFunction, **kwargs):
  #   """Creates a Minuit object with the **kwargs"""
  #   print str(minimizationFunction) + ',',
  #   for name, value in kwargs.items():
  #   #   
  #   return name = value[0] 
  #   # + ', fix_' + name + '=' + str(value[1]) + ', limit_' + name + '=' + str(value[2]) + ', err_' + name + '=' + str(value[3]) + ', strategy=2)'#m.Minuit(minimizationFunction, )
  # 
  # pass dictionary to function like: createMinuit(**theDictionary)
  # linear dispersion coefficient
  # spectral line depth
  # normalization: 
  # 2nd-order dispersion coefficient
  # width of main Gaussian IP
  # width of box IP
  # residuals
  # plot kernel
  # plot best fit between the two
  
  def KernelFunction(self, kernelType='gauss', elements=1000, sigma=50, width=200):
    """returns whatever kernel is being called"""
    if kernelType=='gauss':
      return self.gaussKernel(elements, sigma)
    elif kernelType=='box':
      return self.boxKernel(elements, width)
    elif kernelType=='gaussbox':
      return self.gaussboxKernel(elements, sigma, width)
    else:
      return "Not known."

  def gaussKernel(self, elements, sigma):
    """returns a normalized gaussian using scipy.signal"""
    return ss.gaussian(elements, sigma) / np.sum(ss.gaussian(elements, sigma))
  
  def boxKernel(self, elements, width):
    """docstring for boxKernel"""
    gbox = np.zeros(elements)
    if width > elements/2:
      print "Width has to be less than half the window"
      return 1
    gbox[elements/2 - width:elements/2 + width] = 1.0
    return gbox / np.sum(gbox)

  def gaussboxKernel(self, elements, sigma, boxwidth):
    """returns normalized box + gaussian kernel."""
    kernel = self.boxKernel(elements, boxwidth) + self.gaussKernel(elements, sigma)
    return kernel / np.sum(kernel)

  def shift(self, order, fmultiple, fshift, fsigma, elements):
    """docstring for shift"""
    kernel = self.gaussKernel(elements, fsigma)
    tck = si.splrep(self.Orders[order]['oiow'], np.convolve(kernel, self.Orders[order]['oiof'] * fmultiple, mode='same'))
    overflx = np.average(si.splev(np.hstack(self.Orders[order]['overwav']) + fshift, tck).reshape(np.shape(self.Orders[order]['overwav'])), axis=1)
    return np.sum( (overflx - self.Orders[order]['flx'] / self.Orders[order]['con']) ** 2 / \
                    (self.Orders[order]['err'] / self.Orders[order]['con']) ** 2 )

  def shiftandtilt(self, order, fmultiple, fshift, fsigma, elements, fslope, **kwargs):
    """Fit like shift with the addition of a slope across the order."""
    # continuum slope: wav, flx + slope * (wav - np.average(wav))
    kernel = self.gaussKernel(elements, fsigma)
    tck = si.splrep(self.Orders[order]['oiow'], np.convolve(kernel, (self.Orders[order]['oiof'] * fmultiple) + fslope * (self.Orders[order]['oiow'] - np.average(self.Orders[order]['oiow'])), mode='same'))
    overflx = np.average(si.splev(np.hstack(self.Orders[order]['overwav']) + fshift, tck).reshape(np.shape(self.Orders[order]['overwav'])), axis=1)
    return np.sum( ((overflx - self.Orders[order]['flx'] / self.Orders[order]['con']) / \
                    (self.Orders[order]['err'] / self.Orders[order]['con'])) ** 2 )
  
  def newshiftandtilt(self, order, fmultiple, fshift, fsigma, elements, fslope, **kwargs):
    """trying to smooth, interpolate, and integrate the fit."""
    kernel = self.gaussKernel(elements, fsigma)
    s = si.UnivariateSpline(self.Orders[order]['iow'], np.convolve(kernel, (self.Orders[order]['iof'] * fmultiple) + fslope * (self.Orders[order]['iow'] - np.average(self.Orders[order]['iow'])), mode='same'), s=0)
    overflx = np.array([s.integral(x - self.Orders[order]['mindel']/2.0 + fshift, x + self.Orders[order]['mindel']/2.0 + fshift) for x in self.Orders[order]['wav']])
    # print overflx - self.Orders[order]['flx'] / self.Orders[order]['con']
    return np.sum( ((overflx - self.Orders[order]['flx'] / self.Orders[order]['con']) / \
                    (self.Orders[order]['err'] / self.Orders[order]['con'])) ** 2 )
  
  def gaussboxShift(self, order, fmultiple, fshift, fsigma, elements, fwidth):
    """docstring for gaussboxShift"""
    kernel = self.gaussboxKernel(elements, fsigma, fwidth)
    tck = si.splrep(self.Orders[order]['oiow'], np.convolve(kernel, self.Orders[order]['oiof'] * fmultiple, mode='same'))
    overflx = np.average(si.splev(np.hstack(self.Orders[order]['overwav']) + fshift, tck).reshape(np.shape(self.Orders[order]['overwav'])), axis=1)
    return np.sum( (overflx - self.Orders[order]['flx'] / self.Orders[order]['con']) ** 2 / \
                    (self.Orders[order]['err'] / self.Orders[order]['con']) ** 2 )
  
  def createBinArrays(self, order=7, binSize=350):
    """docstring for createBinArrays"""
    # create possibility of overlap
    lamb = np.average(self.Orders[order]['wav'])
    # TODO check if already exists
    try:
      type(self.fitResults[binSize])
    except:
      self.fitResults[binSize] = {}
    try:
      type(self.Orders[order][binSize])
      return
    except:
      self.Orders[order][binSize] = {}
    # create bin wavelength bounds
    self.Orders[order][binSize]['binEdges'] = np.arange(self.Orders[order]['wav'][0], self.Orders[order]['wav'][-1], lamb * binSize * 1000 / c_light)
    self.Orders[order][binSize]['binEdges'] = np.append(self.Orders[order][binSize]['binEdges'], self.Orders[order]['wav'][-1]) # add last wavelength point to edges)
    oiowTolerance = 2.0
    self.Orders[order][binSize]['bins'] = {}
    for i in range(len(self.Orders[order][binSize]['binEdges']) - 1):
      self.Orders[order][binSize]['bins'][i] = {}
      self.Orders[order][binSize]['bins'][i]['ok'] = (self.Orders[order]['wav'] > self.Orders[order][binSize]['binEdges'][i]) & (self.Orders[order]['wav'] <= self.Orders[order][binSize]['binEdges'][i + 1])
      self.Orders[order][binSize]['bins'][i]['iok'] = (self.Orders[order]['oiow'] > self.Orders[order][binSize]['binEdges'][i] - oiowTolerance) & (self.Orders[order]['oiow'] <= self.Orders[order][binSize]['binEdges'][i + 1] + oiowTolerance)
    pass

  def newCreateBinArrays(self, order=7, binSize=350, overlap=0.5):
    """overlap is the fractional overlap or how much the bin is shifted relative to the binSize. so overlapping by .5 shifts by half binSize; .33 by .33 binSize. """
    lamb = np.average(self.Orders[order]['wav'])
    try:
      type(self.fitResults[binSize])
    except:
      self.fitResults[binSize] = {}
    try:
      type(self.Orders[order][binSize])
      return
    except:
      self.Orders[order][binSize] = {}
    binAngstroms = lamb * binSize * 1000 / c_light
    temp = []
    for x in range(int(1.0/overlap)):
      temp.append(np.arange(self.Orders[order]['wav'][0] + overlap * x * binAngstroms, self.Orders[order]['wav'][-1] + overlap * x * binAngstroms, binAngstroms))

    np.append(temp[0], self.Orders[order]['wav'][-1]) # add last wavelength point to first bin edges array
    oiowTolerance = 2.0
    self.Orders[order][binSize]['bins'] = {}
    COUNTER = 0
    for edgearray in temp:
      for i in range(len(edgearray) - 1):
        self.Orders[order][binSize]['bins'][COUNTER] = {}
        self.Orders[order][binSize]['bins'][COUNTER]['ok'] = (self.Orders[order]['wav'] > edgearray[i]) & (self.Orders[order]['wav'] <= edgearray[i + 1])
        self.Orders[order][binSize]['bins'][COUNTER]['iok'] = (self.Orders[order]['iow'] > edgearray[i] - oiowTolerance) & (self.Orders[order]['iow'] <= edgearray[i + 1] + oiowTolerance)
        COUNTER += 1
    pass
  
  def fullOrderBinShift(self, order=7, binSize=350):
    """docstring for fullOrderBinShift"""
    # TODO check if createBinArrays has been run; if not; run first...
    try:
      type(self.fitResults[binSize][order])
    except:
      self.fitResults[binSize][order] = {}
    try:
      type(self.fitResults[binSize][order]['bins'])
    except:
      self.fitResults[binSize][order]['bins'] = {}
    for bin in self.Orders[order][binSize]['bins']:
      self.fitResults[binSize][order]['bins'][bin] = {}
      self.smallBinShift(order, binSize, bin)
    pass

  def binshiftandtilt(self, order, bin, binSize, fmultiple, fshift, fsigma, elements, fslope, **kwargs):
    """Fit like shift with the addition of a slope across the order."""
    kernel = self.gaussKernel(elements, fsigma)
    ok = self.Orders[order][binSize]['bins'][bin]['ok']
    iok = self.Orders[order][binSize]['bins'][bin]['iok']
    tck = si.splrep(self.Orders[order]['oiow'][iok], np.convolve(kernel, (self.Orders[order]['oiof'][iok] * fmultiple) + fslope * (self.Orders[order]['oiow'][iok] - np.average(self.Orders[order]['oiow'][iok])), mode='same'))
    overflx = np.average(si.splev(np.hstack(self.Orders[order]['overwav'][ok]) + fshift, tck).reshape(np.shape(self.Orders[order]['overwav'][ok])), axis=1)
    return np.sum( ((overflx - self.Orders[order]['flx'][ok] / self.Orders[order]['con'][ok]) / \
                    (self.Orders[order]['err'][ok] / self.Orders[order]['con'][ok])) ** 2 )
  
  def smallBinShift(self, order=7, binSize=350, bin=2, veryVerbose=False, robustSearch=False):
    """docstring for smallBinShift"""
    m = mi.Minuit(self.binshiftandtilt, order=order, binSize=binSize, bin=bin, fix_order=True, fix_binSize=True, fix_bin=True, **self.fitGuess)
    if veryVerbose==True:
      m.printMode=1
    if robustSearch==True:
      print "Robust search. Beginning initial scan..."
      m.scan(("fshift",20,-0.5,0.5))
      print "done."
    try: 
      # print "Finding initial shift/fit", '\n', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") 
      print datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "Finding initial shift/fit for order:", order, "and bin:", bin
      m.migrad()
      self.fitResults[binSize][order]['bins'][bin]['values'] = m.values
      self.fitResults[binSize][order]['bins'][bin]['errors'] = m.errors
      ok = self.Orders[order][binSize]['bins'][bin]['ok']
      elements = self.fitGuess['elements']
      iok = self.Orders[order][binSize]['bins'][bin]['iok']
      lamb = np.average(self.Orders[order]['wav'][ok])
      cal = m.values['fshift'] * c_light / lamb
      calerr = m.errors['fshift'] * c_light / lamb
      midpointFTS = np.argmin( np.abs(self.Orders[order]['oiow'][iok] - lamb) )
      FTSchunk = self.Orders[order]['oiow'][iok][midpointFTS + elements/2] - self.Orders[order]['oiow'][iok][midpointFTS - elements/2]
      FTSsigma = FTSchunk * m.values['fsigma'] / elements # size of sigma in wavelength
      FWHM = 2.0 * np.sqrt( 2.0 * np.log(2.0) ) * FTSsigma # size of FWHM in wavelength
      R = lamb / FWHM
      self.fitResults[binSize][order]['bins'][bin]['avwav'] = lamb
      self.fitResults[binSize][order]['bins'][bin]['cal'] = cal
      self.fitResults[binSize][order]['bins'][bin]['calerr'] = calerr
      self.fitResults[binSize][order]['bins'][bin]['R'] = R
      print datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "finished."
    except:
      print "Serious problem with bin:", bin
    pass
  
  def plotInitialGuess(self, order, fmultiple, fshift, fsigma, elements=1000, sigma=50):
    """docstring for plotInitialGuess"""
    kernel = self.gaussKernel(elements, fsigma)
    tck = si.splrep(self.Orders[order]['oiow'], np.convolve(kernel, self.Orders[order]['oiof'] * fmultiple, mode='same'))
    overflx = np.average(si.splev(np.hstack(self.Orders[order]['overwav']) + fshift, tck).reshape(np.shape(self.Orders[order]['overwav'])), axis=1)
    pl.plot(self.Orders[order]['wav'], self.Orders[order]['flx'] / self.Orders[order]['con'], color="black", linewidth=2.0)
    pl.plot(np.average(self.Orders[order]['overwav'],axis=1), overflx)
    pass
  
  def plotFitResults(self, order, fmultiple, fshift, fsigma, elements=1000, **kwargs):
    """docstring for plotFitResults"""
    kernel = self.gaussKernel(elements, fsigma)
    tck = si.splrep(self.Orders[order]['oiow'], np.convolve(kernel, self.Orders[order]['oiof'] * fmultiple, mode='same'))
    overflx = np.average(si.splev(np.hstack(self.Orders[order]['overwav']) + fshift, tck).reshape(np.shape(self.Orders[order]['overwav'])), axis=1)
    pl.plot(self.Orders[order]['wav'], self.Orders[order]['flx'] / self.Orders[order]['con'], color="black", linewidth=2.0)
    pl.plot(np.average(self.Orders[order]['overwav'],axis=1), overflx)    
    pass
  
  def plotTiltFitResults(self, order, fmultiple, fshift, fsigma, fslope, elements=1000, plotResiduals=False, **kwargs):
    """docstring for plotTiltFitResults"""
    kernel = self.gaussKernel(elements, fsigma)
    tck = si.splrep(self.Orders[order]['oiow'], np.convolve(kernel, (self.Orders[order]['oiof'] * fmultiple) + fslope * (self.Orders[order]['oiow'] - np.average(self.Orders[order]['oiow'])), mode='same'))
    overflx = np.average(si.splev(np.hstack(self.Orders[order]['overwav']) + fshift, tck).reshape(np.shape(self.Orders[order]['overwav'])), axis=1)
    pl.plot(self.Orders[order]['wav'], self.Orders[order]['flx'] / self.Orders[order]['con'], color="black", linewidth=2.0)
    pl.plot(np.average(self.Orders[order]['overwav'],axis=1), overflx)
    if plotResiduals == True:
      pl.plot(self.Orders[order]['wav'], self.Orders[order]['flx'] / self.Orders[order]['con'] - overflx, color="red") # data - model
    pass
  
  def plotBinTiltFitResults(self, order, fmultiple, fshift, fsigma, fslope, binSize, binNumber, elements=1000, plotResiduals=False, **kwargs):
    """docstring for plotBinTiltFitResults"""
    kernel = self.gaussKernel(elements, fsigma)
    ok = self.Orders[order][binSize]['bins'][binNumber]['ok']
    iok = self.Orders[order][binSize]['bins'][binNumber]['iok']
    tck = si.splrep(self.Orders[order]['oiow'][iok], np.convolve(kernel, (self.Orders[order]['oiof'][iok] * fmultiple) + fslope * (self.Orders[order]['oiow'][iok] - np.average(self.Orders[order]['oiow'][iok])), mode='same'))
    overflx = np.average(si.splev(np.hstack(self.Orders[order]['overwav'][ok]) + fshift, tck).reshape(np.shape(self.Orders[order]['overwav'][ok])), axis=1)
    pl.plot(self.Orders[order]['wav'][ok], self.Orders[order]['flx'][ok] / self.Orders[order]['con'][ok], color="black", linewidth=2.0)
    pl.plot(np.average(self.Orders[order]['overwav'][ok], axis=1), overflx)
    if plotResiduals == True: 
      pl.plot(self.Orders[order]['wav'][ok], self.Orders[order]['flx'][ok] / self.Orders[order]['con'][ok] - overflx, color="red")
    pass

  def plotOrderBinTiltFitResults(self):
    """docstring for plotOrderBinTiltFitResults"""
    # TODO test for whether order fit; then plot that order
    for x in range(7):
      ceres.plotBinTiltFitResults(order=7, binNumber=x, **ceres.fitResults[350][7]['bins'][x])
    pass
    
  def expplot(self):
    """docstring for plot"""
    print "working..."
    for x in self.Orders:
      pl.plot(self.Orders[x]['wav'], self.Orders[x]['flx'])
    pl.savefig('ordersinexposure.pdf')
    pl.close()
    pass

  def expftsplot(self):
    """docstring for expftsplot"""
    for x in self.Orders:
      if 'iow' in self.Orders[x]:
        pl.plot(self.Orders[x]['wav'], self.Orders[x]['flx'])
        pl.plot(self.Orders[x]['iow'], self.Orders[x]['iof'])
    pl.savefig('ftsandexposure.pdf')
    pl.close()
    pass
  
  pass

    

def main(argv=None):
  pass

if __name__ == "__main__":
  # sys.exit(main())
  main()
