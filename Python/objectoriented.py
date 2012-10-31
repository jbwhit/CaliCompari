#!/usr/bin/env python

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

import sys
import shutil
import os
import glob
import csv
import pprint 
import numpy as np
import scipy as sp
import subprocess as subprocess
import time
import datetime
from optparse import OptionParser
import argparse
from ConfigParser import RawConfigParser, SafeConfigParser
import random as ra
import tempfile
import itertools
import fileinput
import shlex
import pylab as pl
import pyfits as pf
import cPickle as pickle
import scipy.interpolate as si
import scipy.signal as ss

import minuit as mi 
import scipy.constants as spc

c_light = spc.c

# TODO wavelength cut-out regions
# TODO minimum bin-size requirement
# TODO think about integration interval (+/- 1/2 mindel)
# TODO make a flag that says whether each method is actually run.
# TODO implement logging so everything is reproducible.

# TODO create an 
# exposure dictionary[order][wav/flx/err/pixel/masks/whateverelse]
# fit dictionary
# binFit dictionary
# Find a way to save original exposure header and original ThAr calibration header

help_message = '''
Various limitations: 
Must have an FTS spectrum w/o gaps
Must have a telescope spectrum w/ monotonically increasing wavelength (gaps are OK)
The spacing of the nearest two pixels in the telescope spectrum is used as the pixel size for each order.
'''

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
    self.fitGuess = {}
    self.fitGuess['initial'] = { 'fshift':0.002, 'fix_fshift':False, 'limit_fshift':(-1.0,1.0) ,'err_fshift':0.005 }
    self.fitGuess['initial'].update({ 'fsigma':10.5, 'fix_fsigma':False, 'limit_fsigma':(2.0,30.0) ,'err_fsigma':0.5 })
    self.fitGuess['initial'].update({ 'fmultiple':50.25, 'fix_fmultiple':False, 'limit_fmultiple':(0.1, 100.0) ,'err_fmultiple':0.2 })
    # self.fitGuess['initial'].update({ 'fslope':0.0005, 'fix_fslope':True, 'limit_fslope':(-1.0,1.0) ,'err_fslope':0.05 })
    self.fitGuess['initial'].update({ 'fslope':0.0005, 'fix_fslope':False, 'limit_fslope':(-1.0,1.0) ,'err_fslope':0.05 })
    self.fitGuess['initial'].update({ 'elements':100, 'fix_elements':True })
    self.fitGuess['initial'].update({ 'fwidth':200, 'fix_fwidth':True })
    self.fitGuess['initial'].update({ 'strategy':2 })
    self.fitResults = {}
    if self.exposureFile.split('.')[-1] == 'fits':
      print "A fits exposure file."
      self.Orders = {}
      hdu = pf.open(self.exposureFile)
      self.header = hdu[0].header
      for i,x in enumerate(hdu):
        try:
          type(hdu[i].data)
          self.Orders[i] = {}
          self.Orders[i]['wav'] = x.data[0]
          self.Orders[i]['flx'] = x.data[1]
          self.Orders[i]['err'] = x.data[2]
          self.Orders[i]['pix'] = np.arange(len(self.Orders[i]['wav']))
        except:
          self.exposureHeader = hdu[-1].header
      for field in self.Orders.keys():
        if len(self.Orders[field]) < 1:
          del(self.Orders[field])
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
    try: 
      iow, iof = np.loadtxt(self.calibrationFile)
    except:
      print "Consider saving a faster-loading calibration file."
      iow, iof = np.loadtxt(self.calibrationFile, unpack='True')
    print "Reference FTS wavelength range:", iow[0], iow[-1]
    for x in self.Orders:
      if (self.Orders[x]['wav'][0] > iow[0] + 40.0) & (self.Orders[x]['wav'][-1] < iow[-1] - 150.0):
        try:
          ok = (self.Orders[x]['wav'][0] - 10 < iow) & (self.Orders[x]['wav'][-1] + 10 > iow)
          if len(iow[ok]) > 200:
            self.Orders[x]['iow'] = iow[ok]
            self.Orders[x]['iof'] = iof[ok]
          "Reference spectra worked"
        except:
          print "Order", x, "is outside overlap with reference FTS."
    pass
  
  def cleanup(self,verbose=False):
    """mask out bad regions of the spectra
    Example config file setup. 
    [skylines]
    remove: 
      5589.128    5589.132
      5865.454    5865.459
    """
    parser = SafeConfigParser()
    candidates = glob.glob('config*')
    found = parser.read(candidates)
    wavekill = parser.get('skylines','remove')
    if verbose==True:
      print "Beginning cleanup of data...", datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    errorcutoff = 0.0
    flxcutoff = 0.0
    sncutoff = 10.0
    for x in self.Orders:
      masks = []
      masks.append(self.Orders[x]['err'] > errorcutoff)
      masks.append(self.Orders[x]['flx'] > flxcutoff)
      masks.append(self.Orders[x]['flx']/self.Orders[x]['err'] >= sncutoff)
      for killLine in wavekill.splitlines():
        if len(killLine) > 1:
          masks.append(reduce(np.logical_or, [self.Orders[x]['wav'] < float(killLine.split()[0]), self.Orders[x]['wav'] > float(killLine.split()[1])]))
      self.Orders[x]['mask'] = reduce(np.logical_and, masks)
    pass
  
  def continuumFit(self, knots=10, plot=False, verbose=False):
    """fits a continuum via a spline through the flux values."""
    knots = 10
    edgeTolerance = 0.1
    for x in self.Orders:
      mask = self.Orders[x]['mask']
      self.Orders[x]['con'] = np.zeros(len(self.Orders[x]['wav']))
      s = si.LSQUnivariateSpline(self.Orders[x]['wav'][mask],\
                                self.Orders[x]['flx'][mask],\
                                np.linspace(self.Orders[x]['wav'][mask][0]+edgeTolerance, self.Orders[x]['wav'][mask][-1]-edgeTolerance, knots),\
                                w=self.Orders[x]['err'][mask])
      self.Orders[x]['con'][mask] = s(self.Orders[x]['wav'][mask]) # new array is made -- continuum
    pass
  
  def OverSample(self):
    """sets the minimum spacing in the telescope spectra (mindel) for each order over the whole exposure.
    Rename. """
    for x in self.Orders:
      mask = self.Orders[x]['mask']
      self.Orders[x]['mindel'] = self.Orders[x]['wav'][mask][-1] - self.Orders[x]['wav'][mask][0]
      for i in range(len(self.Orders[x]['wav'][mask]) - 1):
        if self.Orders[x]['mindel'] > self.Orders[x]['wav'][mask][i+1] - self.Orders[x]['wav'][mask][i]: 
          self.Orders[x]['mindel'] = self.Orders[x]['wav'][mask][i+1] - self.Orders[x]['wav'][mask][i]
    pass
  
  def fullExposureShift(self, verbose=False, veryVerbose=False, robustSearch=False, binSize=350):
    """docstring for fullExposureShift"""
    starttime=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    for x in self.Orders:
      if 'iow' in self.Orders[x]:
        print "Working on order: ", x
        self.CreateBinArrays(order=x, binSize=binSize) # new!
        try:
          self.OrderShiftandTilt(order=x, veryVerbose=veryVerbose) # new!
          self.fullOrderBinShift(order=x, binSize=binSize)
        except:
          print "Order or bin failed."
    print "Finished working on exposure."
    print "Started: ", starttime, "Ended: ", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    pass

  def OrderShiftandTilt(self, order=7, verbose=False, veryVerbose=False, robustSearch=False):
    """docstring for dictionaryShift"""
    try:
      type(self.fitResults['order'])
    except:
      self.fitResults['order'] = {}
    try:
      type(self.fitResults['order'][order])
    except:
      self.fitResults['order'][order] = {}
    try:
      m = mi.Minuit(self.shiftandtilt, order=order, fix_order=True, **self.fitGuess['initial'])
      if veryVerbose==True:
        m.printMode=1
      if robustSearch==True:
        print "Robust search. Beginning initial scan..."
        m.scan(("fshift",20,-0.5,0.5))
        print "done."
      # try: 
      print "Finding initial full order shift/fit", '\n', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") 
      m.migrad()
      self.fitResults['order'][order]['values'] = m.values
      try: 
        del self.fitResults['order'][order]['values']['order']
      except:
        pass
      self.fitResults['order'][order]['errors'] = m.errors
    except:
      print "Serious problem with order:", order
    pass

  def gaussKernel(self, elements, sigma):
    """returns a normalized gaussian using scipy.signal"""
    return ss.gaussian(elements, sigma) / np.sum(ss.gaussian(elements, sigma))
  
  def shiftandtilt(self, order, fmultiple, fshift, fsigma, elements, fslope, **kwargs):
    """trying to smooth, interpolate, and integrate the fit."""
    mask = self.Orders[order]['mask']
    kernel = self.gaussKernel(elements, fsigma)
    s = si.UnivariateSpline(self.Orders[order]['iow'], np.convolve(kernel, (self.Orders[order]['iof'] * fmultiple) + fslope * (self.Orders[order]['iow'] - np.average(self.Orders[order]['iow'])), mode='same'), s=0)
    overflx = np.array([s.integral(x - self.Orders[order]['mindel']/2.0 + fshift, x + self.Orders[order]['mindel']/2.0 + fshift) for x in self.Orders[order]['wav'][mask]])
    return np.sum( ((overflx - self.Orders[order]['flx'][mask] / self.Orders[order]['con'][mask])  / \
                    (self.Orders[order]['err'][mask] / self.Orders[order]['con'][mask])) ** 2)


  # TODO add masks to everything!
  def CreateBinArrays(self, order=7, binSize=350, overlap=0.5, iowTolerance=2.0, minPixelsPerBin=100):
    """overlap is the fractional overlap or how much the bin is shifted relative to the binSize. so overlapping by .5 shifts by half binSize; .33 by .33 binSize. """
    mask = self.Orders[order]['mask']
    lamb = np.average(self.Orders[order]['wav'][mask])
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
    mask = self.Orders[order]['mask']
    # TODO make sure bin has 'enough' pixels to be useful.
    for x in range(int(1.0/overlap)):
      temp.append(np.arange(self.Orders[order]['wav'][mask][0] + overlap * x * binAngstroms, self.Orders[order]['wav'][mask][-1] + overlap * x * binAngstroms, binAngstroms))
    np.append(temp[0], self.Orders[order]['wav'][mask][-1]) # add last wavelength point to first bin edges array
    iowTolerance = iowTolerance
    minPixelsPerBin = minPixelsPerBin
    COUNTER = 0
    for edgearray in temp:
      for i in range(len(edgearray) - 1):
        if len(self.Orders[order]['wav'][(self.Orders[order]['wav'] > edgearray[i]) & (self.Orders[order]['wav'] <= edgearray[i + 1])]) > minPixelsPerBin:
          self.Orders[order][binSize][COUNTER] = {}
          self.Orders[order][binSize][COUNTER]['ok'] = (self.Orders[order]['wav'] > edgearray[i]) & (self.Orders[order]['wav'] <= edgearray[i + 1])
          self.Orders[order][binSize][COUNTER]['iok'] = (self.Orders[order]['iow'] > edgearray[i] - iowTolerance) & (self.Orders[order]['iow'] <= edgearray[i + 1] + iowTolerance)
          COUNTER += 1
        else:
          print "Bin ", i, " would have had less than ", minPixelsPerBin, " -- not creating a bin for it."
    pass
  
  def fullOrderBinShift(self, order=7, binSize=350):
    """docstring for fullOrderBinShift"""
    # TODO check if createBinArrays has been run; if not; run first...
    try:
      type(self.fitResults[binSize])
    except:
      self.fitResults[binSize] = {}
    try:
      type(self.fitResults[binSize][order])
    except:
      self.fitResults[binSize][order] = {}
    try:
      type(self.fitGuess['order'])
    except:
      self.fitGuess['order'] = {}
    try:
      type(self.fitGuess['order'][order])
    except:
      self.fitGuess['order'][order] = {}
    try:
      del self.fitGuess['initial']['order']
    except: 
      pass
    self.fitGuess['order'][order] = self.fitGuess['initial']
    self.fitGuess['order'][order].update(self.fitResults['order'][order]['values'])
    self.fitGuess['order'][order].update({ 'elements':int(10.0 * self.fitResults['order'][order]['values']['fsigma']) })
    for singlebin in self.Orders[order][binSize]:
      self.fitResults[binSize][order][singlebin] = {}
      self.smallBinShift(order, binSize, singlebin)
    pass


  def smallBinShift(self, order=7, binSize=350, bin=2, veryVerbose=False, robustSearch=False):
    """docstring for smallBinShift"""
    # TODO check that the full order solution has run.
    try:
      type(self.fitResults['order'][order]['values'])
    except:
      print "It doesn't look like the full order was run... "
    m = mi.Minuit(self.binshiftandtilt, order=order, binSize=binSize, bin=bin, fix_order=True, fix_binSize=True, fix_bin=True, **self.fitGuess['order'][order])
    if veryVerbose==True:
      m.printMode=1
    if robustSearch==True:
      print "Robust search. Beginning initial scan..."
      m.scan(("fshift",20,-0.5,0.5))
      print "done."
    try: 
      print datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "Finding initial shift/fit for order:", order, "and bin:", bin
      m.migrad()
      self.fitResults[binSize][order][bin]['values'] = m.values
      self.fitResults[binSize][order][bin]['errors'] = m.errors
      mask = self.Orders[order]['mask']
      ok = reduce(np.logical_and, [self.Orders[order][binSize][bin]['ok'], mask])
      iok = self.Orders[order][binSize][bin]['iok']
      elements = self.fitResults[binSize][order][bin]['values']['elements']
      lamb = np.average(self.Orders[order]['wav'][ok])
      avpix = np.average(self.Orders[order]['pix'][ok])
      cal = m.values['fshift'] * c_light / lamb
      calerr = m.errors['fshift'] * c_light / lamb
      midpointFTS = np.argmin(np.abs(self.Orders[order]['iow'][iok] - lamb))
      FTSchunk = self.Orders[order]['iow'][iok][midpointFTS + elements/2] - self.Orders[order]['iow'][iok][midpointFTS - elements/2]
      FTSsigma = FTSchunk * m.values['fsigma'] / elements # size of sigma in wavelength
      FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0)) * FTSsigma # size of FWHM in wavelength
      R = lamb / FWHM
      posFTSsigma = FTSsigma + m.errors['fsigma'] / elements # positive error
      negFTSsigma = FTSsigma - m.errors['fsigma'] / elements # negative error
      Rsmall = lamb / (2.0 * np.sqrt(2.0 * np.log(2.0)) * posFTSsigma)
      Rbig = lamb / (2.0 * np.sqrt(2.0 * np.log(2.0)) * negFTSsigma)
      self.fitResults[binSize][order][bin]['avwav'] = lamb
      self.fitResults[binSize][order][bin]['cal'] = cal
      self.fitResults[binSize][order][bin]['calerr'] = calerr
      self.fitResults[binSize][order][bin]['R'] = R
      self.fitResults[binSize][order][bin]['Rsmall'] = R - Rsmall
      self.fitResults[binSize][order][bin]['Rbig'] = Rbig - R
      self.fitResults[binSize][order][bin]['avpix'] = avpix
      print datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "finished."
    except:
      # TODO flag bin as bad.
      print "Serious problem with bin:", bin
    pass

  def binshiftandtilt(self, order, bin, binSize, fmultiple, fshift, fsigma, elements, fslope, **kwargs):
    """Fit like shift with the addition of a slope across the order."""
    mask = self.Orders[order]['mask']
    kernel = self.gaussKernel(elements, fsigma)
    ok = reduce(np.logical_and, [self.Orders[order][binSize][bin]['ok'], mask])
    iok = self.Orders[order][binSize][bin]['iok']
    s = si.UnivariateSpline(self.Orders[order]['iow'][iok], np.convolve(kernel, (self.Orders[order]['iof'][iok] * fmultiple) + fslope * (self.Orders[order]['iow'][iok] - np.average(self.Orders[order]['iow'][iok])), mode='same'), s=0)
    overflx = np.array([s.integral(x - self.Orders[order]['mindel']/2.0 + fshift, x + self.Orders[order]['mindel']/2.0 + fshift) for x in self.Orders[order]['wav'][ok]])
    return np.sum( ((overflx - self.Orders[order]['flx'][ok] / self.Orders[order]['con'][ok]) / \
                    (self.Orders[order]['err'][ok] / self.Orders[order]['con'][ok])) ** 2 )

  def rewritebinshiftandtilt(self, order, bin, binSize, fmultiple, fshift, fsigma, elements, fslope, **kwargs):
    """Same as before but returns a tuple of results"""
    mask = self.Orders[order]['mask']
    kernel = self.gaussKernel(elements, fsigma)
    ok = reduce(np.logical_and, [self.Orders[order][binSize][bin]['ok'], mask])
    iok = self.Orders[order][binSize][bin]['iok']
    s = si.UnivariateSpline(self.Orders[order]['iow'][iok], np.convolve(kernel, (self.Orders[order]['iof'][iok] * fmultiple) + fslope * (self.Orders[order]['iow'][iok] - np.average(self.Orders[order]['iow'][iok])), mode='same'), s=0)
    overflx = np.array([s.integral(x - self.Orders[order]['mindel']/2.0 + fshift, x + self.Orders[order]['mindel']/2.0 + fshift) for x in self.Orders[order]['wav'][ok]])
    return self.Orders[order]['wav'][ok], self.Orders[order]['flx'][ok], self.Orders[order]['err'][ok], self.Orders[order]['con'][ok], self.Orders[order]['pix'][ok], overflx
  
  # TODO save data + masks for each order
  # plot(rewritebinshiftandtilt(fitResults)[0],)
  def prettyResults(self):
    self.Results = {}
    # TODO add pixel info to this. 
    for binSizeKey in self.fitResults.keys():
      if binSizeKey == 'order':
        continue
      else:
        self.Results[binSizeKey] = {}
        for order in self.fitResults[binSizeKey].keys():
          self.Results[binSizeKey][order] = {}
          self.Results[binSizeKey][order]['avwav'] = []
          self.Results[binSizeKey][order]['cal'] = []
          self.Results[binSizeKey][order]['calerr'] = []
          self.Results[binSizeKey][order]['R'] = []
          self.Results[binSizeKey][order]['Rerr'] = []
          for bin in self.fitResults[binSizeKey][order].keys():
            if len(self.fitResults[binSizeKey][order][bin]) > 2:
              self.Results[binSizeKey][order]['avwav'].append(self.fitResults[binSizeKey][order][bin]['avwav'])
              self.Results[binSizeKey][order]['cal'].append(self.fitResults[binSizeKey][order][bin]['cal'])
              self.Results[binSizeKey][order]['calerr'].append(self.fitResults[binSizeKey][order][bin]['calerr'])
              self.Results[binSizeKey][order]['R'].append(self.fitResults[binSizeKey][order][bin]['R'])
              self.Results[binSizeKey][order]['Rerr'].append(np.average([self.fitResults[binSizeKey][order][bin]['Rbig'], self.fitResults[binSizeKey][order][bin]['Rsmall']]))
          shuffle = np.argsort(self.Results[binSizeKey][order]['avwav'])
          self.Results[binSizeKey][order]['avwav'] = np.array(self.Results[binSizeKey][order]['avwav'])[shuffle]
          self.Results[binSizeKey][order]['cal'] = np.array(self.Results[binSizeKey][order]['cal'])[shuffle]
          self.Results[binSizeKey][order]['calerr'] = np.array(self.Results[binSizeKey][order]['calerr'])[shuffle]
          self.Results[binSizeKey][order]['R'] = np.array(self.Results[binSizeKey][order]['R'])[shuffle]
          self.Results[binSizeKey][order]['Rerr'] = np.array(self.Results[binSizeKey][order]['Rerr'])[shuffle]
    pass

  # minor save
  def smallSave(self, filename='small.p'):
    """docstring for smallSave"""
    with open(filename, 'wb') as fp:
      pickle.dump(self.Results, fp)
    pass
  # full save
  def bigSave(self, filename='big.p'):
    """docstring for bigSave"""
    with open(filename, 'wb') as fp:
      pickle.dump(self.fitResults, fp)
    pass

  def saveFIT(self, filename="fit.fits"):
    """docstring for saveFIT"""
    with open(filename, 'wb') as fp:
      pickle.dump(self.Results, fp)
      pickle.dump(self.fitResults, fp)
      # TODO pickle.dump(self.exposureHeader, fp)
      # pickle.dump(nextDict, fp)
      # pickle.dump(bigDict, fp)
    pass

  def loadFIT(self, filename="fit.fits"):
    """docstring for loadFIT"""
    with open(filename, 'rb') as fp:
      self.loadResults = pickle.load(fp)
      self.loadfitResults = pickle.load(fp)
      # d1 = pickle.load(fp)
      # d2 = pickle.load(fp)
      # d3 = pickle.load(fp)
    pass
  # TODO 
  # scienceExposure = dict(HD138527.exposureHeader)
  # if scienceExposure['INSTRUME'] == 'HDS':
  #   print "HDS!"
  # calibrationExposure = dict()
  # def binshiftandtilt(self, order, bin, binSize, fmultiple, fshift, fsigma, elements, fslope, **kwargs):
  #   """Fit like shift with the addition of a slope across the order."""
  #   mask = self.Orders[order]['mask']
  #   kernel = self.gaussKernel(elements, fsigma)
  #   ok = reduce(np.logical_and, [self.Orders[order][binSize][bin]['ok'], mask])
  #   iok = self.Orders[order][binSize][bin]['iok']
  #   s = si.UnivariateSpline(self.Orders[order]['iow'][iok], np.convolve(kernel, (self.Orders[order]['iof'][iok] * fmultiple) + fslope * (self.Orders[order]['iow'][iok] - np.average(self.Orders[order]['iow'][iok])), mode='same'), s=0)
  #   overflx = np.array([s.integral(x - self.Orders[order]['mindel']/2.0 + fshift, x + self.Orders[order]['mindel']/2.0 + fshift) for x in self.Orders[order]['wav'][ok]])
  #   return np.sum( ((overflx - self.Orders[order]['flx'][ok] / self.Orders[order]['con'][ok]) / \
  #                   (self.Orders[order]['err'][ok] / self.Orders[order]['con'][ok])) ** 2 )
  
  def plotbinsinorder(self, order, binSize=350):
    """docstring for plotbinsinorder"""
    for bin in self.fitResults[350][order]:
      elements = self.fitResults[binSize][order][bin]['values']['elements']
      fmultiple = self.fitResults[binSize][order][bin]['values']['fmultiple']
      fshift = self.fitResults[binSize][order][bin]['values']['fshift']
      fsigma = self.fitResults[binSize][order][bin]['values']['fsigma']
      fslope = self.fitResults[binSize][order][bin]['values']['fslope']
      mask = self.Orders[order]['mask']
      kernel = self.gaussKernel(elements, fsigma)
      ok = reduce(np.logical_and, [self.Orders[order][binSize][bin]['ok'], mask])
      iok = self.Orders[order][binSize][bin]['iok']
      s = si.UnivariateSpline(self.Orders[order]['iow'][iok], np.convolve(kernel, (self.Orders[order]['iof'][iok] * fmultiple) + fslope * (self.Orders[order]['iow'][iok] - np.average(self.Orders[order]['iow'][iok])), mode='same'), s=0)
      overflx = np.array([s.integral(x - self.Orders[order]['mindel']/2.0 + fshift, x + self.Orders[order]['mindel']/2.0 + fshift) for x in self.Orders[order]['wav'][ok]])
      pl.plot(self.Orders[order]['wav'][ok], self.Orders[order]['flx'][ok]/self.Orders[order]['con'][ok], color="black", linewidth=2.0)
      pl.plot(self.Orders[order]['wav'][ok], overflx)
    pass
  
  def plotBinFit(self, order, bin, binSize, fmultiple, fshift, fsigma, elements, fslope, **kwargs):
    """docstring for plotBinFit"""
    mask = self.Orders[order]['mask']
    kernel = self.gaussKernel(elements, fsigma)
    ok = reduce(np.logical_and, [self.Orders[order][binSize][bin]['ok'], mask])
    iok = self.Orders[order][binSize][bin]['iok']
    s = si.UnivariateSpline(self.Orders[order]['iow'][iok], np.convolve(kernel, (self.Orders[order]['iof'][iok] * fmultiple) + fslope * (self.Orders[order]['iow'][iok] - np.average(self.Orders[order]['iow'][iok])), mode='same'), s=0)
    overflx = np.array([s.integral(x - self.Orders[order]['mindel']/2.0 + fshift, x + self.Orders[order]['mindel']/2.0 + fshift) for x in self.Orders[order]['wav'][ok]])
    pl.plot(self.Orders[order]['wav'], self.Orders[order]['flx'] / self.Orders[order]['con'], color="black", linewidth=2.0)
    pl.plot(np.average(self.Orders[order]['overwav'],axis=1), overflx)
    pass
    
  # def plotFitResults(self, order, fmultiple, fshift, fsigma, elements=1000, **kwargs):
  #   """docstring for plotFitResults"""
  #   kernel = self.gaussKernel(elements, fsigma)
  #   tck = si.splrep(self.Orders[order]['oiow'], np.convolve(kernel, self.Orders[order]['oiof'] * fmultiple, mode='same'))
  #   overflx = np.average(si.splev(np.hstack(self.Orders[order]['overwav']) + fshift, tck).reshape(np.shape(self.Orders[order]['overwav'])), axis=1)
  #   pl.plot(self.Orders[order]['wav'], self.Orders[order]['flx'] / self.Orders[order]['con'], color="black", linewidth=2.0)
  #   pl.plot(np.average(self.Orders[order]['overwav'],axis=1), overflx)    
  #   pass
  
def main(argv=None):
  pass

if __name__ == "__main__":
  main()
  
# Other ideas
# linear dispersion coefficient
# spectral line depth
# normalization: 
# 2nd-order dispersion coefficient
# width of main Gaussian IP
# width of box IP
# residuals
# plot kernel
# plot best fit between the two

# argsort for sorting bin wavelengths. # 
# d=np.arange(10)
# masks = [d>5, d % 2 == 0, d<8]
# you can use reduce to combine all of them:
# 
# total_mask = reduce(np.logical_and, masks)
# you can also use boolean operators explicitely if you need to manually choose the masks:
# 
# total_mask = masks[0] & masks[1] & masks[2]