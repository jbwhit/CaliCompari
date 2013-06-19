import sys
# import shutil
import os
import glob
import barak
from barak.fitcont import spline_continuum
from barak import convolve
from barak import interp

# import pprint 
import numpy as np
import scipy
import scipy as sp
from scipy import optimize
import scipy.interpolate as si
import scipy.signal as ss
import scipy.constants as spc
# import subprocess as subprocess
import time
import datetime
# from optparse import OptionParser
import argparse
from ConfigParser import RawConfigParser, SafeConfigParser
import random as ra
# import tempfile
import itertools
# import pl as pl
import pyfits as pf
import cPickle as pickle
import matplotlib.pylab as pl
from matplotlib.backends.backend_pdf import PdfPages

pl.rcParams['figure.figsize'] = 16, 8  # that's default image size for this interactive session

import minuit as mi 
# import iminuit as mi 

c_light = spc.c
              
help_message = '''
Various limitations: 
Must have an FTS spectrum w/o gaps
Must have a telescope spectrum w/ monotonically increasing wavelength (gaps are OK)
The spacing of the nearest two pixels in the telescope spectrum is used as the pixel size for each order.
'''

class AutoVivification(dict):
    """Implementation of perl's autovivification feature.
    from: http://stackoverflow.com/questions/651794/whats-the-best-way-to-initialize-a-dict-of-dicts-in-python
    Testing:

    >>> a = AutoVivification()
    >>> a[1][2][3] = 4
    >>> a[1][3][3] = 5
    >>> a[1][2]['test'] = 6
    >>> print a
    Output:

    {1: {2: {'test': 6, 3: 4}, 3: {3: 5}}}
    """
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

#
def slope_to_array(slope, array):
    import numpy as np
    return np.linspace(-slope, slope, num=len(array)) + array
    
    
class Exposure(object):
    """An oject class that contains the data for quasar absorption spectroscopy study.
    
    An exposure has: 
        orders which have
            pixels with corresponding
            wavelength
            flux
            error values. 
            
        fit.
        """
    def __init__(self, arcFile='', reductionProgram='', calibrationType='', calibrationFile='', exposureFile=''):
        """docstring for __init__"""
        super(Exposure, self).__init__()
        self.arcFile = arcFile # a calibration Arc File
        self.exposureFile = exposureFile # a calibration Arc File
        self.reductionProgram = reductionProgram # reduction software used
        self.calibrationType = calibrationType # Calibration type: iodine, asteroid, none
        self.calibrationFile = calibrationFile # Calibration File
        self.fitGuess = AutoVivification()
        self.fit_starting = AutoVivification()
        self.fit_starting['initial'] = {}
        self.fit_starting['initial'].update({'shift':-0.3, 'fix_shift':False, 'limit_shift':(-1.0,1.0), 'err_shift':0.05})
        self.fit_starting['initial'].update({'slope':-0.002, 'fix_slope':False, 'limit_slope':(-2.0,2.0), 'err_slope':0.1})
        self.fit_starting['initial'].update({'sigma':8.102, 'fix_sigma':False, 'limit_sigma':(1.0,25.0), 'err_sigma':0.5})
        self.fit_starting['initial'].update({'multiple':1.37, 'fix_multiple':False, 'limit_multiple':(0.1,20.0), 'err_multiple':0.1})
        self.fit_starting['initial'].update({'offset':0.002, 'fix_offset':False, 'limit_offset':(-2.0,2.0), 'err_offset':0.3})
        self.fit_starting['initial'].update({'elements':100, 'fix_elements':True})
        self.fit_starting['initial'].update({'minuit':0, 'fix_minuit':True})
        self.fit_starting['initial'].update({'strategy':2})
        # self.fit_starting['initial'].update({'offset':0.002, 'fix_offset':False, 'limit_offset':(-1.0,1.0)})
        # Basic strategy: convolve iof; then, shift/tilt/multiply/and offset it. 
        self.fitGuess['initial'] = { 'fshift':0.002, 'fix_fshift':False, 'limit_fshift':(-1.0,1.0) ,'err_fshift':0.005 }
        self.fitGuess['initial'].update({ 'fsigma':10.5, 'fix_fsigma':False, 'limit_fsigma':(2.0,30.0) ,'err_fsigma':0.05 })
        self.fitGuess['initial'].update({ 'fmultiple':50.25, \
                                          'fix_fmultiple':False, \
                                          'limit_fmultiple':(0.1, 100.0), \
                                          'err_fmultiple':0.02 })
        self.fitGuess['initial'].update({ 'fslope':0.0005, 'fix_fslope':False, 'limit_fslope':(-3.0, 3.0) ,'err_fslope':0.05 })
        self.fitGuess['initial'].update({ 'elements':100, 'fix_elements':True })
        # self.fitGuess['initial'].update({ 'fwidth':200, 'fix_fwidth':True })
        self.fitGuess['initial'].update({ 'strategy':2 })
        # TODO check if fit hits against any of the imposed limits. 
        self.fitResults = AutoVivification()
        if self.exposureFile.split('.')[-1] == 'fits':
            print "A fits exposure file."
            self.Orders = {}
            hdu = pf.open(self.exposureFile)
            self.header = hdu[0].header
            for index, table in enumerate(hdu):
                try:
                    type(hdu[index].data)
                    self.Orders[index] = {}
                    self.Orders[index]['wav'] = table.data[0]
                    self.Orders[index]['flx'] = table.data[1]
                    self.Orders[index]['err'] = table.data[2]
                    self.Orders[index]['pix'] = np.arange(len(self.Orders[index]['wav']))
                except:
                    self.exposureHeader = hdu[-1].header
            for field in self.Orders.keys():
                if len(self.Orders[field]) < 1:
                    del(self.Orders[field])
        else:
            print "Not a fits file.", self.exposureFile
        pass

    def loadReferenceSpectra(self):
        """docstring for loadReferenceSpectra"""
        try: 
            iow, iof = np.loadtxt(self.calibrationFile)
        except:
            print "Consider saving a faster-loading calibration file."
            iow, iof = np.loadtxt(self.calibrationFile, unpack='True')
        print "Reference FTS wavelength range:", iow[0], iow[-1]
        for order in self.Orders:
            if (self.Orders[order]['wav'][0] > iow[0] + 40.0) & (self.Orders[order]['wav'][-1] < iow[-1] - 150.0):
                try:
                    ok = (self.Orders[order]['wav'][0] - 10 < iow) & (self.Orders[order]['wav'][-1] + 10 > iow)
                    if len(iow[ok]) > 200:
                        self.Orders[order]['iow'] = iow[ok].copy()
                        self.Orders[order]['iof'] = iof[ok].copy()
                except:
                    print "Order", order, "is outside overlap with reference FTS."
        pass
    
    def cleanup(self, verbose=False):
        """mask out bad regions of the spectra
        Example config file setup. 
        [skylines]
        remove: 
            5589.128        5589.132
            5865.454        5865.459
        """
        parser = SafeConfigParser()
        candidates = glob.glob('config*')
        found = parser.read(candidates)
        try:
            wavekill = parser.get('skylines','remove')
        except:
            print "Warning: not removing skylines (if you want this create a config.wavekill file)."
        if verbose==True:
            print "Beginning cleanup of data...", datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
        errorcutoff = 0.0
        flxcutoff = 0.0
        sncutoff = 10.0
        for order in self.Orders:
            masks = []
            masks.append(self.Orders[order]['err'] > errorcutoff)
            masks.append(self.Orders[order]['flx'] > flxcutoff)
            masks.append(np.select([self.Orders[order]['err'] > 0], [self.Orders[order]['flx']/self.Orders[order]['err'] >= sncutoff]))
            try:
                for killLine in wavekill.splitlines():
                    if len(killLine) > 1:
                        masks.append(reduce(np.logical_or, \
                                    [self.Orders[order]['wav'] < float(killLine.split()[0]), \
                                        self.Orders[order]['wav'] > float(killLine.split()[1])]))
            except:
                pass
            self.Orders[order]['mask'] = reduce(np.logical_and, masks)
        pass
    
    def continuumFit(self, knots=10, plot=False, verbose=False):
        """fits a continuum via a spline through the flux values."""
        knots = 10
        edgeTolerance = 0.1
        for order in self.Orders:
            mask = self.Orders[order]['mask']
            self.Orders[order]['con'] = np.zeros(len(self.Orders[order]['wav']))
            s = si.LSQUnivariateSpline(self.Orders[order]['wav'][mask],\
                                                                self.Orders[order]['flx'][mask],\
                                                                np.linspace(self.Orders[order]['wav'][mask][0]+edgeTolerance,\
                                                                self.Orders[order]['wav'][mask][-1]-edgeTolerance, knots),\
                                                                w=self.Orders[order]['err'][mask])
            self.Orders[order]['con'][mask] = s(self.Orders[order]['wav'][mask]) # new array is made -- continuum
        pass
    
    def continuum_fit_2(self, knots=10, nsig=4.0):
        """barak implementation"""
        knots = knots
        for order in self.Orders:
            self.Orders[order]['con'] = np.zeros_like(self.Orders[order]['wav'])
            mask = self.Orders[order]['mask']
            self.Orders[order]['con'][mask] = spline_continuum(self.Orders[order]['wav'][mask], 
                                                    self.Orders[order]['flx'][mask], 
                                                    self.Orders[order]['err'][mask], 
                                                    np.linspace(self.Orders[order]['wav'][mask][0], self.Orders[order]['wav'][mask][-1], knots,), 
                                                    nsig=nsig)[0]
        pass
        
    
    def OverSample(self):
        """sets the minimum spacing in the telescope spectra (mindel) 
        for each order over the whole exposure.
        Rename. """
        for order in self.Orders:
            mask = self.Orders[order]['mask']
            self.Orders[order]['mindel'] = self.Orders[order]['wav'][mask][-1] - self.Orders[order]['wav'][mask][0]
            for i in range(len(self.Orders[order]['wav'][mask]) - 1):
                adjacent_difference = self.Orders[order]['wav'][mask][i+1] - self.Orders[order]['wav'][mask][i]
                if self.Orders[order]['mindel'] > adjacent_difference: 
                    self.Orders[order]['mindel'] = adjacent_difference
        pass
    
    def fullExposureShift(self, verbose=False, veryVerbose=False, robustSearch=False, binSize=350):
        """docstring for fullExposureShift"""
        starttime=datetime.datetime.now()
        #.strftime("%Y-%m-%d %H:%M:%S")
        for order in self.Orders:
            if 'iow' in self.Orders[order]:
                print "Working on order: ", order
                self.CreateBinArrays(order=order, binSize=binSize) # new!
                try:
                    self.OrderShiftandTilt(order=order, veryVerbose=veryVerbose) # new!
                    self.fullOrderBinShift(order=order, binSize=binSize)
                except:
                    print "Order or bin failed."
        print "Took:", str(datetime.datetime.now() - starttime), " to finish exposure."
        pass

    def OrderShiftandTilt(self, order=7, verbose=False, veryVerbose=False, robustSearch=False):
        """docstring for dictionaryShift"""
        try:
            m = mi.Minuit(self.shiftandtilt, order=order, fix_order=True, **self.fitGuess['initial'])
            if veryVerbose==True:
                m.printMode=1
            if robustSearch==True:
                print "Robust search. Beginning initial scan..."
                m.scan(("fshift", 20, -0.5, 0.5))
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

    def full_order_shift_scale(self, order=7, verbose=False, veryVerbose=False, robustSearch=False):
        """docstring for dictionaryShift"""
        try:
            m = mi.Minuit(self.order_shift_and_scale, order=order, fix_order=True, **self.fit_starting['initial'])
            if veryVerbose==True:
                m.printMode=1
            if robustSearch==True:
                print "Robust search. Beginning initial scan..."
                m.scan(("fshift", 20, -0.5, 0.5))
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
        return np.sum( ((overflx - self.Orders[order]['flx'][mask] / self.Orders[order]['con'][mask])    / \
                                        (self.Orders[order]['err'][mask] / self.Orders[order]['con'][mask])) ** 2)

    def order_shift_and_scale(self, order, multiple, shift, sigma, slope, offset, minuit, **kwargs):
        """trying to smooth, interpolate, and integrate the fit."""
        mask = self.Orders[order]['mask']
        # ok = reduce(np.logical_and, [self.Orders[order][binSize][bin]['ok'], mask])
        # iok = self.Orders[order][binSize][bin]['iok']
        iow = self.Orders[order]['iow']
        iof = self.Orders[order]['iof']
        wav = self.Orders[order]['wav'][mask]
        flx = self.Orders[order]['flx'][mask]
        err = self.Orders[order]['err'][mask]
        con = self.Orders[order]['con'][mask]
        pix = self.Orders[order]['pix'][mask]
        # TODO
        #print slope, len(wav), len(iow), len(iof), sigma, offset
        overflx = multiple * slope_to_array(slope, interp.interp_Akima(wav + shift, iow, convolve.convolve_constant_dv(iow, iof, vfwhm=sigma))) + offset        
        chi_square = np.sum((overflx - flx/con)**2 / (err/con)**2 ) / len(wav)
        print chi_square, slope, sigma, offset, shift, multiple
        if minuit == 0:
            return chi_square
        else:
            return chi_square, wav, flx/con, err/con, pix, overflx

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
        self.fitGuess['order'][order] = self.fitGuess['initial']
        self.fitGuess['order'][order].update(self.fitResults['order'][order]['values'])
        self.fitGuess['order'][order].update({ 'elements':int(10.0 * self.fitResults['order'][order]['values']['fsigma']) })
        print self.fitResults['order'][order]
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
        m = mi.Minuit(self.bin_shift_and_tilt, order=order, binSize=binSize, bin=bin, fix_order=True, fix_binSize=True, fix_bin=True, chi2=1, fix_chi2=True, **self.fitGuess['order'][order])
        # m = mi.Minuit(self.binChiSquare, order=order, binSize=binSize, bin=bin, fix_order=True, fix_binSize=True, fix_bin=True, **self.fitGuess['order'][order])
        if veryVerbose==True:
            m.printMode=1
        if robustSearch==True:
            print "Robust search. Beginning initial scan..."
            m.scan(("fshift",20,-0.5,0.5))
            print "done."
        try: 
            print datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "Finding initial shift/fit for order:", order, "and bin:", bin
            m.migrad()
            # m.minos()
            self.fitResults[binSize][order][bin]['values'] = m.values
            self.fitResults[binSize][order][bin]['errors'] = m.errors # todo m.merrors 
            mask = self.Orders[order]['mask']
            ok = reduce(np.logical_and, [self.Orders[order][binSize][bin]['ok'], mask])
            iok = self.Orders[order][binSize][bin]['iok']
            elements = self.fitResults[binSize][order][bin]['values']['elements']
            lamb = np.average(self.Orders[order]['wav'][ok]) # 
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
            self.fitResults[binSize][order][bin]['converged'] = True
            print datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "finished."
        except:
            self.fitResults[binSize][order][bin]['converged'] = False
            print "Serious problem with bin:", bin
        pass
        
    # TODO offset; shift; tilt; multiple; sigma
    def bin_shift_and_tilt(self, order, bin, binSize, fmultiple, fshift, fsigma, elements, fslope, chi2=1, **kwargs):
        """Fit like shift with the addition of a slope across the order."""
        mask = self.Orders[order]['mask']
        ok = reduce(np.logical_and, [self.Orders[order][binSize][bin]['ok'], mask])
        iok = self.Orders[order][binSize][bin]['iok']
        iow = self.Orders[order]['iow'][iok]
        iof = self.Orders[order]['iof'][iok]
        wav = self.Orders[order]['wav'][ok]
        flx = self.Orders[order]['flx'][ok]
        err = self.Orders[order]['err'][ok]
        con = self.Orders[order]['con'][ok]
        pix = self.Orders[order]['pix'][ok]
        half_pixel = self.Orders[order]['mindel']/2.0
        # s = si.UnivariateSpline(iow, np.convolve(kernel, (iof * fmultiple) + fslope * (iow - np.average(iow)), mode='same'), s=0)
        # spline = interp.AkimaSpline(iow, iof)
        # hodor = interp.interp_Akima(wavdict[0], iow, iof)
        # overflx = np.array([s.integral(x - half_pixel + fshift, x + half_pixel + fshift) for x in wav])
        # overflx = interp.interp_Akima(wav + fshift, iow, np.convolve(kernel, iof * fmultiple + fslope*(iow-np.average(iow))))
        convolved = convolve.convolve_constant_dv(iow, iof * fmultiple + fslope * (iow - np.average(iow)), vfwhm=fsigma)
        overflx = interp.interp_Akima(wav + fshift, iow, convolved)
        if chi2 == 1:
            return np.sum( ((overflx - flx / con) / (err / con)) ** 2 )
        else:
            return np.sum( ((overflx - flx / con) / (err / con)) ** 2 ), wav, flx/con, err/con, pix, overflx        
    
    # def bin_shift_and_tilt(self, order, bin, binSize, fmultiple, fshift, fsigma, elements, fslope, chi2=1, **kwargs):
    #     """Fit like shift with the addition of a slope across the order."""
    #     mask = self.Orders[order]['mask']
    #     kernel = self.gaussKernel(elements, fsigma)
    #     ok = reduce(np.logical_and, [self.Orders[order][binSize][bin]['ok'], mask])
    #     iok = self.Orders[order][binSize][bin]['iok']
    #     iow = self.Orders[order]['iow'][iok]
    #     iof = self.Orders[order]['iof'][iok]
    #     wav = self.Orders[order]['wav'][ok]
    #     flx = self.Orders[order]['flx'][ok]
    #     err = self.Orders[order]['err'][ok]
    #     con = self.Orders[order]['con'][ok]
    #     pix = self.Orders[order]['pix'][ok]
    #     half_pixel = self.Orders[order]['mindel']/2.0
    #     s = si.UnivariateSpline(iow, np.convolve(kernel, (iof * fmultiple) + fslope * (iow - np.average(iow)), mode='same'), s=0)
    #     overflx = np.array([s.integral(x - half_pixel + fshift, x + half_pixel + fshift) for x in wav])
    #     if chi2 == 1:
    #         return np.sum( ((overflx - flx / con) / (err / con)) ** 2 )
    #     else:
    #         return np.sum( ((overflx - flx / con) / (err / con)) ** 2 ), wav, flx/con, err/con, pix, overflx        
    # 
    # TODO save data + masks for each order
    # plot(rewritebin_shift_and_tilt(fitResults)[0],)
    def prettyResults(self):
        self.Results = {}
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
                    self.Results[binSizeKey][order]['avpix'] = []
                    self.Results[binSizeKey][order]['converged'] = []
                    for bin in self.fitResults[binSizeKey][order].keys():
                        if len(self.fitResults[binSizeKey][order][bin]) > 2:
                            self.Results[binSizeKey][order]['avwav'].append(self.fitResults[binSizeKey][order][bin]['avwav'])
                            self.Results[binSizeKey][order]['cal'].append(self.fitResults[binSizeKey][order][bin]['cal'])
                            self.Results[binSizeKey][order]['calerr'].append(self.fitResults[binSizeKey][order][bin]['calerr'])
                            self.Results[binSizeKey][order]['R'].append(self.fitResults[binSizeKey][order][bin]['R'])
                            self.Results[binSizeKey][order]['Rerr'].append(np.average([self.fitResults[binSizeKey][order][bin]['Rbig'], self.fitResults[binSizeKey][order][bin]['Rsmall']]))
                            self.Results[binSizeKey][order]['avpix'].append(self.fitResults[binSizeKey][order][bin]['avpix'])
                            self.Results[binSizeKey][order]['converged'].append(self.fitResults[binSizeKey][order][bin]['converged'])
                    shuffle = np.argsort(self.Results[binSizeKey][order]['avwav'])
                    self.Results[binSizeKey][order]['avwav'] = np.array(self.Results[binSizeKey][order]['avwav'])[shuffle]
                    self.Results[binSizeKey][order]['cal'] = np.array(self.Results[binSizeKey][order]['cal'])[shuffle]
                    self.Results[binSizeKey][order]['calerr'] = np.array(self.Results[binSizeKey][order]['calerr'])[shuffle]
                    self.Results[binSizeKey][order]['R'] = np.array(self.Results[binSizeKey][order]['R'])[shuffle]
                    self.Results[binSizeKey][order]['Rerr'] = np.array(self.Results[binSizeKey][order]['Rerr'])[shuffle]
                    self.Results[binSizeKey][order]['avpix'] = np.array(self.Results[binSizeKey][order]['avpix'])[shuffle]
                    self.Results[binSizeKey][order]['converged'] = np.array(self.Results[binSizeKey][order]['converged'])[shuffle]
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
        pass

    # TODO
    # scienceExposure = dict(HD138527.exposureHeader)
    # if scienceExposure['INSTRUME'] == 'HDS':
    #     print "HDS!"
    # calibrationExposure = dict()

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

