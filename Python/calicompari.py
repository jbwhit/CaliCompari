#!/usr/bin/env python
# Copyright 2013 Jonathan Whitmore
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
import os
import glob
import barak
from barak.fitcont import spline_continuum
from barak import convolve
from barak import interp

import numpy as np
import scipy
import scipy as sp
# from scipy import optimize
import scipy.interpolate as si
import scipy.signal as ss
import scipy.constants as spc
import time
import datetime
import argparse
from ConfigParser import RawConfigParser, SafeConfigParser
import random as ra
import itertools
import pyfits as pf
import cPickle as pickle
import gzip
import matplotlib.pylab as pl
from matplotlib.backends.backend_pdf import PdfPages

pl.rcParams['figure.figsize'] = 16, 8  # that's default image size for this interactive session

import iminuit as mi

c_light = spc.c

help_message = '''
Various limitations:
Must have an FTS spectrum w/o gaps
Must have a telescope spectrum w/ monotonically increasing wavelength per order (gaps are OK)
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

def slope_to_array(slope, wavelength_array):
    """Return a wavelength array that is modified by a slope
    across the array."""
    import numpy as np
    midpoint = np.average([wavelength_array[0], wavelength_array[-1]])
    return slope * (wavelength_array - midpoint) + wavelength_array

def slope_to_array_old(slope, array):
    """Previous implementation."""
    import numpy as np
    return np.linspace(-slope, slope, num=len(array)) + array

def load_exposure(filename='big.gz'):
    """load_exposure recreates the final fit """
    with gzip.open(filename, 'rb') as file_handle:
        loadexposure = pickle.load(file_handle)
    return loadexposure

def save_exposure(object, filename='big.gz'):
    """Save to a large file all the information from the supercalibration
    fitting process. Includes data, and headers."""
    with gzip.open(filename, 'wb') as file_handle:
        pickle.dump(object, file_handle, pickle.HIGHEST_PROTOCOL)
    pass

def save_small_exposure(expo, filename="small.gz"):
    """docstring"""
    small_dictionary = {}
    small_dictionary['Results'] = expo.Results
    try:
        small_dictionary['arc_header'] = expo.arc_header
    except:
        small_dictionary['arc_header'] = ""
    try:
        small_dictionary['flux_header'] = expo.flux_header
    except:
        small_dictionary['flux_header'] = ""
    try:
        small_dictionary['header'] = expo.header
    except:
        small_dictionary['header'] = ""
    try:
        small_dictionary['exposure_file'] = expo.exposure_file
    except:
        small_dictionary['exposure_file'] = expo.exposureFile
    small_dictionary['safe_orders'] = expo.safe_orders
    with gzip.open(filename, 'wb') as file_handle:
        pickle.dump(small_dictionary, file_handle, pickle.HIGHEST_PROTOCOL)
    pass

def modify_small_exposure(small_dictionary, filename="small.gz"):
    """after modifying a small_dictionary file, save the resulting output."""
    with gzip.open(filename, 'wb') as file_handle:
        pickle.dump(small_dictionary, file_handle, pickle.HIGHEST_PROTOCOL)
    pass

def cleaned_data(filename_expo):
    """Modifies hand_tweak-ed file to contain only best data."""
    infile, expo = filename_expo
    upper_error_bound = expo["hand_tweak"]["upper_error_bound"]
    upper_wavelength_cutoff = expo["hand_tweak"]["upper_wavelength_cutoff"]
    badorders = expo["hand_tweak"]["badorders"]
    orderbegin = expo["hand_tweak"]["orderbegin"]
    orderend = expo["hand_tweak"]["orderend"]
    offset = expo["hand_tweak"]["offset"]
    minimum_number_of_chunks = expo["hand_tweak"]["minimum_number_of_chunks"]
    expo['cleaned'] = {}
    expo['cleaned'][500] = {}
    for order in [x for x in expo['safe_orders'] if x not in badorders]:
        if order in expo['Results'][500]:
            mask = (expo['Results'][500][order]['calerr'] < upper_error_bound) & (expo['Results'][500][order]['avwav'] < upper_wavelength_cutoff)
            expo['cleaned'][500][order] = {}
            if np.sum(mask) > minimum_number_of_chunks:
                expo['cleaned'][500][order]['avwav'] = expo['Results'][500][order]['avwav'][mask][orderbegin:orderend]
                expo['cleaned'][500][order]['cal'] = expo['Results'][500][order]['cal'][mask][orderbegin:orderend]
                expo['cleaned'][500][order]['calerr'] = expo['Results'][500][order]['calerr'][mask][orderbegin:orderend]
                expo['cleaned'][500][order]['avpix'] = expo['Results'][500][order]['avpix'][mask][orderbegin:orderend]
                expo['cleaned'][500][order]['R'] = expo['Results'][500][order]['R'][mask][orderbegin:orderend]
                expo['cleaned'][500][order]['Rerr'] = expo['Results'][500][order]['Rerr'][mask][orderbegin:orderend]
    pass

def hand_tweak( filename_expo,
                help=False,
                color="blue",
                linewidth=2.0,
                clobber=False,
                plot=True,
                vbuffer=800.0,
                *args,
                **kwargs):
    """Organized way of hand-tweaking the final calicompari results.
    The input is a already final fit with a dictionary that defines the modifications
    that allows for a reasonable fit."""
    tempx = []
    tempy = []
    infile, expo = filename_expo
    if ("hand_tweak" in expo.keys()) and (clobber == False):
        print "Turn clobber=True to modify this."
        return
    if "hand_tweak" not in expo.keys():
        print "Must first create hand_tweak dictionary. Some defaults."
        print """expo["hand_tweak"] = {}
        expo["hand_tweak"]["upper_error_bound"] = 200.0
        expo["hand_tweak"]["upper_wavelength_cutoff"] = 7600.0
        expo["hand_tweak"]["badorders"] = []
        expo["hand_tweak"]["orderbegin"] = 0
        expo["hand_tweak"]["orderend"] = -1
        expo["hand_tweak"]["offset"] = 0.
        expo["hand_tweak"]["minimum_number_of_chunks"] = 5"""

    upper_error_bound = expo["hand_tweak"]["upper_error_bound"]
    upper_wavelength_cutoff = expo["hand_tweak"]["upper_wavelength_cutoff"]
    badorders = expo["hand_tweak"]["badorders"]
    orderbegin = expo["hand_tweak"]["orderbegin"]
    orderend = expo["hand_tweak"]["orderend"]
    offset = expo["hand_tweak"]["offset"]
    minimum_number_of_chunks = expo["hand_tweak"]["minimum_number_of_chunks"]

    for order in [x for x in expo['safe_orders'] if x not in badorders]:
        if order in expo['Results'][500]:
            mask = (expo['Results'][500][order]['calerr'] < upper_error_bound) & (expo['Results'][500][order]['avwav'] < upper_wavelength_cutoff)
            if np.sum(mask) > minimum_number_of_chunks:
                if plot==True:
                    pl.errorbar(expo['Results'][500][order]['avwav'][mask][orderbegin:orderend],
                         expo['Results'][500][order]['cal'][mask][orderbegin:orderend] + offset,
                         expo['Results'][500][order]['calerr'][mask][orderbegin:orderend], color=color, linewidth=linewidth)
                tempx.append(np.average(expo['Results'][500][order]['avwav'][mask][orderbegin:orderend]))
                tempy.append(np.average(expo['Results'][500][order]['cal'][mask][orderbegin:orderend] + offset,
                    weights=1.0/(expo['Results'][500][order]['calerr'][mask][orderbegin:orderend])**2))
    tempx = np.hstack(tempx)
    tempy = np.hstack(tempy)
    # vbuffer = 800.0
    p, res, _, _, _ = np.polyfit(tempx, tempy, 1, full=True)
    wbuffer = 50.0
    pwav = np.arange(np.min(tempx) - wbuffer, np.max(tempx) + wbuffer)
    slope = p[0]
    intercept = p[1]
    expo["hand_tweak"]["calc_slope"] = slope
    expo["hand_tweak"]["calc_intercept"] = intercept
    if plot==True:
        pl.plot(pwav, slope * pwav + intercept, color="black", label="slope: " + str(round(slope * 1000,2)) + " m/s/1000 A")
        pl.ylim(np.average(tempy) - vbuffer, np.average(tempy) + vbuffer)
        pl.legend()
        pl.title(infile)
        pl.xlabel("Wavelength (Angstroms)", fontsize=20.0)
        pl.xticks(fontsize=20.0)
        pl.ylabel("v_shift (m/s)", fontsize=20.0)
        pl.yticks(fontsize=20.0)
    print "Filename: ", infile
    print "Slope: ", str(round(slope * 1000, 2)) + " m/s/1000 A"
    print "Current setup: "
    print " Discard data points with error larger than: ", upper_error_bound
    print " Exclude wavelenths greater than: ", upper_wavelength_cutoff
    print " Chunks each order between indices: ", orderbegin, orderend
    print " Removed orders: ", [order for order in badorders]
    print " Offset: ", offset
    print " Min # chunks required / order: ", minimum_number_of_chunks
    try:
        instrument = expo['flux_header']['INSTRUME']
        if instrument == "UVES":
            key = [key for key in expo['flux_header'].keys() if "WLEN" in str(key)][0]
            center_wavelength = expo['flux_header'][key]
            expo["hand_tweak"]["center_wavelength"] = center_wavelength
            expo["hand_tweak"]["center_offset"] = slope * center_wavelength + intercept
    except:
        pass
    if help == True:
        print """Some help:
        expo["hand_tweak"]["upper_error_bound"] = 200.0
        expo["hand_tweak"]["upper_wavelength_cutoff"] = 7600.0
        expo["hand_tweak"]["badorders"] = []
        expo["hand_tweak"]["orderbegin"] = 0
        expo["hand_tweak"]["orderend"] = -1
        expo["hand_tweak"]["offset"] = 0.
        expo["hand_tweak"]["minimum_number_of_chunks"] = 5"""
    pass

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
    def __init__(self, arcFile='', reduction_program='', calibration_type='', calibration_file='', exposure_file='', header_file=''):
        """docstring for __init__"""
        super(Exposure, self).__init__()
        self.arcFile = arcFile # a calibration Arc File
        self.exposure_file = exposure_file # a calibration Arc File
        self.reduction_program = reduction_program # reduction software used
        self.calibration_type = calibration_type # Calibration type: iodine, asteroid, none
        self.calibration_file = calibration_file # Calibration File

        self.header_file = header_file # science and arc file headers
        if self.header_file:
            with gzip.open(self.header_file, 'rb') as file_handle:
                loadheader = pickle.load(file_handle)
            self.arc_header, self.flux_header = loadheader[0], loadheader[1]

        self.fitGuess = AutoVivification()
        self.fit_starting = AutoVivification()
        self.fit_starting['initial'] = {}
        self.fit_starting['initial'].update({'shift':-0.003, 'fix_shift':False, 'limit_shift':(-1.5, 1.5), 'error_shift':0.03})
        self.fit_starting['initial'].update({'slope':-0.002, 'fix_slope':False, 'limit_slope':(-2.0, 2.0), 'error_slope':0.04})
        self.fit_starting['initial'].update({'sigma':3.102, 'fix_sigma':False, 'limit_sigma':(1.0, 10.0), 'error_sigma':0.2})
        self.fit_starting['initial'].update({'multiple':1.37, 'fix_multiple':False, 'limit_multiple':(0.1, 20.0), 'error_multiple':0.03})
        self.fit_starting['initial'].update({'offset':0.002, 'fix_offset':False, 'limit_offset':(-2.0, 2.0), 'error_offset':0.03})
        self.fit_starting['initial'].update({'minuit':0, 'fix_minuit':True})
        # self.fit_starting['initial'].update({'elements':100, 'fix_elements':True}) # UNDO THIS if minuit!
        # self.fit_starting['initial'].update({'set_strategy':2}) # UNDO THIS if minuit!
        self.fitResults = AutoVivification()
        if self.exposure_file.split('.')[-1] == 'fits':
            print "A fits exposure file."
            self.Orders = {}
            hdu = pf.open(self.exposure_file)
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
                    print "deleting", field
                    del(self.Orders[field])
        else:
            print "Not a fits file.", self.exposure_file
        pass

    def load_reference_spectra(self):
        """docstring for load_reference_spectra"""
        try:
            iow, iof = np.loadtxt(self.calibration_file)
        except:
            print "Consider saving a faster-loading calibration file."
            iow, iof = np.loadtxt(self.calibration_file, unpack='True')
        print "Reference FTS wavelength range:", iow[0], iow[-1]
        self.safe_orders = []
        for order in self.Orders:
            self.safe_orders.append(order)
            if (self.Orders[order]['wav'][0] > iow[0] + 40.0) & (self.Orders[order]['wav'][-1] < iow[-1] - 150.0):
                try:
                    ok = (self.Orders[order]['wav'][0] - 10 < iow) & (self.Orders[order]['wav'][-1] + 10 > iow)
                    if len(iow[ok]) > 200:
                        self.Orders[order]['iow'] = iow[ok].copy()
                        self.Orders[order]['iof'] = iof[ok].copy()
                except:
                    print "Order", order, "is outside overlap with reference FTS."
        for order in self.Orders:
            """removes orders"""
            try:
                len(self.Orders[order]['iow'])
                print order, "is safe."
            except:
                self.safe_orders.remove(order)
                print order, "was removed."
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
        for order in self.safe_orders:
            masks = []
            masks.append(self.Orders[order]['err'] > errorcutoff)
            masks.append(self.Orders[order]['flx'] > flxcutoff)
            masks.append(np.select([self.Orders[order]['err'] > 0],
                                   [np.nan_to_num(self.Orders[order]['flx']/self.Orders[order]['err']) >= sncutoff]))
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

    def continuum_fit(self, knots=10, plot=False, verbose=False):
        """fits a continuum via a spline through the flux values."""
        knots = 10
        edgeTolerance = 0.1
        remove_orders = []
        for order in self.safe_orders:
            mask = self.Orders[order]['mask']
            if np.sum(mask) < 100:
                remove_orders.append(order)
        for order in remove_orders:
            self.safe_orders.remove(order)
            print "Removing from safe_orders: ", order
        for order in self.safe_orders:
            mask = self.Orders[order]['mask']
            self.Orders[order]['con'] = np.zeros_like(self.Orders[order]['wav'])
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
        for order in self.safe_orders:
            self.Orders[order]['con'] = np.zeros_like(self.Orders[order]['wav'])
            mask = self.Orders[order]['mask']
            self.Orders[order]['con'][mask] = spline_continuum(self.Orders[order]['wav'][mask],
                                                    self.Orders[order]['flx'][mask],
                                                    self.Orders[order]['err'][mask],
                                                    np.linspace(self.Orders[order]['wav'][mask][0], self.Orders[order]['wav'][mask][-1], knots,),
                                                    nsig=nsig)[0]
        pass


    def oversample(self):
        """sets the minimum spacing in the telescope spectra (mindel)
        for each order over the whole exposure.
        Rename. """
        for order in self.safe_orders:
            mask = self.Orders[order]['mask']
            self.Orders[order]['mindel'] = self.Orders[order]['wav'][mask][-1] - self.Orders[order]['wav'][mask][0]
            for i in range(len(self.Orders[order]['wav'][mask]) - 1):
                adjacent_difference = self.Orders[order]['wav'][mask][i+1] - self.Orders[order]['wav'][mask][i]
                if self.Orders[order]['mindel'] > adjacent_difference:
                    self.Orders[order]['mindel'] = adjacent_difference
        pass


    def full_order_shift_scale(self, order=7, verbose=False, veryVerbose=False,
                                    robustSearch=False,
                                    first_guesses=None):
        """docstring for dictionaryShift
        first_guesses needs to look something like:
        first_guesses = {}
        first_guesses.update({'shift':-0.003,
                              'fix_shift':False,
                              'limit_shift':(-1.5, 1.5),
                              'error_shift':0.03})
        first_guesses.update({'slope':-0.002,
                              'fix_slope':False,
                              'limit_slope':(-2.0, 2.0),
                              'error_slope':0.04})
        first_guesses.update({'sigma':3.102,
                              'fix_sigma':False,
                              'limit_sigma':(1.0, 10.0),
                              'error_sigma':0.2})
        first_guesses.update({'multiple':1.37,
                              'fix_multiple':False,
                              'limit_multiple':(0.1, 20.0),
                              'error_multiple':0.03})
        first_guesses.update({'offset':0.002,
                              'fix_offset':False,
                              'limit_offset':(-2.0, 2.0),
                              'error_offset':0.03})
        first_guesses.update({'minuit':0, 'fix_minuit':True})
        """
        if first_guesses is None:
            first_guesses = self.fit_starting['initial']
        try:
            m = mi.Minuit(self.order_shift_and_scale_Akima, order=order, fix_order=True, **first_guesses)
            if veryVerbose==True:
                m.printMode=1
            if robustSearch==True:
                print "Robust search. Beginning initial scan..."
                m.scan(("fshift", 20, -0.5, 0.5))
                print "done."
            print "Finding initial full order shift/fit", '\n', datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            m.set_strategy(2)
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


    def order_shift_and_scale_Akima(self, order, multiple, shift, sigma, slope, offset, minuit, **kwargs):
        """trying to smooth, interpolate, and integrate the fit."""
        mask = self.Orders[order]['mask']
        iow = self.Orders[order]['iow']
        iof = self.Orders[order]['iof']
        wav = self.Orders[order]['wav'][mask]
        flx = self.Orders[order]['flx'][mask]
        err = self.Orders[order]['err'][mask]
        con = self.Orders[order]['con'][mask]
        pix = self.Orders[order]['pix'][mask]
        overflx = multiple * slope_to_array(slope, interp.interp_Akima(wav + shift, iow, convolve.convolve_constant_dv(iow, iof, vfwhm=sigma))) + offset
        chi_square = np.sum((overflx - flx/con)**2 / (err/con)**2)
        if minuit == 0:
            return chi_square
        else:
            return chi_square, wav, flx/con, err/con, pix, overflx

    def order_shift_and_scale_spline(self, order, multiple, shift, sigma, slope, offset, minuit, **kwargs):
        """trying to smooth, interpolate, and integrate the fit."""
        mask = self.Orders[order]['mask']
        iow = self.Orders[order]['iow']
        iof = self.Orders[order]['iof']
        wav = self.Orders[order]['wav'][mask]
        flx = self.Orders[order]['flx'][mask]
        err = self.Orders[order]['err'][mask]
        con = self.Orders[order]['con'][mask]
        pix = self.Orders[order]['pix'][mask]
        overflx = multiple * slope_to_array(slope, interp.interp_spline(wav + shift, iow, convolve.convolve_constant_dv(iow, iof, vfwhm=sigma))) + offset
        chi_square = np.sum((overflx - flx/con)**2 / (err/con)**2)
        if minuit == 0:
            return chi_square
        else:
            return chi_square, wav, flx/con, err/con, pix, overflx

    def create_bin_arrays(self, order=7, binSize=350, overlap=0.5, iowTolerance=2.0, minPixelsPerBin=100):
        """overlap is the fractional overlap or how much the bin is shifted relative to the binSize.
        so overlapping by .5 shifts by half binSize; .33 by .33 binSize. """
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

    def full_order_bin_shift_and_scale(self, order=7, binSize=350):
        self.fit_starting['order'][order] = self.fit_starting['initial']
        self.fit_starting['order'][order].update(self.fitResults['order'][order]['values'])
        for singlebin in self.Orders[order][binSize]:
            self.fitResults[binSize][order][singlebin] = {}
            self.small_bin_shift(order, binSize, singlebin)
        pass

    def small_bin_shift(self, order=7, binSize=350, singlebin=2, veryVerbose=False, robustSearch=False):
        """docstring for smallBinShift"""
        # TODO check that the full order solution has run.
        try:
            type(self.fitResults['order'][order]['values'])
        except:
            print "It doesn't look like the full order was run... "
        m = mi.Minuit(self.bin_shift_and_tilt_Akima, order=order, binSize=binSize, singlebin=singlebin, fix_order=True, fix_binSize=True, fix_singlebin=True, **self.fit_starting['order'][order])
        if veryVerbose==True:
            m.printMode=1
        if robustSearch==True:
            print "Robust search. Beginning initial scan..."
            m.scan(("fshift", 20, -0.5, 0.5))
            print "done."
        try:
            print datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "Finding initial shift/fit for order:", order, "and bin:", singlebin
            m.set_strategy(2)
            m.migrad()
            self.fitResults[binSize][order][singlebin]['values'] = m.values
            self.fitResults[binSize][order][singlebin]['errors'] = m.errors
            mask = self.Orders[order]['mask']
            ok = reduce(np.logical_and, [self.Orders[order][binSize][singlebin]['ok'], mask])
            iok = self.Orders[order][binSize][singlebin]['iok']
            wav = self.Orders[order]['wav'][ok]
            pix = self.Orders[order]['pix'][ok]
            lamb = np.average(wav)
            avpix = np.average(pix)
            cal = m.values['shift'] * c_light / lamb
            calerr = m.errors['shift'] * c_light / lamb
            R = c_light / m.values['sigma'] / 1000.
            Rerr = c_light / m.errors['sigma'] / 1000.
            self.fitResults[binSize][order][singlebin]['avwav'] = lamb
            self.fitResults[binSize][order][singlebin]['cal'] = cal
            self.fitResults[binSize][order][singlebin]['calerr'] = calerr
            self.fitResults[binSize][order][singlebin]['R'] = R
            self.fitResults[binSize][order][singlebin]['Rerr'] = Rerr
            self.fitResults[binSize][order][singlebin]['avpix'] = avpix
            self.fitResults[binSize][order][singlebin]['converged'] = True
            self.fitResults[binSize][order][singlebin]['values']['minuit'] = 1
            chisq, wav, nflx, nerr, pix, overflx = self.bin_shift_and_tilt_Akima(**self.fitResults[binSize][order][singlebin]['values'])
            self.fitResults[binSize][order][singlebin]['chisq'] = chisq
            self.fitResults[binSize][order][singlebin]['wav'] = wav
            self.fitResults[binSize][order][singlebin]['nflx'] = nflx
            self.fitResults[binSize][order][singlebin]['nerr'] = nerr
            self.fitResults[binSize][order][singlebin]['pix'] = pix
            self.fitResults[binSize][order][singlebin]['overflx'] = overflx
            print datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "finished."
        except:
            self.fitResults[binSize][order][singlebin]['converged'] = False
            print "Serious problem with bin:", singlebin
        pass

    def bin_shift_and_tilt_Akima(self, order, singlebin, binSize, multiple, shift, sigma, slope, offset, minuit, **kwargs):
        """trying to smooth, interpolate, and integrate the fit."""
        mask = self.Orders[order]['mask']
        ok = reduce(np.logical_and, [self.Orders[order][binSize][singlebin]['ok'], mask])
        iok = self.Orders[order][binSize][singlebin]['iok']
        iow = self.Orders[order]['iow'][iok]
        iof = self.Orders[order]['iof'][iok]
        wav = self.Orders[order]['wav'][ok]
        flx = self.Orders[order]['flx'][ok]
        err = self.Orders[order]['err'][ok]
        con = self.Orders[order]['con'][ok]
        pix = self.Orders[order]['pix'][ok]
        overflx = multiple * slope_to_array(slope, interp.interp_Akima(wav + shift, iow, convolve.convolve_constant_dv(iow, iof, vfwhm=sigma))) + offset
        chi_square = np.sum((overflx - flx/con)**2 / (err/con)**2)
        if minuit == 0:
            return chi_square
        else:
            return chi_square, wav, flx/con, err/con, pix, overflx

    def bin_shift_and_tilt_spline(self, order, singlebin, binSize, multiple, shift, sigma, slope, offset, minuit, **kwargs):
        """trying to smooth, interpolate, and integrate the fit."""
        mask = self.Orders[order]['mask']
        ok = reduce(np.logical_and, [self.Orders[order][binSize][singlebin]['ok'], mask])
        iok = self.Orders[order][binSize][singlebin]['iok']
        iow = self.Orders[order]['iow'][iok]
        iof = self.Orders[order]['iof'][iok]
        wav = self.Orders[order]['wav'][ok]
        flx = self.Orders[order]['flx'][ok]
        err = self.Orders[order]['err'][ok]
        con = self.Orders[order]['con'][ok]
        pix = self.Orders[order]['pix'][ok]
        overflx = multiple * slope_to_array(slope, interp.interp_spline(wav + shift, iow, convolve.convolve_constant_dv(iow, iof, vfwhm=sigma))) + offset
        chi_square = np.sum((overflx - flx/con)**2 / (err/con)**2)
        if minuit == 0:
            return chi_square
        else:
            return chi_square, wav, flx/con, err/con, pix, overflx

    def make_pretty_results(self):
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
                            self.Results[binSizeKey][order]['Rerr'].append(self.fitResults[binSizeKey][order][bin]['Rerr'])
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

