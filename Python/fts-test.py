#!/usr/bin/env python
# encoding: utf-8
"""
fts-test.py

Created by Jonathan Whitmore on 2010-05-26.
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
# import matplotlib
# matplotlib.use('PDF')
import matplotlib.pylab as pl
import sys
import os
import glob
import datetime
from optparse import OptionParser
from ConfigParser import RawConfigParser



filename='config'
config = RawConfigParser()
config.read(filename)
telescope = config.get('Information','telescope')
astro_object = "FTS"
spectrograph = config.get('Information','spectrograph')
# bin_size = float(config.get('Parameters','bin_size'))  # 300.0 # in km/s
# step = float(config.get('Parameters','step_size')) #  50.0 # km/s
bin_size = 350.0
step = 100.0
gauss_elements = 60.0
sigma = float(config.get('Parameters','sigma'))

if spectrograph=="UVES": 
    print "Instrument: ", spectrograph
    FTS_FILE = config.get('Information','path_to_UVES_FTS')
if spectrograph=="HIRES":
    print "Instrument: ", spectrograph
    FTS_FILE = config.get('Information','path_to_HIRES_FTS')

print "Telescope: ", telescope
print "Object: ", astro_object
print "FTS File: ", FTS_FILE
print "Bin size: ", bin_size, "km/s", "Step size: ", step, "km/s"
print "# of Elements: ", gauss_elements, "Initial sigma: ", sigma

# =================
# = Parse Options =
# =================

core_test = 0

parser = OptionParser()
parser.add_option("-e","--exposure", 
                  action="store", type="string", dest="exponame", 
                  help="name exposures", metavar="EXPO")
parser.add_option("-c","--core", 
                action="store", type="string", dest="corename", 
                help="Dual Core Advantage", metavar="CORE")

(options, args) = parser.parse_args()

# continuum_data = str(options.exponame) 

# Optional use of dual core stuff
try: core_test += float(options.corename) 
except: pass

# ======================
# = Set Data Variables =
# ======================

# Read in exposures that have iodine exposures
i2exposures = ['01'] 


# Run with flag -c 1    or --core 1
core1 = i2exposures[:len(i2exposures)/2]
core2 = i2exposures[len(i2exposures)/2:]

if core_test == 1:
    exposures_to_analyze = core1
    print "Exposures: ", exposures_to_analyze
if core_test == 2:
    exposures_to_analyze = core2
    print "Exposures: ", exposures_to_analyze
if core_test != 1 and core_test !=2: 
    exposures_to_analyze = core1 + core2
    print "Exposures: ", exposures_to_analyze


# =============
# = Constants =
# =============
c_light = 299792458.

# ====================
# = Read in FTS Data =
# ====================
iodine_reader = csv.reader(open("/Users/jonathan/Code/Spectrograph/FTS/vlt.20091012.ascii"), delimiter=' ')
iow = []
iof = []
for row in iodine_reader:
    iow.append(float(row[0]))
    iof.append(float(row[1]))

iow = np.array(iow)
iof = np.array(iof)

iodine2_reader = csv.reader(open("/Users/jonathan/Code/Spectrograph/FTS/hires.20100526.ascii"), delimiter=' ')
iow2 = []
iof2 = []
for row in iodine2_reader:
    iow2.append(float(row[0]))
    iof2.append(float(row[1]))

iow2 = np.array(iow2)
iof2 = np.array(iof2)


    
for core_exposure in i2exposures:
    continuum_data = str("Continuum/" + astro_object + "." + core_exposure + "*")
    DATA_LIST = glob.glob("Continuum/HR1996.01.l.01*")

    today = datetime.date.today()
    logfile = open("Logs/log." + str(today) + ".ascii", 'w')

    for data_file in DATA_LIST:
        starting = 5015. 
        ending = 6015. 
        wav = iow2[np.where(iow2[np.where(iow2 < ending)] > starting)]
        flx = iof2[np.where(iow2[np.where(iow2 < ending)] > starting)]
        
        # ================================================
        # = Logic Test to see if file overlaps FTS data =
        # ================================================
        if(wav[0] < iow[0]+10 or wav[-1] > iow[-1] - 10):
            print data_file, "is outside of overlap"
            continue
        
        # =================================
        # = Slice FTS to manageable size  =
        # =================================
        slice_iow = iow[np.where(iow > wav[0] - 10)]
        slice_iow = slice_iow[np.where(slice_iow < wav[-1] + 10)]
        slice_iof = iof[np.where(iow > wav[0] - 10)]
        slice_iof = slice_iof[np.where(slice_iow < wav[-1] + 10)]
    
        starting_flx = flx
    
        # ======================================
        # = Initial Shift; Kernel Size finding =
        # ======================================
        def  initial_shift(fmultiple,fshift,fsigma):
            better_flx = starting_flx * fmultiple
            better_g = ss.gaussian(gauss_elements,fsigma)
            better_kernel = better_g / sum(better_g)
            better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                                   slice_iof, mode='same'))
            better_y2 = si.splev(wav+fshift,better_tck)
            return sum((better_y2 - better_flx) ** 2 / \
                       (fmultiple ) ** 2 )
    
        m = mi.Minuit(initial_shift,fmultiple=0.82,fshift=0.01,\
                             fsigma=15.,strategy=2)

        # m.printMode=1
        m.migrad()
        #m.minos()
        order_velocity = c_light / wav[len(wav)/2]
        print "\nShift: ", m.values["fshift"] * order_velocity
        print "Sigma: ", m.values["fsigma"], m.errors["fsigma"]
        
        # Warning if constraining fitting function
        if gauss_elements/m.values["fsigma"] < 7.5:
            print "Use a larger gaussian window\nYour current array is: ", gauss_elements
            print "Your current sigma is: ", m.values["fsigma"]
        
        print >>logfile, wav[len(wav)/2], m.values["fshift"] * order_velocity, \
                m.errors["fshift"] * order_velocity, data_file
    
        # TODO Print a log that makes sense 
        
        # ========================================
        # = Create the bins and structure needed =
        # ========================================
    
        # Velocity from mid-wavelength of the order 
        whit_index = []
        whit_step = wav[0] * step * 1000.0 / c_light
        for i in np.arange(wav[0],wav[-1], whit_step ):
            temp = i + i * bin_size * 1000.0 / c_light
            if len(np.where(wav[np.where(wav < temp)] > i )[0]) > 100 :
                whit_index.append(np.where(wav[np.where(wav < temp)[0]] > i ))
        
    
        def shiftperbin(fshift,i,fmultiple,fsigma):
            better_flx = starting_flx[whit_index[int(i)]] * fmultiple
            better_g = ss.gaussian(gauss_elements,fsigma)
            better_kernel = better_g / sum(better_g)
            better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                                   slice_iof, mode='same'))
            better_y2 = si.splev(wav[whit_index[int(i)]]+fshift,better_tck)
            return sum((better_y2 - better_flx) ** 2 / \
                        (fmultiple ) ** 2)

    
        wav_bin_av = []
        calib_fit = []
        calib_err = []
        resolution_fit = []
        resolution_err = []
        bin_call = []
    
        for i in np.arange(len(whit_index)):
            if(len(wav[whit_index[i]]) < 100):
                continue
            bin_call = mi.Minuit(shiftperbin,\
                fshift=m.values["fshift"],\
                i=i,\
                fix_i=True,\
                fmultiple=m.values["fmultiple"],\
                fix_fmultiple=True,\
                fsigma=m.values["fsigma"],\
                # fix_fsigma=True,\ #test commented 
                strategy=2)
            bin_call.migrad()
            lamb = np.average(wav[whit_index[i]])
            print lamb, c_light * bin_call.values["fshift"] / lamb, bin_call.values["fshift"], bin_call.errors["fshift"]
            midpoint = np.argmin( np.abs( [ slice_iow - lamb ] ) )
        
            FWHM_prefactor = 2.0 * np.sqrt( 2.0 * np.log( 2.0 ) )
            FWHM = FWHM_prefactor * bin_call.values["fsigma"] 
            FWHM_err1 = FWHM_prefactor * ( bin_call.values["fsigma"] + bin_call.errors["fsigma"] )
            FWHM_err2 = FWHM_prefactor * ( bin_call.values["fsigma"] - bin_call.errors["fsigma"] )
            
            del_lamb = ( slice_iow[midpoint + gauss_elements / 2 ]
                        - slice_iow[ midpoint - gauss_elements / 2] ) * FWHM / gauss_elements
            del_lamb_err1 = ( slice_iow[midpoint + gauss_elements / 2 ]
                            - slice_iow[ midpoint - gauss_elements / 2] ) * FWHM_err1 / gauss_elements
            del_lamb_err2 = ( slice_iow[midpoint + gauss_elements / 2 ]
                            - slice_iow[ midpoint - gauss_elements / 2] ) * FWHM_err2 / gauss_elements
            wav_bin_av.append(lamb)
            calib_fit.append(bin_call.values["fshift"])    
            calib_err.append(bin_call.errors["fshift"])
            # print lamb, del_lamb, bin_call.values["fsigma"]
            resolution_fit.append( lamb / del_lamb )
            resolution_err.append( np.abs( lamb/del_lamb_err1 - lamb/del_lamb_err2) /2.0 ) 
    
        wav_bin_av = np.array(wav_bin_av) 
        calib_fit = np.array(calib_fit)
        calib_err = np.array(calib_err)
        resolution_fit = np.array(resolution_fit)
        resolution_err = np.array(resolution_err)
        
        firstsplit = data_file.split("/")
        secondsplit = firstsplit[-1].split(".")
        naming = "." + secondsplit[1] + "." + secondsplit[2] + "." + secondsplit[3] 
        
        outfilename = "Test."
        outfilename += astro_object
        outfilename += naming 
        outfilename += ".ascii" 

        resolutionfile = "Test2." 
        resolutionfile += astro_object
        resolutionfile += naming 
        resolutionfile += ".ascii" 
        
        
        velocity_fit = c_light * calib_fit / wav_bin_av
        velocity_err = c_light * calib_err / wav_bin_av

        f = open(outfilename, 'w')
        for j in range(len(wav_bin_av)):
            print >>f, wav_bin_av[j], velocity_fit[j], velocity_err[j], wav[len(wav)/2], m.values["fshift"] * order_velocity, \
                    m.errors["fshift"] * order_velocity

        f.close()
        g = open(resolutionfile,'w')
        for j in range(len(wav_bin_av)):
            print >>g, wav_bin_av[j], resolution_fit[j], resolution_err[j] 
        
        g.close()
    
    
        if(len(calib_fit) > 0):
            plotfilename = "Plots/Resolution/" + astro_object + naming + ".pdf" 
        
            titlestr = "Exposure: " + secondsplit[0] + " Chip: " + secondsplit[1] + " Order: " + secondsplit[2]
        
            pl.subplot(311)
            pl.ylabel("Normalized Flux")
            pl.title(titlestr)
            pl.plot(wav,flx * m.values["fmultiple"])
            temp = pl.xlim()
            pl.subplot(312)
            pl.ylabel("v_cal (m/s)")
            pl.errorbar(wav_bin_av,velocity_fit,yerr=velocity_err)
            pl.plot(wav_bin_av,np.zeros(len(wav_bin_av)))
            pl.xlim(temp)
            pl.subplot(313)
            pl.ylabel("R value")
            # pl.ylabel("Signal/Noise")
            # pl.plot(wav,flx/err,wav,np.zeros(len(wav)))
            # pl.plot(wav_bin_av,resolution_fit)
            pl.errorbar(wav_bin_av,resolution_fit,yerr=resolution_err)
            pl.xlim(temp)
            # plotfilename += "pdf"
            pl.savefig(plotfilename)
            pl.close()

    
    logfile.close()

    # Run from uves/HR1996