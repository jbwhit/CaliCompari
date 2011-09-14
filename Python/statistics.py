#!/usr/bin/env python
# encoding: utf-8
"""
statistics.py

Created by Jonathan Whitmore on 2010-05-13.
Modified: 
"""
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


parser = OptionParser()
parser.add_option("-e","--exposure", 
                  action="store", type="string", dest="exponame", 
                  help="name exposures", metavar="EXPO")
(options, args) = parser.parse_args()


# =============
# = Constants =
# =============
c_light = 299792458.

# =============
# = Variables =
# =============
clipping_limit = 3.0

for expo in CALIBRATION_DATA:
    COUNTER = COUNTER + 1
    plotx.append([])
    ploty.append([])
    ploterr.append([])
    plotdev.append([])
    DATA_LIST = glob.glob(expo)
    # DATA_LIST = glob.glob(CALIBRATION_DATA)

    today = datetime.date.today()

    wav     = []
    cal     = []
    err     = []
    tot_wav = []
    tot_cal = []
    tot_err = []
    name = []
    deviation = []
    clipped_deviation = []
    i = 0

    for data_file in DATA_LIST:
        # 1. bin_wavelength 2. shift in m/s 3. error on shift 
        # 4. midpoint wavelength 5. full order calibration 6. error
        wav.append([])
        cal.append([])
        err.append([])
        tot_wav.append([])
        tot_cal.append([])
        tot_err.append([])
        deviation.append([])
        clipped_deviation.append([])
        wav[i], cal[i], err[i], tot_wav[i], tot_cal[i], tot_err[i] = np.loadtxt(data_file,unpack='True')
        tot_wav[i] = np.average(tot_wav[i])
        tot_cal[i] = np.average(tot_cal[i])
        tot_err[i] = np.std(cal[i])
        deviation[i] = np.max(cal[i]) - np.min(cal[i])
        clipped_err = cal[i][np.where(np.abs(err[i] - np.average(err[i])) < clipping_limit * np.std(err[i]))]
        clipped_deviation[i] = np.max(clipped_err) - np.min(clipped_err)
        first_split = data_file.split("/")
        second_split = first_split[-1].split(".")
        filename = str("Ex: " + second_split[0] + " Chip: " + second_split[1] + " Order: " + second_split[2])
        name.append(second_split[0]+"." + second_split[1]+"." + second_split[2])
        print >>f, name[i], np.average(cal[i]), np.std(cal[i]), deviation[i], clipped_deviation[i]
        i = i + 1
        
    # TODO make consistent: err and std 

    plotx[COUNTER] = tot_wav
    ploty[COUNTER] = tot_cal
    ploterr[COUNTER] = tot_err
    plotdev[COUNTER] = clipped_deviation
    # print name
    
plotfilename = "Plots/Calibration/scatter.pdf"
titlestr = "HR1996 -- VLT I2 Exposures"
pl.title(titlestr)
pl.scatter(plotx[0],plotdev[0], color='b', label="Ex1")
pl.scatter(plotx[1],plotdev[1], color='g', label="Ex2")
pl.scatter(plotx[2],plotdev[2], color='r', label="Ex3")
pl.scatter(plotx[3],plotdev[3], color='k', label="Ex7")
pl.scatter(plotx[4],plotdev[4], color='purple', label="Ex8")
pl.scatter(plotx[5],plotdev[5], color='gray', label="Ex9")
pl.xlabel("Wavelength in Angstroms")
pl.ylabel("v_cal (m/s)")
pl.legend(loc='lower right')
pl.ylim((0,500))
pl.savefig(plotfilename)
pl.close()

plotfilename = "Plots/Calibration/errorbar.pdf"
pl.title(titlestr)
pl.errorbar(plotx[0],ploty[0], yerr=ploterr[0], color='b', label="Ex1")
pl.errorbar(plotx[1],ploty[1], yerr=ploterr[1], color='g', label="Ex2")
pl.errorbar(plotx[2],ploty[2], yerr=ploterr[2], color='r', label="Ex3")
pl.errorbar(plotx[3],ploty[3], yerr=ploterr[3], color='k', label="Ex7")
pl.errorbar(plotx[4],ploty[4], yerr=ploterr[4], color='purple', label="Ex8")
pl.errorbar(plotx[5],ploty[5], yerr=ploterr[5], color='gray', label="Ex9")
pl.xlabel("Wavelength in Angstroms")
pl.ylabel("v_cal (m/s)")
pl.legend( loc='lower right')
pl.ylim((-2000,2000))
pl.savefig(plotfilename)
pl.close()

plotfilename = "Plots/Calibration/std.pdf"
pl.title("Standard Deviations")
pl.plot(plotx[0],ploterr[0], color='b', label="Ex1")
pl.plot(plotx[1],ploterr[1], color='g', label="Ex2")
pl.plot(plotx[2],ploterr[2], color='r', label="Ex3")
pl.plot(plotx[3],ploterr[3], color='k', label="Ex7")
pl.plot(plotx[4],ploterr[4], color='purple', label="Ex8")
pl.plot(plotx[5],ploterr[5], color='gray', label="Ex9")
pl.xlabel("Wavelength in Angstroms")
pl.ylabel("std of v_cal (m/s)")
pl.legend( loc='upper left')
pl.ylim((0,200))
pl.savefig(plotfilename)
pl.close()
# TODO Lower and upper chips -- separate statistics

