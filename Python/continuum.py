#!/usr/bin/env python
# encoding: utf-8
"""
continuum.py

Created by Jonathan Whitmore on 2010-05-03.
Modified: 2010-05-12
"""
version = 1.0

from config_calibration import *


# =============
# = Constants =
# =============

c_light = 299792458.

# ========================
# = Set Data Directories =
# ========================

ToContinuumList = glob.glob(RawMaskedDir + "*")

print "Beginning continuum fitting of data...", datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")

for ContinuumDataFile in ToContinuumList:
    wav = []
    flx = []
    err = []
    pix = []

    ContinuumDataReader = csv.reader(open(ContinuumDataFile), delimiter=' ')
    for row in ContinuumDataReader:
        wav.append(float(row[0]))
        flx.append(float(row[1]))
        err.append(float(row[2]))
        if len(row) > 3:
            pix.append(float(row[3]))

    wav = np.array(wav)
    flx = np.array(flx)
    err = np.array(err)
    pix = np.array(pix)
    if len(wav) < 100:
        continue
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
    
    firstsplit = ContinuumDataFile.split("/")
    secondsplit = firstsplit[-1].split(".")
    naming = "." + secondsplit[1] + "." + secondsplit[2] + "." + secondsplit[3] 
    OutputContinuumFile = OutputContinuumDir + astro_object + naming + ".ascii"
    
    # Making a fundamental change here. Printing wav, flx, err, contin fit, pixel 
    OutputContinuumFileHandle = open(OutputContinuumFile, 'w')
    for j in range(len(x0)):
        if len(wav) == len(pix):
            # print >>OutputContinuumFileHandle, x0[j], continuum_fit_y[j], continuum_fit_err[j], pix[j]
            print >>OutputContinuumFileHandle, x0[j], flx[j], err[j], continuum_points[j], pix[j]
            continue
        # print >>OutputContinuumFileHandle, x0[j], continuum_fit_y[j], continuum_fit_err[j]
        print >>OutputContinuumFileHandle, x0[j], flx[j], err[j], continuum_points[j]

    OutputContinuumFileHandle.close()
    # ==================================
    # = Plot PDF and PS of Continuum Fit =
    # ==================================
    PdfContinuumFile = PdfContinuumDir + astro_object + naming + ".pdf" 
    PsContinuumFile = PsContinuumDir + astro_object + naming + ".eps" 
    titlestr = astro_object + "Ex: " + secondsplit[1] + " Ch: " + secondsplit[2] + " Or: " + secondsplit[3]
    pl.ylabel("Flux")
    pl.title(titlestr)
    pl.plot(wav, flx, x0,interpolate.splev(x0,(part1,p2,splineorder)), 'r-')
    pl.savefig(PdfContinuumFile)
    pl.savefig(PsContinuumFile)
    pl.close()

print "Finished continuum fitting.", datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")