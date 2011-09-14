#!/usr/bin/env python
# encoding: utf-8
"""
cleanup.py

Created by Jonathan Whitmore on 2010-05-26.
Modified: 
"""
version = 1.0

from config_calibration import *

# =======================================
# = Make List for Data to be Calibrated =
# =======================================

today = datetime.date.today()
CleanupLogFile = "Logs/cleanup." + str(today) + ".ascii"
CleanupLogFileHandle = open(CleanupLogFile, 'a')    

print "Beginning cleanup of data...", datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")

ToCleanList = glob.glob(RawCopyDir + "*")

for CleanDataFile in ToCleanList:
    CleanDataReader = csv.reader(open(CleanDataFile), delimiter=' ')
    #   Reset variables
    wav = []
    flx = []
    err = []
    pix = []
    #------------------------------------------
    #   Read in this iteration's file
    #------------------------------------------
    for row in CleanDataReader:
        wav.append(float(row[0]))
        try: flx.append(float(row[1]))
        except: 
            print "Remove ",CleanDataFile, "it has serious problems"
            break
        err.append(float(row[2]))
        if len(row) > 3:
            pix.append(float(row[3]))
        
    wav = np.array(wav)
    flx = np.array(flx)
    err = np.array(err)
    pix = np.array(pix)
    
    # all data points that have the required S/N and have an error value larger than zero.
    tcleanwav = wav[np.where(err > error_lower)]
    tcleanflx = flx[np.where(err > error_lower)]
    tcleanerr = err[np.where(err > error_lower)]
    if len(wav) == len(pix):
        tcleanpix = pix[np.where(err > error_lower)]
    
    cleanwav = tcleanwav[np.where(tcleanflx/tcleanerr > sn_lower)]
    cleanflx = tcleanflx[np.where(tcleanflx/tcleanerr > sn_lower)]
    cleanerr = tcleanerr[np.where(tcleanflx/tcleanerr > sn_lower)]
    if len(wav) == len(pix):
        cleanpix = tcleanpix[np.where(tcleanflx/tcleanerr > sn_lower)]

    # cleanwav = wav[np.where(flx[np.where(err > error_lower)]/err[np.where(err > error_lower)] > sn_lower)]
    # cleanflx = flx[np.where(flx[np.where(err > error_lower)]/err[np.where(err > error_lower)] > sn_lower)]
    # cleanerr = err[np.where(flx[np.where(err > error_lower)]/err[np.where(err > error_lower)] > sn_lower)]
    # if len(wav) == len(pix):
    #     cleanpix = pix[np.where(flx[np.where(err > error_lower)]/err[np.where(err > error_lower)] > sn_lower)]
        
    
    # Remove odd discrepancies near the edge
    cleanwav = cleanwav[edgebuffer:-edgebuffer]
    cleanflx = cleanflx[edgebuffer:-edgebuffer]
    cleanerr = cleanerr[edgebuffer:-edgebuffer]
    if len(wav) == len(pix):
        cleanpix = cleanpix[edgebuffer:-edgebuffer]
    
    # Sometimes this makes things too small, so punt early 
    if len(cleanwav) < 200:
        continue

    # killwav is cleanwav with forbidden wavelengths killed
    killwav = []
    killflx = []
    killerr = []
    killpix = []
    
    for j in range(len(cleanwav)):
        test = 0
        # Check if the wavelength is in any of the forbidden windows
        for k in range(len(begin_kill_array)):
            if cleanwav[j] > begin_kill_array[k] and cleanwav[j] < end_kill_array[k]:
                test = 1
        if test == 0: 
            killwav.append(cleanwav[j])
            killflx.append(cleanflx[j])
            killerr.append(cleanerr[j])
            if len(wav) == len(pix):
                killpix.append(cleanpix[j])
    
    killwav = np.array(killwav)
    killflx = np.array(killflx)
    killerr = np.array(killerr)
    killpix = np.array(killpix)
    print >>CleanupLogFileHandle, cleanwav[0], cleanwav[-1], len(cleanwav) - len(killwav)

    firstsplit = CleanDataFile.split("/")
    secondsplit = firstsplit[-1].split(".")
    naming = "." + secondsplit[0] + "." + secondsplit[1] + "." + secondsplit[2] 
    RawMaskedFile = RawMaskedDir + astro_object + naming + ".ascii" 
    
    RawMaskedFileHandle = open(RawMaskedFile, 'w')
    for j in range(len(killwav)):
        if len(killwav) == len(killpix):
            print >>RawMaskedFileHandle, killwav[j], killflx[j], killerr[j], killpix[j]
            continue
        print >>RawMaskedFileHandle, killwav[j], killflx[j], killerr[j] 
    
    RawMaskedFileHandle.close()
    

print "Finished cleaning up. ", datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")