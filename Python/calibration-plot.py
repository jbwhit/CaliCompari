#!/usr/bin/env python
# encoding: utf-8
"""
calibration-plot.py

Created by Jonathan Whitmore on 2010-05-27.
Modified: 2010-06-01
"""
version = 1.0

from config_calibration import *

c_light = 299792458.

# Create dictionaries for wavelength, calibration, and error
avwav = {}
avcal = {}
averr = {}
avpix = {}
avres = {}
avres_err = {}

# Create dictionaries for wavelength, flux, continuum, etc. 
wav = {}
flx = {}
err = {}
con = {}
pix = {}

ExpList = exposures_to_analyze
ExpChipList = []
ExpChipOrderList = []

for core_exposure in exposures_to_analyze:
    ExpStr = core_exposure
    print core_exposure
    tempexpwav = []
    tempexpcal = []
    tempexperr = []
    tempexpres = []
    tempexpres_err = []
    tempexppix = []
    
    tempconexpwav = []
    tempconexpflx = []
    tempconexperr = []
    tempconexpcon = []
    tempconexppix = []
    for chip in chips:
        print chip
        tempexpchipwav = []
        tempexpchipcal = []
        tempexpchiperr = []
        tempexpchipres = []
        tempexpchipres_err = []
        tempexpchippix = []
        ExpChipStr = ExpStr + "." + chip
        ExpChipList.append(ExpChipStr)
        CalibrationFileList = glob.glob(OutputCalibrationDir + astro_object + "." + ExpChipStr + "*")
        for CalibrationFile in CalibrationFileList:
            ExpChipOrderStr = '.'.join(CalibrationFile.split("/")[-1].split(".")[1:4])
            ExpChipOrderList.append(ExpChipOrderStr)
            tempexpchipordwav = []
            tempexpchipordcal = []
            tempexpchiporderr = []
            tempexpchipordres = []
            tempexpchipordres_err = []
            tempexpchipordpix = []            
            
            # print CalibrationFile
            
            for row in csv.reader(open(CalibrationFile),delimiter=' '):
                tempexpchipordwav.append(float(row[0]))
                tempexpchipordcal.append(float(row[1]))
                tempexpchiporderr.append(float(row[2]))
                tempexpchipordres.append(float(row[3]))
                tempexpchipordres_err.append(float(row[4]))
                if len(row) > 5:
                    tempexpchipordpix.append(float(row[5]))
            avwav[ExpChipOrderStr] = np.array(tempexpchipordwav)
            avcal[ExpChipOrderStr] = np.array(tempexpchipordcal)
            averr[ExpChipOrderStr] = np.array(tempexpchiporderr)
            avres[ExpChipOrderStr] = np.array(tempexpchipordres)
            avres_err[ExpChipOrderStr] = np.array(tempexpchipordres_err)
            if len(tempexpchipordpix) > 0:
              avpix[ExpChipOrderStr] = np.array(tempexpchipordpix)
            
            tempexpchipwav.append(np.array(tempexpchipordwav))
            tempexpchipcal.append(np.array(tempexpchipordcal))
            tempexpchiperr.append(np.array(tempexpchiporderr))
            tempexpchipres.append(np.array(tempexpchipordres))
            tempexpchipres_err.append(np.array(tempexpchipordres_err))
            if len(tempexpchipordpix) > 0:
              tempexpchippix.append(np.array(tempexpchipordpix))
        
        # For exposure names assign the avwav 
        print len(tempexpchipwav)
        
        avwav[ExpChipStr] = np.concatenate(tempexpchipwav)
        avcal[ExpChipStr] = np.concatenate(tempexpchipcal)
        averr[ExpChipStr] = np.concatenate(tempexpchiperr)
        avres[ExpChipStr] = np.concatenate(tempexpchipres)
        avres_err[ExpChipStr] = np.concatenate(tempexpchipres_err)
        avpix[ExpChipStr] = np.concatenate(tempexpchippix)
        
        tempexpwav.append(np.concatenate(tempexpchipwav))
        tempexpcal.append(np.concatenate(tempexpchipcal))
        tempexperr.append(np.concatenate(tempexpchiperr))
        tempexpres.append(np.concatenate(tempexpchipres))
        tempexpres_err.append(np.concatenate(tempexpchipres_err))
        tempexppix.append(np.concatenate(tempexpchippix))
        
        tempconexpchipwav = []
        tempconexpchipflx = []
        tempconexpchiperr = []
        tempconexpchipcon = []
        tempconexpchippix = []
        ContinuumFileList = glob.glob(OutputContinuumDir + astro_object + "." + ExpChipStr + "*")

        for ContinuumFile in ContinuumFileList:
            ExpChipOrderStr = '.'.join(ContinuumFile.split("/")[-1].split(".")[1:4])
            tempconexpchipordwav = []
            tempconexpchipordflx = []
            tempconexpchiporderr = []
            tempconexpchipordcon = []
            tempconexpchipordpix = []
            
            for row in csv.reader(open(ContinuumFile),delimiter=' '):
                tempconexpchipordwav.append(float(row[0]))
                tempconexpchipordflx.append(float(row[1]))
                tempconexpchiporderr.append(float(row[2]))
                tempconexpchipordcon.append(float(row[3]))
                if len(row) > 4:
                    tempconexpchipordpix.append(float(row[4]))
            wav[ExpChipOrderStr] = np.array(tempconexpchipordwav)
            flx[ExpChipOrderStr] = np.array(tempconexpchipordflx)
            err[ExpChipOrderStr] = np.array(tempconexpchiporderr)
            con[ExpChipOrderStr] = np.array(tempconexpchipordcon)
            pix[ExpChipOrderStr] = np.array(tempconexpchipordpix)
        
            tempconexpchipwav.append(tempconexpchipordwav)
            tempconexpchipflx.append(tempconexpchipordflx)
            tempconexpchiperr.append(tempconexpchiporderr)
            tempconexpchipcon.append(tempconexpchipordcon)
            tempconexpchippix.append(tempconexpchipordpix)

        wav[ExpChipStr] = np.concatenate(tempconexpchipwav)
        flx[ExpChipStr] = np.concatenate(tempconexpchipflx)
        err[ExpChipStr] = np.concatenate(tempconexpchiperr)
        con[ExpChipStr] = np.concatenate(tempconexpchipcon)
        pix[ExpChipStr] = np.concatenate(tempconexpchippix)

        tempconexpwav.append(np.concatenate(tempconexpchipwav))
        tempconexpflx.append(np.concatenate(tempconexpchipflx))
        tempconexperr.append(np.concatenate(tempconexpchiperr))
        tempconexpcon.append(np.concatenate(tempconexpchipcon))
        tempconexppix.append(np.concatenate(tempconexpchippix))
        
        # TODO         Print exp.chip avcal, std, (weightedstd), 
            # sigma-clipped std, S/N, Res, np.median(avres)
        
    
    avwav[ExpStr] = np.concatenate(tempexpwav)
    avcal[ExpStr] = np.concatenate(tempexpcal)
    averr[ExpStr] = np.concatenate(tempexperr)
    avres[ExpStr] = np.concatenate(tempexpres)
    avres_err[ExpStr] = np.concatenate(tempexpres_err)
    avpix[ExpStr] = np.concatenate(tempexppix)
    
    wav[ExpStr] = np.concatenate(tempconexpwav)
    flx[ExpStr] = np.concatenate(tempconexpflx)
    err[ExpStr] = np.concatenate(tempconexperr)
    con[ExpStr] = np.concatenate(tempconexpcon)
    pix[ExpStr] = np.concatenate(tempconexppix)
    
    
def plot_histogram(key1,key2,key3):
    """docstring for histogram"""
    pl.hist(avcal[key1],bins=100,range=[-500,1000],histtype="step",color="blue",alpha=0.5)
    pl.hist(avcal[key2],bins=100,range=[-500,1000],histtype="step",color="green",alpha=0.2)
    pl.hist(avcal[key3],bins=100,range=[-500,1000],histtype="step",color="red")
    pl.hist(avcal[key1],bins=100,range=[-500,1000],histtype="step",color="blue",alpha=0.5,label="Expo 1",fill=True)
    pl.hist(avcal[key2],bins=100,range=[-500,1000],histtype="step",color="green",alpha=0.2,label="Expo 2",fill=True)
    pl.hist(avcal[key3],bins=100,range=[-500,1000],histtype="step",color="red",label="Expo 3",fill=True)
    pl.xlabel("Calibration Shift in m/s")
    pl.ylabel("Counts")
    pl.title("Histogram of vshift by Exposure")
    pl.legend()
    # pl.show()
    PdfFilename = PdfDir + "Analysis/" + "histogram.pdf"
    PsFilename = PsDir + "Analysis/" + "histogram.eps"
    pl.savefig(PdfFilename)
    pl.savefig(PsFilename)
    print "If you used transparency, use GIMP to turn pdf into eps w/ transparency."
    pass
    
AnalysisFile = "Summary/Analysis." + astro_object + ".ascii"
AnalysisFileHandle = open(AnalysisFile, 'w')
print >>AnalysisFileHandle, "key, np.average(avcal[key]), np.std(avcal[key]), \
    np.average(averr[key]), np.median(avres[key]), np.average(avres_err[key]), \
    np.median(flx[key]/err[key]), weighted_av_and_std(avcal[key],averr[key])"
for key in ExpChipList:
    print >>AnalysisFileHandle, key, np.average(avcal[key]), np.std(avcal[key]), \
        np.average(averr[key]), np.median(avres[key]), np.average(avres_err[key]), \
        np.median(flx[key]/err[key]), weighted_av_and_std(avcal[key],averr[key])

AnalysisFileHandle.close()


# pl.hist(avcal['01'],bins=100,range=[-500,1000],histtype="step",color="blue",alpha=0.2,label="Expo 3")
# pl.hist(avcal['01'],bins=100,range=[-500,1000],histtype="step",color="blue",alpha=0.2,label="Expo 3",fill=True)

# pl.savefig('Histogram1.eps')
# pl.close()

# -----
#     
# pdfplotfilename = PdfDir + astro_object + "." + "all" + ".pdf" 
# psplotfilename = PsDir + astro_object + "." + "all" + ".eps" 
# titlestr = astro_object + " All Exposures"
# pl.ylabel("v_cal (m/s)")
# pl.xlabel("Wavelength in Angstroms")
# pl.title(titlestr)
# for i in range(len(exposures_to_analyze)):
#     colorname = "blue"
#     label = None
#     if i == 0:
#         colorname = "blue"
#         label = exposures_to_analyze[i]
#     if i == 1:
#         colorname = "green"
#         label = exposures_to_analyze[i]
#     if i == 2:
#         colorname = "red"
#         label = exposures_to_analyze[i]
#     if i == 3:
#         colorname = "black"
#         label = exposures_to_analyze[i]
#     if i == 4:
#         colorname = "purple"
#         label = exposures_to_analyze[i]
#     if i == 5:
#         colorname = "magenta"
#         label = exposures_to_analyze[i]
#     for j in range(len(expo_wav[i])):
#         pl.errorbar(expo_wav[i][j],expo_cal[i][j],yerr=expo_err[i][j],color=colorname)
#         pl.errorbar(expo_wav[i][j],expo_cal[i][j],yerr=expo_err[i][j],color=colorname,label=label)
# 
# if pl.ylim()[-1] > 3000.0:
#     pl.ylim(pl.ylim()[0],3000.0)
# pl.legend()
# pl.savefig(pdfplotfilename)
# pl.savefig(psplotfilename)
# pl.close()
# 
# slope_list = []
# slopefile = open("Summary/slope.ascii", 'w')
# for i in range(len(exposures_to_analyze)):
#     slope_list.append(linregress(plot_wav[i], plot_cal[i])[0])
#     print >>slopefile,exposures_to_analyze[i], linregress(plot_wav[i], plot_cal[i])[0], "m/s per Angstrom"
# print >>slopefile, "Average: ", np.average(slope_list), "m/s per Angstrom"
# slopefile.close()
# 
# if len(pix) == len(wav):
#     pixpdfplotfilename = PdfDir + astro_object + "." + "pix" + ".pdf" 
#     pixpsplotfilename = PsDir + astro_object + "." + "pix" + ".eps" 
#     titlestr = astro_object + " Pixels Exposures"
#     pl.ylabel("v_cal (m/s)")
#     pl.xlabel("Pixel Number")
#     pl.title(titlestr)
#     for i in range(len(exposures_to_analyze)):
#         colorname = "blue"
#         if i == 0:
#             colorname = "blue"
#         if i == 1:
#             colorname = "green"
#         if i == 2:
#             colorname = "red"
#         for j in range(len(expo_wav[i])):
#             pl.scatter(expo_pix[i][j],expo_cal[i][j],color=colorname)
#             # pl.errorbar(expo_pix[i][j],expo_cal[i][j],yerr=expo_err[i][j],color=colorname)
#     pl.savefig(pixpdfplotfilename)
#     pl.savefig(pixpsplotfilename)
#     pl.close()
# 
#     # Test if the OLD directory exists -- if so, run: 
# if len(glob.glob("OLD")) > 0:
#     
#     old_wav = []
#     old_flx = []
#     old_err = []
#     # This is a one-off w/ 3 exposures
#     old_wav.append([])
#     old_flx.append([])
#     old_err.append([])
#     old_wav.append([])
#     old_flx.append([])
#     old_err.append([])
#     old_wav.append([])
#     old_flx.append([])
#     old_err.append([])
# 
#     old_wav[0], old_flx[0], old_err[0], old1,old2,old3 = np.loadtxt('OLD/old.01.ascii', unpack=True)
#     old_wav[1], old_flx[1], old_err[1], old1,old2,old3 = np.loadtxt('OLD/old.02.ascii', unpack=True)
#     old_wav[2], old_flx[2], old_err[2], old1,old2,old3 = np.loadtxt('OLD/old.03.ascii', unpack=True)
# 
#     resultsfile = open("Logs/old.slope.ascii", 'w')
#     slope_list = []
#     for i in range(len(old_wav)):
#         slope_list.append(linregress(old_wav[i], old_flx[i])[0])
#         print >>resultsfile, linregress(old_wav[i], old_flx[i])[0], "m/s"
#     print >>resultsfile, "Average: ", np.average(slope_list), "m/s"
#     resultsfile.close()
