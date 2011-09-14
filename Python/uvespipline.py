#!/usr/bin/env python
# encoding: utf-8
"""
uvespipeline.py

Created by Jonathan Whitmore on 2010-05-03.
Modified: 2010-05-31
Comparing the wavelength calibrations with actual shifts in 
the ascii files Michael Murphy handed us.
"""
version = 1.0

from config_calibration import *

# =============
# = Constants =
# =============
c_light = 299792458.

def initial_shift(fmultiple,fshift,fsigma):
    better_flx = starting_flx * fmultiple
    better_kernel = normal_gaussian(gauss_elements,fsigma)
    better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                           slice_iof, mode='same'))
    better_y2 = si.splev(wav+fshift,better_tck)
    return sum((better_y2 - better_flx) ** 2 / \
               (fmultiple * err / con) ** 2 )

#
def fullordershift(fmultiple,fshift,fsigma):
    better_flx = starting_flx * fmultiple
    better_kernel = normal_gaussian(gauss_elements,fsigma)
    better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                           slice_iof, mode='same'))
    better_y2 = si.splev(wav+fshift,better_tck)
    return (better_y2 - better_flx) ** 2 / \
               (fmultiple * err / con) ** 2 

    # spline_tck = si.slrep(slice_iow, \
    #                         np.convolve(normal_gaussian(gauss_elements,fsigma),\
    #                         slice_iof, mode='same'))
    # test_flux = si.splev(wav + fshift, spline_tck)
    # return (test_flux - flx * fmultiple) ** 2 /\
    #             (fmultiple * err / con) ** 2
    
    
def continuum_shift(fmultiple,fshift,fsigma,cmultiple):
    better_flx = flx / (con + con * cmultiple) * fmultiple
    better_kernel = normal_gaussian(gauss_elements,fsigma)
    better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                           slice_iof, mode='same'))
    better_y2 = si.splev(wav+fshift,better_tck)
    return sum((better_y2 - better_flx) ** 2 / \
               (fmultiple * err / (con + con * cmultiple)) ** 2 )

#
def continuumshiftarray(fmultiple,fshift,fsigma,cmultiple):
    better_flx = flx / (con + con * cmultiple) * fmultiple
    better_kernel = normal_gaussian(gauss_elements,fsigma)
    better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                           slice_iof, mode='same'))
    better_y2 = si.splev(wav+fshift,better_tck)
    return (better_y2 - better_flx) ** 2 / \
               (fmultiple * err / (con + con * cmultiple)) ** 2 


def second_shift(fmultiple,fshift,fsigma,fskew):
    better_flx = starting_flx * fmultiple
    better_kernel = normal_skew_gaussian(gauss_elements,fsigma,fskew)
    better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                           slice_iof, mode='same'))
    better_y2 = si.splev(wav+fshift,better_tck)
    return sum((better_y2 - better_flx) ** 2 / \
               (fmultiple * err / con) ** 2 )

#
def skewshiftarray(fmultiple,fshift,fsigma,fskew):
    better_flx = starting_flx * fmultiple
    better_kernel = normal_skew_gaussian(gauss_elements,fsigma,fskew)
    better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                           slice_iof, mode='same'))
    better_y2 = si.splev(wav+fshift,better_tck)
    return (better_y2 - better_flx) ** 2 / \
               (fmultiple * err / con) ** 2 


# bin-shifting fmultiple,fshift,fsigma
def shiftperbin(i,fmultiple,fshift,fsigma):
    better_flx = starting_flx[whit_index[int(i)]] * fmultiple
    better_kernel = normal_gaussian(gauss_elements,fsigma)
    better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                           slice_iof, mode='same'))
    better_y2 = si.splev(wav[whit_index[int(i)]]+fshift, better_tck)
    return sum((better_y2 - better_flx) ** 2 / \
                (fmultiple * starting_err[whit_index[int(i)]] ) ** 2)
                # (fmultiple * err[whit_index[int(i)]] / con[whit_index[int(i)]]) ** 2)

#
def shiftperbinarray(i,fmultiple,fshift,fsigma):
    better_flx = starting_flx[whit_index[int(i)]] * fmultiple
    better_kernel = normal_gaussian(gauss_elements,fsigma)
    better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                           slice_iof, mode='same'))
    better_y2 = si.splev(wav[whit_index[int(i)]]+fshift, better_tck)
    return (better_y2 - better_flx) ** 2 / \
                (fmultiple * err[whit_index[int(i)]] / con[whit_index[int(i)]]) ** 2

#   
def plotresidualshiftperbin(fshift,i,fmultiple,fsigma):
    better_flx = starting_flx[whit_index[int(i)]] * fmultiple
    better_kernel = normal_gaussian(gauss_elements,fsigma)
    better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                           slice_iof, mode='same'))
    better_y2 = si.splev(wav[whit_index[int(i)]]+fshift, better_tck)
    pl.plot((better_y2 - better_flx) )
    return

#   
def plotsidebyside(i):
    pl.plot(wav[fitbinindex],better_y2,wav[fitbinindex],flx[fitbinindex] * fitbinmultiple)
    
def plotcontinuum_shift(fmultiple,fshift,fsigma,cmultiple):
    better_flx = flx / (con + con * cmultiple) * fmultiple
    better_kernel = normal_gaussian(gauss_elements,fsigma)
    better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                           slice_iof, mode='same'))
    better_y2 = si.splev(wav+fshift,better_tck)
    return sum((better_y2 - better_flx) ** 2 / \
               (fmultiple * err / (con + con * cmultiple)) ** 2 )

def skewshiftperbin(fshift,i,fmultiple,fsigma,fskew):
    better_flx = starting_flx[whit_index[int(i)]] * fmultiple
    better_kernel = normal_skew_gaussian(gauss_elements,fsigma,fskew)
    better_tck = si.splrep(slice_iow, np.convolve(better_kernel,\
                           slice_iof, mode='same'))
    better_y2 = si.splev(wav[whit_index[int(i)]]+fshift, better_tck)
    return sum((better_y2 - better_flx) ** 2 / \
                (fmultiple * starting_err[whit_index[int(i)]]) ** 2)

#

CalibrationSummaryDir = "Summary/calibration."
AsciiSpectrographDir = "SpectrographAscii/"
# ====================
# = Read in FTS Data =
# ====================

today = datetime.date.today()

print "Beginning calibration analysis of data..."

output1 = open("SpectrographAscii/output1.ascii",'w')
output2 = open("SpectrographAscii/output2.ascii",'w')

for core_exposure in exposures_to_analyze:
    CalibrationDataList = glob.glob("Differences/" + core_exposure + "*")
    print datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")

    for CalibrationDataFile in CalibrationDataList:
        #   Iterate over the Data files
        
        pix = []
        diff1 = []
        diff2 = []
        diff3 = []
        wavelength = []
        
        CalibrationDataReader = csv.reader(open(CalibrationDataFile), delimiter=' ')
        for row in CalibrationDataReader:
            pix.append(float(row[0]))
            diff1.append(float(row[1]))
            diff2.append(float(row[2]))
            diff3.append(float(row[3]))
            wavelength.append(float(row[4]))
        
        pix = np.array(pix)
        diff1 = np.array(diff1)
        diff2 = np.array(diff2)
        diff3 = np.array(diff3)
        wavelength = np.array(wavelength)
        
        firstsplit = CalibrationDataFile.split("/")
        secondsplit = firstsplit[-1].split(".")
        naming = "." + secondsplit[0] + "." + secondsplit[1] + "." + secondsplit[2] 
        nodotnaming = secondsplit[0] + "." + secondsplit[1] + "." + secondsplit[2] 
        

        deg4 = "HE0515"
        deg5 = "HE0515-deg5"
        deg6 = "HE0515-deg6"
        key4 = "4" + naming
        key5 = "5" + naming
        key6 = "6" + naming
        key45 = "45" + naming
        key46 = "46" + naming
        key56 = "56" + naming
        
        tempwav = []
        tempcal = []
        temperr = []
        CalFile = "CalibrationWhitmore/" + deg4 + naming + ".ascii"

        for row in csv.reader(open(CalFile), delimiter=' '):
            tempwav.append(float(row[0]))
            tempcal.append(float(row[1]))
            temperr.append(float(row[2]))
        
        wav = {}
        cal = {}
        err = {}
        wav[key4] = np.array(tempwav)
        cal[key4] = np.array(tempcal)
        err[key4] = np.array(temperr)
        tempwav = []
        tempcal = []
        temperr = []
        CalFile = "CalibrationWhitmore/" + deg5 + naming + ".ascii"
        for row in csv.reader(open(CalFile), delimiter=' '):
            tempwav.append(float(row[0]))
            tempcal.append(float(row[1]))
            temperr.append(float(row[2]))
        
        wav[key5] = np.array(tempwav)
        cal[key5] = np.array(tempcal)
        err[key5] = np.array(temperr)
        tempwav = []
        tempcal = []
        temperr = []
        CalFile = "CalibrationWhitmore/" + deg6 + naming + ".ascii"
        for row in csv.reader(open(CalFile), delimiter=' '):
            tempwav.append(float(row[0]))
            tempcal.append(float(row[1]))
            temperr.append(float(row[2]))
        
        wav[key6] = np.array(tempwav)
        cal[key6] = np.array(tempcal)
        err[key6] = np.array(temperr)        
        
        
        PdfCalibrationFile = "PDF/" + astro_object + naming + ".pdf" 
        PsCalibrationFile = "PS/" + astro_object + naming + ".eps" 
        titlestr = "Exposure: " + secondsplit[0] + " Chip: " + secondsplit[1] + " Order: " + secondsplit[2]
        pl.ylabel("Difference in m/s")
        pl.title(titlestr)
        pl.scatter(wav[key4],cal[key5]-cal[key4],color="blue",label="4-5")
        pl.scatter(wav[key4],cal[key6]-cal[key4],color="green",label="4-6")
        pl.scatter(wav[key4],cal[key6]-cal[key5],color="red",label="5-6")
        tempy = pl.ylim()
        pl.plot(wavelength, diff1,label="deg4 - deg5")
        pl.plot(wavelength, diff2,label="deg4 - deg6")
        pl.plot(wavelength, diff3,label="deg6 - deg5")
        pl.ylim(tempy)
        pl.savefig(PdfCalibrationFile)
        pl.savefig(PsCalibrationFile)
        pl.ylim((-100,100))
        pl.savefig("PDF/100." + astro_object + naming + ".pdf" )
        pl.savefig("PS/100." + astro_object + naming + ".eps" )
        pl.close()
        
        vdiff = {}
        tempvdiff45 = []
        tempvdiff46 = []
        tempvdiff56 = []

        interp45 = []
        interp46 = []
        interp56 = []
        # I use three points because when I find the min point, I don't know if I'm on the left or right
        # clearly, possible to fix if it mattered.
        for thingy in range(len(wav[key4])):
            temp1 = np.argmin(np.abs(wavelength - wav[key4][thingy]))
            temp45 = linregress(wavelength[temp1-1:temp1+2], diff1[temp1-1:temp1+2])
            temp46 = linregress(wavelength[temp1-1:temp1+2], diff2[temp1-1:temp1+2])
            temp56 = linregress(wavelength[temp1-1:temp1+2], diff3[temp1-1:temp1+2])
            interp45.append(temp45[0]*wav[key4][thingy] + temp45[1])
            interp46.append(temp46[0]*wav[key4][thingy] + temp46[1])
            interp56.append(temp56[0]*wav[key4][thingy] + temp56[1])
            tempvdiff45.append( cal[key4][thingy] - cal[key5][thingy] + (temp45[0]*wav[key4][thingy] + temp45[1]) )
            tempvdiff46.append( cal[key4][thingy] - cal[key6][thingy] + (temp46[0]*wav[key4][thingy] + temp46[1]) )
            tempvdiff56.append( cal[key5][thingy] - cal[key6][thingy] + (temp56[0]*wav[key4][thingy] + temp56[1]) )
        
        vdiff[key45] = np.array(tempvdiff45)
        vdiff[key46] = np.array(tempvdiff46)
        vdiff[key56] = np.array(tempvdiff56)
        
        Output45File = AsciiSpectrographDir + key45 + ".ascii"
        Output46File = AsciiSpectrographDir + key46 + ".ascii"
        Output56File = AsciiSpectrographDir + key56 + ".ascii"

        Output45FileHandle = open(Output45File, 'w')
        Output46FileHandle = open(Output46File, 'w')
        Output56FileHandle = open(Output56File, 'w')
        
        for thingy in range(len(wav[key4])):
            print >>Output45FileHandle, wav[key4][thingy], \
                cal[key4][thingy], err[key4][thingy], cal[key5][thingy], err[key5][thingy], \
                interp45[thingy], vdiff[key45][thingy]
            print >>Output46FileHandle, wav[key4][thingy], \
                cal[key4][thingy], err[key4][thingy], cal[key6][thingy], err[key6][thingy], \
                interp46[thingy], vdiff[key46][thingy]
            print >>Output56FileHandle, wav[key4][thingy], \
                cal[key5][thingy], err[key5][thingy], cal[key6][thingy], err[key6][thingy], \
                interp56[thingy], vdiff[key56][thingy]            

        Output45FileHandle.close()
        Output46FileHandle.close()
        Output56FileHandle.close()
        
        print >>output1, key45, np.std(vdiff[key45]), np.std(cal[key4]), np.std(cal[key5]), np.std(cal[key4] - cal[key5]), \
            np.sqrt(.5 * (np.std(cal[key4])**2 + np.std(cal[key5])**2 - np.std(cal[key4] - cal[key5])**2 ))
            
        print >>output1, key46, np.std(vdiff[key46]), np.std(cal[key4]), np.std(cal[key6]), np.std(cal[key4] - cal[key6]), \
            np.sqrt(.5 * (np.std(cal[key4])**2 + np.std(cal[key6])**2 - np.std(cal[key4] - cal[key6])**2 ))
        
        print >>output1, key56, np.std(vdiff[key56]), np.std(cal[key5]), np.std(cal[key6]), np.std(cal[key5] - cal[key6]), \
            np.sqrt(.5 * (np.std(cal[key5])**2 + np.std(cal[key6])**2 - np.std(cal[key5] - cal[key6])**2 ))
       
        print >>output2, nodotnaming, np.sqrt(.5 * (np.std(cal[key4])**2 + np.std(cal[key5])**2 - np.std(cal[key4] - cal[key5])**2 )), \
            np.sqrt(.5 * (np.std(cal[key4])**2 + np.std(cal[key6])**2 - np.std(cal[key4] - cal[key6])**2 )), \
            np.sqrt(.5 * (np.std(cal[key5])**2 + np.std(cal[key6])**2 - np.std(cal[key5] - cal[key6])**2 ))
        
     #    print naming
     #    if naming == ".01.l.09":
     #        break
     #        
     # break
        
output1.close()
output2.close()

print "Finished. ", datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
