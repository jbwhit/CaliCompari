#!/usr/bin/env python
# encoding: utf-8
"""
template-fitting.py

Created by Jonathan Whitmore on 2010-05-03.
Modified: 2010-11-19
"""
version = 1.0

from config_calibration import *

# v0065, success = leastsq(error_function_array, v0, args=(List0065, List0065), maxfev=10000)

c_light = 299792458.

# Create dictionaries for wavelength, calibration, and error
avwav = {}
avcal = {}
averr = {}
avpix = {}
avres = {}
avres_err = {}

header = {}

ExpList = i2exposures
ExpChipList = []
ExpChipOrderList = []

from fileparsing import *

# for Exp in ExpList:
#   # print Exp
#   tempexpwav = []
#   tempexpcal = []
#   tempexperr = []
#   tempexpres = []
#   tempexpres_err = []
#   tempexppix = []
#   
#   header[Exp] = {}
#   headerfile = glob.glob("Headers/*" + Exp + ".header")
#   try:
#     headerReader = csv.reader(open(headerfile[0]), delimiter=' ')
#     for row in headerReader:
#       header[Exp][row[0]] = row[1]
#   except:
#     print "Problem with header file for ", Exp
#     
#   for chip in chips:
#     # print chip
#     tempexpchipwav = []
#     tempexpchipcal = []
#     tempexpchiperr = []
#     tempexpchipres = []
#     tempexpchipres_err = []
#     tempexpchippix = []
#     ExpChip = Exp + "." + chip
#     ExpChipList.append(ExpChip)
#     CalibrationFileList = glob.glob(OutputCalibrationDir + "*" + ExpChip + "*")
#     for CalibrationFile in CalibrationFileList:
#       ExpChipOrder = '.'.join(CalibrationFile.split("/")[-1].split(".")[2:5])
#       # print ExpChipOrder, CalibrationFile
#       ExpChipOrderList.append(ExpChipOrder)
#       tempexpchipordwav = []
#       tempexpchipordcal = []
#       tempexpchiporderr = []
#       tempexpchipordres = []
#       tempexpchipordres_err = []
#       tempexpchipordpix = []            
#       for row in csv.reader(open(CalibrationFile),delimiter=' '):
#         tempexpchipordwav.append(float(row[0]))
#         tempexpchipordcal.append(float(row[1]))
#         tempexpchiporderr.append(float(row[2]))
#         tempexpchipordres.append(float(row[3]))
#         tempexpchipordres_err.append(float(row[4]))
#         if len(row) > 5:
#           tempexpchipordpix.append(float(row[5]))
#       avwav[ExpChipOrder] = np.array(tempexpchipordwav)
#       avcal[ExpChipOrder] = np.array(tempexpchipordcal)
#       averr[ExpChipOrder] = np.array(tempexpchiporderr)
#       avres[ExpChipOrder] = np.array(tempexpchipordres)
#       avres_err[ExpChipOrder] = np.array(tempexpchipordres_err)
#       if len(tempexpchipordpix) > 0:
#         avpix[ExpChipOrderStr] = np.array(tempexpchipordpix)
#       tempexpchipwav.append(np.array(tempexpchipordwav))
#       tempexpchipcal.append(np.array(tempexpchipordcal))
#       tempexpchiperr.append(np.array(tempexpchiporderr))
#       tempexpchipres.append(np.array(tempexpchipordres))
#       tempexpchipres_err.append(np.array(tempexpchipordres_err))
#       if len(tempexpchipordpix) > 0:
#         tempexpchippix.append(np.array(tempexpchipordpix))
#     avwav[ExpChip] = np.concatenate(tempexpchipwav)
#     avcal[ExpChip] = np.concatenate(tempexpchipcal)
#     averr[ExpChip] = np.concatenate(tempexpchiperr)
#     avres[ExpChip] = np.concatenate(tempexpchipres)
#     avres_err[ExpChip] = np.concatenate(tempexpchipres_err)
#     if len(tempexpchippix) > 1:
#       avpix[ExpChip] = np.concatenate(tempexpchippix)
#     
#     tempexpwav.append(np.concatenate(tempexpchipwav))
#     tempexpcal.append(np.concatenate(tempexpchipcal))
#     tempexperr.append(np.concatenate(tempexpchiperr))
#     tempexpres.append(np.concatenate(tempexpchipres))
#     tempexpres_err.append(np.concatenate(tempexpchipres_err))
#     if len(tempexpchippix) > 1:
#       tempexppix.append(np.concatenate(tempexpchippix))
# 
#   avwav[Exp] = np.concatenate(tempexpwav)
#   avcal[Exp] = np.concatenate(tempexpcal)
#   averr[Exp] = np.concatenate(tempexperr)
#   avres[Exp] = np.concatenate(tempexpres)
#   avres_err[Exp] = np.concatenate(tempexpres_err)
#   if len(tempexppix) > 1:
#     avpix[Exp] = np.concatenate(tempexppix)

# This creates lists of all orders in exposures
List = {}
for Exp in ExpList:
  List[Exp] = []
  for i in glob.glob("Calibration/*" + "." + Exp + "*"):
    List[Exp].append(".".join(i.split("/")[1].split(".")[1:4]))

Stars = ['0064','0065','1093','1094','2092','2093']

floatlist = []
def headerfloat():
  """docstring for headerfloat"""
  for index,exposure in enumerate(i2exposures):
    for key in header[exposure]:
      try:
        header[exposure][key] = float(header[exposure][key])
        if index == 0:
          floatlist.append(key)
      except:
        pass
  pass

ParameterList = []
def uniqueheaderfloat():
  """docstring for headerfloat"""
  temp = []
  for key in header[i2exposures[0]]:
    try:
      header[i2exposures[0]][key] = float(header[i2exposures[0]][key])
      print "worked", key
    except:
      print "borked"
      continue
    if all([header[exposure][key] for exposure in i2exposures]):
      ParameterList.append(key)
      print "Added to Parameter List"
  pass

AngleList = ['XDANGL',
'AZ',
'EL',
'PARANG',
'IROT2ANG',
'ECHANGL']

def headerscatter():
  for key in ParameterList:
    pl.close()
    temp = []
    for exposure in i2exposures:
      temp.append(header[exposure][key])
    if all(x == temp[0] for x in temp):
      continue
    for exposure in i2exposures:
      pl.scatter(header[exposure][key],clipfit[exposure][2])
      pl.title(key)
    pl.savefig("Plots/Correlation/" + key + ".pdf")
  pl.close()
  pass

def starheaderscatter():
  for key in ParameterList:
    pl.close()
    temp = []
    for exposure in Stars:
      temp.append(header[exposure][key])
    if all(x == temp[0] for x in temp):
      continue
    for exposure in Stars:
      pl.scatter(header[exposure][key],clipfit[exposure][2])
      pl.title(key)
    pl.savefig("Plots/StarCorrelation/" + key + ".pdf")
  pl.close()
  pass

def Anglescatter():
  for key in AngleList:
    for exposure in i2exposures:
      pl.scatter(np.sin(np.radians(header[exposure][key])), clipfit[exposure][2])
    pl.savefig("Plots/Angle/" + key + ".pdf")
    pl.close()

def StarAnglescatter():
  for key in AngleList:
    for exposure in Stars:
      pl.scatter(np.sin(np.radians(header[exposure][key])), clipfit[exposure][2])
    pl.savefig("Plots/StarAngle/" + key + ".pdf")
    pl.close()

def shiftguess(exposure):
  """input a key, and this returns the linear fit for the shiftguess
  for all the orders in that exposure."""
  temp = []
  for key in List[exposure]:
    temp.append([np.average(avwav[key]),np.average(avcal[key])])
  initialshift[exposure] = linregress(temp)[0],linregress(temp)[1]
  pass

def maskedshiftguess(exposure):
  """input a key, and this returns the linear fit for the shiftguess
  for all the orders in that exposure."""
  temp = []
  for key in List[exposure]:
    temp.append([np.average(avwav[key]),ma.average(maskedcal[key])])
  maskedinitialshift[exposure] = linregress(temp)[0],linregress(temp)[1]
  pass
  
def multguess(exposure):
  """input a key, and this returns the linear fit for the shiftguess
  for all the orders in that exposure."""
  temp = []
  for key in List[exposure]:
    temp.append([np.average(avwav[key]),np.max(avcal[key]) - np.min(avcal[key])])
  initialmult[exposure] = linregress(temp)[0],linregress(temp)[1]
  pass

def maskedmultguess(exposure):
  """input a key, and this returns the linear fit for the shiftguess
  for all the orders in that exposure."""
  temp = []
  for key in List[exposure]:
    temp.append([np.average(avwav[key]),np.max(avcal[key]) - np.min(avcal[key])])
  maskedinitialmult[exposure] = linregress(temp)[0],linregress(temp)[1]
  pass

def clipmultguess(exposure):
  temp = []
  for key in List[exposure]:
    temp.append([np.average(avwav[key]),np.max(avcal[key]) - np.min(avcal[key])])
  initialmult[exposure] = linregress(temp)[0],linregress(temp)[1]
  pass

def templateguess(exposure):
  """Initial guess for template of an exposure"""
  shiftguess(exposure)
  multguess(exposure)
  templist = []
  for key in List[exposure]:
    templist.append((avcal[key] - ((initialshift[exposure][0]) * np.average(avwav[key]) + initialshift[exposure][1])) / (initialmult[exposure][0] * np.average(avwav[key]) + initialmult[exposure][1]))
  initialtemplate[exposure] = np.array(np.average(templist,0))
  pass

def maskedtemplateguess(exposure):
  """Initial guess for template of an exposure"""
  maskedshiftguess(exposure)
  maskedmultguess(exposure)
  templist = []
  for key in List[exposure]:
    templist.append((avcal[key] - ((maskedinitialshift[exposure][0]) * np.average(avwav[key]) + maskedinitialshift[exposure][1])) / (maskedinitialmult[exposure][0] * np.average(avwav[key]) + maskedinitialmult[exposure][1]))
  maskedinitialtemplate[exposure] = np.array(np.average(templist,0))
  pass

def plottemplateguess(template, exposure):
  """docstring for plottemplateguess"""
  print "whatever"
  for key in List[exposure]:
    pl.plot(avwav[key],avcal[key])
    pl.plot(avwav[key], initialtemplate[exposure] * (initialmult[exposure][0] * np.average(avwav[key]) + initialmult[exposure][1]) + ((initialshift[exposure][0]) * np.average(avwav[key]) + initialshift[exposure][1]) )
  pass

def newTemplate(parameterArray, exposure):
  """docstring for newTemplate"""
  templist = []
  for key in List[exposure]:
    templist.append((avcal[key] - ((parameterArray[2]) * np.average(avwav[key]) + parameterArray[3])) / (parameterArray[0] * np.average(avwav[key]) + parameterArray[1]))
  return np.array(np.average(templist,0))

def maskednewTemplate(parameterArray, exposure):
  """docstring for newTemplate"""
  templist = []
  for key in List[exposure]:
    templist.append((maskedcal[key] - ((parameterArray[2]) * np.average(avwav[key]) + parameterArray[3])) / (parameterArray[0] * np.average(avwav[key]) + parameterArray[1]))
  return ma.array(ma.average(templist,0))

def fullerrorFunction(parameterArray, exposure):
  """docstring for fullerrorFunction"""
  difference = []
  fittingtemplate = newTemplate(parameterArray, exposure)
  for key in List[exposure]:
    difference.append(avcal[key] - (fittingtemplate * (parameterArray[0] * np.average(avwav[key]) + parameterArray[1]) + ((parameterArray[2]) * np.average(avwav[key]) + parameterArray[3])))
  difference = np.concatenate(difference)
  clipdiff = difference[sigmaClippedElements(difference, sigma_clip)]
  return clipdiff / averr[exposure][sigmaClippedElements(difference, sigma_clip)]

def maskedfullerrorFunction(parameterArray, exposure):
  """docstring for fullerrorFunction"""
  difference = []
  fittingtemplate = maskednewTemplate(parameterArray, exposure)
  for key in List[exposure]:
    difference.append(maskedcal[key] - ma.masked_where(ma.getmask(maskedcal[key]), (fittingtemplate * (parameterArray[0] * np.average(avwav[key]) + parameterArray[1]) + ((parameterArray[2]) * np.average(avwav[key]) + parameterArray[3]))))
  difference = ma.concatenate(difference)
  # clipdiff = difference[sigmaClippedElements(difference, sigma_clip)]
  return difference / ma.masked_where(ma.getmask(maskedcal[key]), averr[exposure])

def sigmaClipped(inarray, sigma_clip):
  """docstring for sigmaClipped"""
  return inarray[np.where(np.abs(inarray - np.average(inarray)) < sigma_clip * np.std(inarray))]

def sigmaClippedElements(inarray, sigma_clip):
  """docstring for sigmaClipped"""
  return np.where(np.abs(inarray - np.average(inarray)) < sigma_clip * np.std(inarray))

def clippednewTemplate(parameterArray, exposure):
  templist = []
  for key in List[exposure]:
    templist.append((maskedcal[key] - ma.masked_where(ma.getmask(maskedcal[key]), ma.zeros(len(maskedcal[key])) + (parameterArray[2]) * np.average(avwav[key]) + parameterArray[3])) / (parameterArray[0] * np.average(avwav[key]) + parameterArray[1]))
  return ma.array(ma.average(templist,0))

def clippederrorFunction(parameterArray, exposure):
  difference = []
  fittingtemplate = clippednewTemplate(parameterArray, exposure)
  for key in List[exposure]:
    difference.append(maskedcal[key] - ma.masked_where(ma.getmask(maskedcal[key]),fittingtemplate * (parameterArray[0] * np.average(avwav[key]) + parameterArray[1]) + ((parameterArray[2]) * np.average(avwav[key]) + parameterArray[3])))
  bob = ma.vstack(difference)
  return ma.average(ma.abs(bob),1)

def plotdata(exposure):
  """docstring for plotdata"""
  for key in List[exposure]:  
    pl.plot(avwav[key],avcal[key])
  pass

def plotfit(parameterArray, exposure):
  """docstring for plotfit"""
  fittingtemplate =  clippednewTemplate(parameterArray, exposure)
  for key in List[exposure]:  
    pl.plot(avwav[key], fittingtemplate * (parameterArray[0] * np.average(avwav[key]) + parameterArray[1]) + ((parameterArray[2]) * np.average(avwav[key]) + parameterArray[3]))
    pl.plot(avwav[key],avcal[key])
  pl.show()
  pass

def clippedplotfit(parameterArray, exposure):
  """docstring for plotfit"""
  fittingtemplate =  clippednewTemplate(parameterArray, exposure)
  for key in List[exposure]:  
    pl.plot(avwav[key], fittingtemplate * (parameterArray[0] * np.average(avwav[key]) + parameterArray[1]) + ((parameterArray[2]) * np.average(avwav[key]) + parameterArray[3]))
    pl.plot(avwav[key],avcal[key])
  pass

def plotparametermult(parameterArray, exposure):
  """docstring for plotparameters"""
  xval = []
  yval = []
  for key in List[exposure]:
    xval.append(np.average(avwav[key]))
    yval.append(parameterArray[0] * np.average(avwav[key]) + parameterArray[1])
  pl.plot(xval,yval)
  pass

def plotmedian(exposure):
  """docstring for plotmedian"""
  xval = []
  yval = []
  for key in List[exposure]:
    xval.append(np.average(avwav[key]))
    yval.append(medianshift[exposure][0] * np.average(avwav[key]) + medianshift[exposure][1])
  pl.plot(xval,yval)  
  pass

def plotparametershift(parameterArray, exposure):
  """docstring for plotparameters"""
  xval = []
  yval = []
  for key in List[exposure]:
    xval.append(np.average(avwav[key]))
    yval.append(parameterArray[2] * np.average(avwav[key]) + parameterArray[3])
  pl.plot(xval,yval)
  pass

def plotresiduals(parameterArray, exposure):
  """docstring for plotresiduals"""
  fittingtemplate = newTemplate(parameterArray, exposure)
  for key in List[exposure]:
    pl.plot(avwav[key], avcal[key] - (fittingtemplate * (parameterArray[0] * np.average(avwav[key]) + parameterArray[1]) + ((parameterArray[2]) * np.average(avwav[key]) + parameterArray[3])))
  pass

def plotdifference(parameterArray, parameterArray2, exposure):
  """docstring for plotdifference"""
  fittingtemplate = newTemplate(parameterArray, exposure)
  fittingtemplate2 = newTemplate(parameterArray2, exposure)
  for key in List[exposure]:
    pl.plot(avwav[key], (fittingtemplate * (parameterArray[0] * np.average(avwav[key]) + parameterArray[1]) + ((parameterArray[2]) * np.average(avwav[key]) + parameterArray[3])) - (fittingtemplate2 * (parameterArray2[0] * np.average(avwav[key]) + parameterArray2[1]) + ((parameterArray2[2]) * np.average(avwav[key]) + parameterArray2[3])))  
  pass

def plotNoClip():
  """docstring for plotNoClip"""
  for exposure in i2exposures:
    plotfit(linearfit[exposure], exposure)
    pl.savefig("Plots/initial." + exposure + ".pdf")
    pl.close()
    plotfit(fit[exposure], exposure)
    pl.savefig("Plots/fit." + exposure + ".pdf")
    pl.close()
    plotresiduals(linearfit[exposure], exposure)
    pl.savefig("Plots/initialresiduals." + exposure + ".pdf")
    pl.close()
    plotresiduals(fit[exposure], exposure)
    pl.savefig("Plots/fitresiduals." + exposure + ".pdf")
    pl.close()
    plotdifference(fit[exposure],linearfit[exposure],exposure)
    pl.savefig("Plots/fitdifference." + exposure + ".pdf")
    pl.close()
    shiftmedian(exposure)
    plotdata(exposure)
    plotparametershift(linearfit[exposure],exposure)
    # plotparametershift(fit[exposure],exposure)
    plotmedian(exposure)
    tempstr = str(round(medianshift[exposure][0],4)) + "   " + str(round(initialshift[exposure][0],4)) + "   " + str(round(medianshift[exposure][0]/initialshift[exposure][0],4))
    # + "   " + str(round(fit[exposure][2],4))
    pl.title(tempstr)
    pl.savefig("Plots/slopes." + exposure + ".pdf")
    pl.close()
  pass

def residualgoodness(templateExposure, dataExposure):
  """docstring for residualgoodness"""
  fittingtemplate = newTemplate(fit[templateExposure], templateExposure)
  templist = []
  for key in List[dataExposure]:
    templist.append(avcal[key] - (fittingtemplate * (fit[templateExposure][0] * np.average(avwav[key]) + fit[templateExposure][1]) + ((fit[templateExposure][2]) * np.average(avwav[key]) + fit[templateExposure][3])))
  return np.concatenate(templist)

def initialresidualgoodness(templateExposure, dataExposure):
  """docstring for residualgoodness"""
  fittingtemplate = newTemplate(linearfit[templateExposure], templateExposure)
  templist = []
  for key in List[dataExposure]:
    templist.append(avcal[key] - (fittingtemplate * (linearfit[templateExposure][0] * np.average(avwav[key]) + linearfit[templateExposure][1]) + ((linearfit[templateExposure][2]) * np.average(avwav[key]) + linearfit[templateExposure][3])))
  return np.concatenate(templist)
  
# TODO Don't make 13 and 50 magic numbers
initialshift = {}
maskedinitialshift = {}
initialmult = {}
maskedinitialmult = {}
medianshift = {}
initialtemplate = {}
maskedinitialtemplate = {}
sigma_clip = 3.0

# clipping_limit = 2.5

linearfit = {}
maskedlinearfit = {}
# This loads up the inital guesses dictionaries
fit = {}
maskedfit = {}

def first():
  """docstring for first"""
  for exposure in i2exposures:
    for key in List[exposure]:
      templateguess(exposure)
      linearfit[exposure] = [initialmult[exposure][0],initialmult[exposure][1],initialshift[exposure][0],initialshift[exposure][1]]
    fit[exposure], success = leastsq(fullerrorFunction, linearfit[exposure], args=(exposure), maxfev=10000)
  pass

def maskedfirst():
  third()
  for exposure in i2exposures:
    for key in List[exposure]:
      maskedtemplateguess(exposure)
      maskedlinearfit[exposure] = [maskedinitialmult[exposure][0],maskedinitialmult[exposure][1],maskedinitialshift[exposure][0],maskedinitialshift[exposure][1]]
    print exposure
    maskedfit[exposure], success = leastsq(maskedfullerrorFunction, maskedlinearfit[exposure], args=(exposure), maxfev=10000)
  
# what's going on with 0065???

fittemplate = {}
def second():
  """docstring for second"""
  for exposure in i2exposures:
    fittemplate[exposure] = newTemplate(fit[exposure], exposure)
  pass

maskedcal = {}
matrixcal = {}
def third():
  """docstring for third"""
  for exposure in i2exposures:
    count = 0.0
    maskedcal[exposure] = avcal[exposure].reshape(13,50).transpose()
    for index in range(len(maskedcal[exposure])):
      av = np.average(maskedcal[exposure][index])
      std = np.std(maskedcal[exposure][index])
      temp1 = ma.masked_outside(maskedcal[exposure][index], av - std * sigma_clip, av + std * sigma_clip)
      maskedcal[exposure + str(index)] = temp1
      count += ma.count_masked(maskedcal[exposure + str(index)])
    matrixcal[exposure] = ma.column_stack([maskedcal[thing] for thing in [exposure + str(index) for index in range(len(maskedcal[exposure]))]])
    for index, key in enumerate(List[exposure]):
      maskedcal[key] = matrixcal[exposure][index]
    print exposure, count
  pass

clipfit = {}
halfclipfit = {}

def fourth():
  """docstring for fourth"""
  for exposure in i2exposures:
    print exposure
    for key in List[exposure]:
      templateguess(exposure)
      linearfit[exposure] = [initialmult[exposure][0],initialmult[exposure][1],initialshift[exposure][0],initialshift[exposure][1]]
    clipfit[exposure], success = leastsq(clippederrorFunction,  linearfit[exposure], args=(exposure), maxfev=10000)
  pass


  # for index in range(50):
  #     slice.append(linregress(avwav[exposure].reshape(13,50).transpose()[index],avcal[exposure].reshape(13,50).transpose()[index])[0])
  #     sliceshift.append(linregress(avwav[exposure].reshape(13,50).transpose()[index],avcal[exposure].reshape(13,50).transpose()[index])[1])


# TODO trouble exposures
def tweak():
  """docstring for tweak"""
  tweak = clipfit['1098']
  for exposure in ['1099']:
    for key in List[exposure]:
      templateguess(exposure)
      linearfit[exposure] = [initialmult[exposure][0],initialmult[exposure][1],initialshift[exposure][0],initialshift[exposure][1]]
    clipfit[exposure], success = leastsq(clippederrorFunction, tweak, args=(exposure), maxfev=10000)
  tweak = clipfit['0068']
  for exposure in ['0069']:
    for key in List[exposure]:
      templateguess(exposure)
      linearfit[exposure] = [initialmult[exposure][0],initialmult[exposure][1],initialshift[exposure][0],initialshift[exposure][1]]
    clipfit[exposure], success = leastsq(clippederrorFunction, tweak, args=(exposure), maxfev=10000)
  pass

def plotclippedfit():
  """docstring for plotclippedfit"""
  for exposure in i2exposures:
    clippedplotfit(clipfit[exposure], exposure)
    pl.savefig("Plots/clipped." + exposure + ".pdf")
    pl.close()
  pass

# Analysis

def clippedresidualgoodness(templateExposure, dataExposure, parameterArray):
  """docstring for residualgoodness"""
  fittingtemplate =  clippednewTemplate(clipfit[templateExposure], templateExposure)
  templist = []
  for key in List[dataExposure]:
    templist.append(maskedcal[key] - ma.masked_where(ma.getmask(maskedcal[key]),fittingtemplate * (parameterArray[0] * np.average(avwav[key]) + parameterArray[1]) + ((parameterArray[2]) * np.average(avwav[key]) + parameterArray[3])))
  return ma.concatenate(templist)

def clippedprint():
  print "Expo & Self ", 
  for exposure in i2exposures:
    print " & ", exposure, 
  print "\\\\ \\hline"
  print " \\hline"
  for exposure in i2exposures:
    linerow =     [ma.std(clippedresidualgoodness(exposure,secondexposure,clipfit[exposure])) for secondexposure in i2exposures]
    stringrow = ''
    for item in linerow:
      stringrow +=  " & " + str(np.round(item, 1))
    print exposure, " & ", np.round(ma.std(maskedcal[exposure]),1), stringrow,
    print "\\\\ \\hline"
  pass

def clippedplotdata(exposure):
  for key in List[exposure]:
    pl.plot(avwav[key],maskedcal[key], color="black")
  pl.title("Exposure: " + exposure)
  pl.xlabel("Wavelength in Angstroms")
  pl.ylabel("Calibration Shift in m/s")
  pl.savefig("Plots/data" + exposure + ".pdf")
  pl.close()
  pass

def clippedslope(parameterArray, exposure):
  xval = []
  yval = []
  for key in List[exposure]:
    pl.plot(avwav[key],maskedcal[key], color="black")
    xval.append(np.average(avwav[key]))
    yval.append(parameterArray[2] * np.average(avwav[key]) + parameterArray[3])
  pl.plot(xval,yval,color="green")  
  pl.title("Exposure: " + exposure + " Slope: " + str(round(clipfit[exposure][2],3)) + " m/s per Angstrom")
  pl.xlabel("Wavelength in Angstroms")
  pl.ylabel("Calibration Shift in m/s")
  pl.vlines(5025.,ma.average(maskedcal[exposure]) - 650, ma.average(maskedcal[exposure])+650,linewidth=2)
  pl.savefig("Plots/slope" + exposure + ".pdf")
  pl.close()
  pass

def clippedfunction(function):
  """docstring for clippedfunction"""
  for exposure in i2exposures:
    function(clipfit[exposure], exposure)
  pass


def plotclippedresiduals(parameterArray, exposure, save=False):
  fittingtemplate = newTemplate(parameterArray, exposure)
  for key in List[exposure]:
    pl.plot(avwav[key], maskedcal[key] - (fittingtemplate * (parameterArray[0] * np.average(avwav[key]) + parameterArray[1]) + ((parameterArray[2]) * np.average(avwav[key]) + parameterArray[3])),color="red")
  pl.plot(np.linspace(pl.xlim()[0],pl.xlim()[1]), np.zeros(len(np.linspace(pl.xlim()[0],pl.xlim()[1]))))
  pl.title("Residuals for Exposure: " + exposure)
  pl.xlabel("Wavelength in Angstroms")
  pl.ylabel("Calibration Shift in m/s")
  if save == True:
    pl.savefig("Plots/residuals" + exposure + ".pdf")
    pl.close()
  pass

def plotapplytemplate(templateExposure, dataExposure, parameterArray, save=False):
  """docstring for plotapplytemplate"""
  fittingtemplate =  clippednewTemplate(clipfit[templateExposure], templateExposure)
  for key in List[dataExposure]:
    pl.plot(avwav[key], maskedcal[key] - ma.masked_where(ma.getmask(maskedcal[key]),fittingtemplate * (parameterArray[0] * np.average(avwav[key]) + parameterArray[1]) + ((parameterArray[2]) * np.average(avwav[key]) + parameterArray[3])),color="red")
    pl.plot(avwav[key],maskedcal[key], color="black")
  pl.plot(np.linspace(pl.xlim()[0],pl.xlim()[1]), np.zeros(len(np.linspace(pl.xlim()[0],pl.xlim()[1]))))
  pl.title("Templ: " + templateExposure + " Expo: " + dataExposure + " std: " + str(np.round(ma.std(maskedcal[dataExposure]),1)) + " new-std: " + str(np.round(ma.std(clippedresidualgoodness(templateExposure, dataExposure, parameterArray)))))
  pl.xlabel("Wavelength in Angstroms")
  pl.ylabel("Calibration Shift in m/s")
  pl.vlines(5025.,ma.average(maskedcal[dataExposure]) - 650, ma.average(maskedcal[dataExposure])+650,linewidth=2)
  if save == True:
    pl.savefig("Plots/Template/apply" + templateExposure + dataExposure  + ".pdf")
    pl.close()
  pass

def overplot(exposure):
  for key in List[exposure]:
    pl.plot(maskedcal[key])
  pl.xlabel("Width of an Order")
  pl.ylabel("Calibration Shift in m/s")
  pass

  # def clippedresidualgoodness(templateExposure, dataExposure, parameterArray):
  #   """docstring for residualgoodness"""
  #   fittingtemplate =  clippednewTemplate(clipfit[templateExposure], templateExposure)
  #   templist = []
  #   for key in List[dataExposure]:
  #     templist.append(maskedcal[key] - ma.masked_where(ma.getmask(maskedcal[key]),fittingtemplate * (parameterArray[0] * np.average(avwav[key]) + parameterArray[1]) + ((parameterArray[2]) * np.average(avwav[key]) + parameterArray[3])))
  #   return ma.concatenate(templist)

#
# for exposure in i2exposures:
#   plotclippedresiduals(clipfit[exposure], exposure)
# for exposure in i2exposures:
#   clippedslope(clipfit[exposure], exposure)
# for exposure in i2exposures:
#   clippedplotdata(exposure)
# for exposure in i2exposures:
#   for secondexposure in i2exposures:
#     plotapplytemplate(exposure,secondexposure,clipfit[exposure])
# Example 1 -- One Order
# for exposure in ['0064.k.063']:
#   pl.plot(avwav[exposure],maskedcal[exposure],color='black')
#   pl.xlabel("Wavelength in Angstroms")
#   pl.ylabel("Calibration Shift in m/s")
#   pl.title("Exposure 0064; Order 63")
#   pl.vlines(avwav[exposure][0] - 10.5,ma.average(maskedcal[exposure]) - 650, ma.average(maskedcal[exposure])+650,linewidth=2)
#   pl.savefig('Plots/example1.pdf')
#   pl.close()
# 
# 
# # Example 2 -- One Order
# for exposure in ['0064']:
#   for key in List[exposure]:
#     pl.plot(avwav[key],maskedcal[key], color="black")    
#   pl.xlabel("Wavelength in Angstroms")
#   pl.ylabel("Calibration Shift in m/s")
#   pl.title("Exposure 0064")
#   pl.vlines(5025.,ma.average(maskedcal[exposure]) - 650, ma.average(maskedcal[exposure])+650,linewidth=2)
#   pl.savefig('Plots/example2.pdf')
#   pl.close()

# Example 3 -- Template

# TODO Analyze Monte Carlo of slopes, etc. 
# TODO include errors in fitting
# TODO sigma clip the overall number for "goodness of fit"
# TODO add all extra information -- temperature, alt, az, el, etc. 
# TODO plot cross correlations of extra information with slopes and variables
# TODO plot residuals of all cross template plots
# TODO categorize how much better the template makes things


# # Create dictionaries for wavelength, flux, continuum, etc. 
# wav = {}
# flx = {}
# err = {}
# con = {}
# pix = {}


# def errorFunction(parameterArray, exposure):
#   """docstring for errorFunction"""
#   difference = []
#   fittingtemplate = newTemplate(parameterArray, exposure)
#   for key in List[exposure]:
#     difference.append(avcal[key] - (fittingtemplate * (parameterArray[0] * np.average(avwav[key]) + parameterArray[1]) + ((parameterArray[2]) * np.average(avwav[key]) + parameterArray[3])))
#   return difference
  # return np.concatenate(difference)

# def newerrorFunction(parameterArray, exposure):
#   """docstring for errorFunction"""
#   difference = []
#   fittingtemplate = newTemplate(parameterArray, exposure)
#   for key in List[exposure]:
#     difference.append(avcal[key] - (fittingtemplate * (parameterArray[0] * np.average(avwav[key]) + parameterArray[1]) + ((parameterArray[2]) * np.average(avwav[key]) + parameterArray[3])))
#   bob = np.vstack(difference)
#   # print bob
#   clippedaverage = []
#   clippederror = []
#   for element in range(len(bob)):
#     clippedaverage.append(np.average(sigmaClipped(bob[:, element], clipping_limit) / averr[exposure][sigmaClippedElements(bob[:, element], clipping_limit)] ))
#     print sigmaClipped(bob[:,element],clipping_limit)
#   return np.array(clippedaverage)
# TODO Fix the above -- why doesn't the sigma clipping work?

# def shiftmedian(exposure):
#   """docstring for shiftmedian"""
#   temp = []
#   for key in List[exposure]:
#     temp.append([np.average(avwav[key]),np.median(avcal[key])])
#   medianshift[exposure] = linregress(temp)[0],linregress(temp)[1]
#   pass
# 
