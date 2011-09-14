#!/usr/bin/env python
# encoding: utf-8
"""
monte-firsttry.py

Created by Jonathan Whitmore on 2011-01-27.
"""

from config_calibration import *

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

for Exp in ExpList:
  # print Exp
  tempexpwav = []
  tempexpcal = []
  tempexperr = []
  tempexpres = []
  tempexpres_err = []
  tempexppix = []
  
  header[Exp] = {}
  try:
    headerfile = "Headers/" + Exp + ".header"
    headerReader = csv.reader(open(headerfile), delimiter=' ')
    for row in headerReader:
        header[Exp][row[0]] = row[1]
  except:
    print "Problem with header file for ", Exp
    
  for chip in chips:
    # print chip
    tempexpchipwav = []
    tempexpchipcal = []
    tempexpchiperr = []
    tempexpchipres = []
    tempexpchipres_err = []
    tempexpchippix = []
    ExpChip = Exp + "." + chip
    ExpChipList.append(ExpChip)
    CalibrationFileList = glob.glob("Calibration/*" + ExpChip + "*")
    for CalibrationFile in CalibrationFileList:
      ExpChipOrder = '.'.join(CalibrationFile.split("/")[-1].split(".")[1:4])
      # print ExpChipOrder, CalibrationFile
      ExpChipOrderList.append(ExpChipOrder)
      tempexpchipordwav = []
      tempexpchipordcal = []
      tempexpchiporderr = []
      tempexpchipordres = []
      tempexpchipordres_err = []
      tempexpchipordpix = []            
      for row in csv.reader(open(CalibrationFile),delimiter=' '):
        tempexpchipordwav.append(float(row[0]))
        tempexpchipordcal.append(float(row[1]))
        tempexpchiporderr.append(float(row[2]))
        tempexpchipordres.append(float(row[3]))
        tempexpchipordres_err.append(float(row[4]))
        if len(row) > 5:
          tempexpchipordpix.append(float(row[5]))
      avwav[ExpChipOrder] = np.array(tempexpchipordwav)
      avcal[ExpChipOrder] = np.array(tempexpchipordcal)
      averr[ExpChipOrder] = np.array(tempexpchiporderr)
      avres[ExpChipOrder] = np.array(tempexpchipordres)
      avres_err[ExpChipOrder] = np.array(tempexpchipordres_err)
      if len(tempexpchipordpix) > 0:
        avpix[ExpChipOrderStr] = np.array(tempexpchipordpix)
      tempexpchipwav.append(np.array(tempexpchipordwav))
      tempexpchipcal.append(np.array(tempexpchipordcal))
      tempexpchiperr.append(np.array(tempexpchiporderr))
      tempexpchipres.append(np.array(tempexpchipordres))
      tempexpchipres_err.append(np.array(tempexpchipordres_err))
      if len(tempexpchipordpix) > 0:
        tempexpchippix.append(np.array(tempexpchipordpix))
    avwav[ExpChip] = np.concatenate(tempexpchipwav)
    avcal[ExpChip] = np.concatenate(tempexpchipcal)
    averr[ExpChip] = np.concatenate(tempexpchiperr)
    avres[ExpChip] = np.concatenate(tempexpchipres)
    avres_err[ExpChip] = np.concatenate(tempexpchipres_err)
    if len(tempexpchippix) > 1:
      avpix[ExpChip] = np.concatenate(tempexpchippix)
    
    tempexpwav.append(np.concatenate(tempexpchipwav))
    tempexpcal.append(np.concatenate(tempexpchipcal))
    tempexperr.append(np.concatenate(tempexpchiperr))
    tempexpres.append(np.concatenate(tempexpchipres))
    tempexpres_err.append(np.concatenate(tempexpchipres_err))
    if len(tempexpchippix) > 1:
      tempexppix.append(np.concatenate(tempexpchippix))

  avwav[Exp] = np.concatenate(tempexpwav)
  avcal[Exp] = np.concatenate(tempexpcal)
  averr[Exp] = np.concatenate(tempexperr)
  avres[Exp] = np.concatenate(tempexpres)
  avres_err[Exp] = np.concatenate(tempexpres_err)
  if len(tempexppix) > 1:
    avpix[Exp] = np.concatenate(tempexppix)

# This creates lists of all orders in exposures
List = {}
for Exp in ExpList:
  List[Exp] = []
  for i in glob.glob("Calibration/*" + "." + Exp + "*"):
    List[Exp].append(".".join(i.split("/")[1].split(".")[1:4]))

        

QValueFile = config.get('Information','path_to_q_values')

q_reader = csv.reader(open(QValueFile), delimiter=' ')
trans = []
lwav = []
fval = []
gamma = []
mass = []
qval = []

for row in q_reader:
  trans.append(row[0])
  lwav.append(float(row[1]))
  fval.append(float(row[2]))
  gamma.append(float(row[3]))
  mass.append(float(row[4]))
  qval.append(float(row[5]))

qval = np.array(qval)    
lwav_ang = np.array(lwav)
lwav_cm = lwav_ang * 1e-8

alpha_x_axis = 2.0 * c_light * lwav_cm * qval

# redshifted = lwav * (1 + z) * 1.0e8
tolerance = 1.5

def gaussian_setup(sigma):
  tempstr = "g." + str(int(sigma))
  wav[tempstr] = np.arange(2900,10000,0.5)
  cal[tempstr] = np.random.standard_normal(len(wav[tempstr])) * sigma
  xvalues[tempstr] = []
  alphavalues[tempstr] = []
  zvalues[tempstr] = []
  slopevalues[tempstr] = []
  transxvalues[tempstr] = []
  transalphavalues[tempstr] = []
  transzvalues[tempstr] = []
  transslopevalues[tempstr] = []

zvaluesdictionary = {}
calibrationShiftsdictionary = {}
slopedictionary = {}
x_valdictionary = {}

Stars = ['0064','0065','1093','1094','2092','2093']

for exposure in ExpList:
  x_valdictionary[exposure] = []
  zvaluesdictionary[exposure] = []
  calibrationShiftsdictionary[exposure] = []
  slopedictionary[exposure] = []

avwav['flat'] = np.linspace(5000., 6000., 700.)
avcal['flat'] = np.zeros(700.)

avwav['slope+3'] = np.linspace(3000., 9000., 2000.)
avcal['slope+3'] = np.linspace(0., 2000., 2000.)

avwav['slope-3'] = np.linspace(3000., 9000., 2000.)
avcal['slope-3'] = np.linspace(0., -2000., 2000.)

Extras = ['flat', 'slope+3', 'slope-3']

Slopes = ['slope-3','slope+3']

for exposure in Extras:
  x_valdictionary[exposure] = []
  zvaluesdictionary[exposure] = []
  calibrationShiftsdictionary[exposure] = []
  slopedictionary[exposure] = []

zlow = 0.2
zhigh = 3.7

def test():
  exposure = '0064'
  mc_add(exposure, 5, 5)

def mc_add(exposure, numlines, systems, exact=False):
  """adds to memory if >= numlines if exact = False
  adds to memory == numlines if exact = True"""
  calibrationShifts = []
  x_val = []
  z_record = []
  slope =[]
  while len(calibrationShifts) < systems:
    average_vshift = []
    save_redindex = []
    z = np.random.uniform(zlow,zhigh)
    redshifted = lwav_ang * (1. + z)
    for redindex in np.where(redshifted[np.where(redshifted > np.min(avwav[exposure]))] < np.max(avwav[exposure]))[0]:
      for index, wave in enumerate(avwav[exposure]):
        if np.min(np.abs(wave - redshifted[redindex]) < tolerance):
          wavelengthindex = ra.choice(np.where(np.abs(avwav[exposure] - redshifted[redindex]) < tolerance)[0])
          average_vshift.append(avcal[exposure][wavelengthindex])
          save_redindex.append(redindex)
          break
    vshift = np.array(average_vshift)
    if len(vshift) >= numlines:
      if exact == False:
        chooseindex = range(len(vshift))
      if exact == True:
        chooseindex = ra.sample(range(len(vshift)),numlines)
      zvaluesdictionary[exposure].append(z)
      vshiftlist = [vshift[i] for i in chooseindex]
      calibrationShifts.append(np.array(vshiftlist))
      alphalist = [alpha_x_axis[save_redindex][i] for i in chooseindex]
      x_valdictionary[exposure].append(np.array(alphalist))
      slopedictionary[exposure].append(linregress(alphalist, vshiftlist)[0])
  pass

def realtimehist(exposure, numlines, systems, maxruns):
  for dummy in range(maxruns):
    mc_add(exposure, numlines, systems)
    if dummy % 100 == 0:
      print dummy
      pl.hist(slopedictionary[exposure])
      pl.title(exposure + " Av: " + str(np.round(np.average(slopedictionary[exposure]),7)) + " Std: " + str(np.round(np.std(slopedictionary[exposure]),7)) + " N: " + str(len(slopedictionary[exposure])))
      pl.savefig("Plots/Histograms/" + exposure + "." + str(numlines) + "." + str(len(slopedictionary[exposure])) + ".pdf")
      pl.close()
  pass

def mcexposures(numlines, systems, maxruns):
  """docstring for mcexposures"""
  for exposure in ExpList:
    realtimehist(exposure, numlines, systems, maxruns)
  pass

def printmc(exposure, numlines, systems):
  """docstring for printmc"""
  MCFile = "MonteCarlo/min." + str(exposure) + "." + str(numlines) + ".ascii"
  MCFileHandle = open(MCFile, 'a')
  mc_add(exposure, numlines, systems)
  print >>MCFileHandle, slopedictionary[exposure][-1], len(x_valdictionary[exposure][-1])
  pass

def transprintmc(exposure, numlines, systems):
  """docstring for printmc"""
  NtransMCFile = "MonteCarlo/exact." + str(exposure) + "." + str(numlines) + ".ascii"
  NtransMCFileHandle = open(NtransMCFile, 'a')
  mc_add(exposure, numlines, systems, exact=True)
  print >>NtransMCFileHandle, slopedictionary[exposure][-1], len(x_valdictionary[exposure][-1])
  pass

xvalues = {}
alphavalues = {}
zvalues = {}
slopevalues = {}

transxvalues = {}
transalphavalues = {}
transzvalues = {}
transslopevalues = {}

for key in avwav:
  xvalues[key] = []
  alphavalues[key] = []
  zvalues[key] = []
  slopevalues[key] = []
  transxvalues[key] = []
  transalphavalues[key] = []
  transzvalues[key] = []
  transslopevalues[key] = []

# 

sysp = {}
sysn = {}
try:
  posslope, possys = np.loadtxt('MonteCarlo/exact.slope+3.3.ascii', unpack=True)
  negslope, negsys = np.loadtxt('MonteCarlo/exact.slope-3.3.ascii', unpack=True)


def loadup():
  for systemsize in np.arange(25, 500, 25):
    sysp[systemsize] = [np.average(posslope[i:i+systemsize]) for i in xrange(0, len(posslope), systemsize)]
    sysn[systemsize] = [np.average(negslope[i:i+systemsize]) for i in xrange(0, len(negslope), systemsize)]
  pass
  
def scatterp(save=False):
  """docstring for scatter"""
  for systemsize in np.arange(25, 500, 25):
    pl.scatter(systemsize, np.std(sysp[systemsize]))
    pl.title("STD vs. Ntran")
    pl.xlabel("Ntran")
    pl.ylabel("STD in m/s")
  if save == True:
    pl.savefig('Plots/positive.slope.pdf')
    pl.close()
  pass
  
def scattern(save=False):
  """docstring for scatter"""
  for systemsize in np.arange(25, 500, 25):
    pl.scatter(systemsize, np.std(sysn[systemsize]))
    pl.title("STD vs. Ntran")
    pl.xlabel("Ntran")
    pl.ylabel("STD in m/s")
  if save == True:
    pl.savefig('Plots/negative.slope.pdf')
    pl.close()
  pass


def addafewexposures():
  for systems in [3, 4, 5, 6]:
    print systems, "Begins at: ", datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    for first in range(100):
      for exposure in Stars:
        transprintmc(exposure, systems, 1)
    print systems, "Ended at: ", datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
  pass

def addafewslopes():
  for systems in [3, 4]:
    print systems, "Begins at: ", datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    for first in range(100):
      for exposure in Slopes:
        transprintmc(exposure, systems, 1)
    print systems, "Ended at: ", datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
  pass



def Ntranstd(key,mintotalruns,Ntrans,systemsize):
    remainingruns = mintotalruns - len(transslopevalues[key])
    if remainingruns > 0:
        add_N_trans_runs(key,remainingruns,Ntrans,1)
    
    tempslope = np.concatenate(transslopevalues[key])
    temparray = [np.average(tempslope[i:i+systemsize]) for i in xrange(0,len(tempslope), systemsize)]
    print "Ntrans: ", Ntrans, "\nSystems per Experiment: ", systemsize, "\n# of Experiments: ", len(temparray)
    print "Average: ", np.average(temparray), "m/s Std: ", np.std(temparray), "m/s"
    pass 
    

def gauss_Nmin_std(mintotalruns,Nmin,systemsize,sigma):
    tempstr = "g." + str(int(sigma))
    gaussian_setup(sigma)
    add_N_min_runs(tempstr,mintotalruns,Nmin,1)
    tempslope = np.concatenate(slopevalues[tempstr])
    temparray = [np.average(tempslope[i:i+systemsize]) for i in xrange(0,len(tempslope), systemsize)]
    print "Nmin: ", Nmin, "\nSystems per Experiment: ", systemsize, "\n# of Experiments: ", len(temparray)
    print "Sigma: ", sigma, "Average: ", np.average(temparray), "m/s Std: ", np.std(temparray), "m/s"
    pass 


def ntransfunction():    
    NTranFileHandle = open(OutputMonteCarloDir + "ntran.txt",'w')
    print >>NTranFileHandle, "#sigma, Nmin, average, std"
    NTranFileHandle.close()
    Nminlist = range(2,10,2)
    SigmaList = [20,40,80,120,150,300]
    addruns = 5000
    SystemList = [143]
    for sigma in SigmaList:
        print "\n",datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S"), sigma
        NTranFileHandle = open(OutputMonteCarloDir + "ntran.txt",'a')
        for Nmin in Nminlist:
            tempstr = "g." + str(int(sigma))
            gaussian_setup(sigma)
            add_N_trans_runs(tempstr,addruns,Nmin,1)
            tempslope = np.concatenate(transslopevalues[tempstr])
            print Nmin,
            for systemsize in SystemList:
                temparray = [np.average(tempslope[i:i+systemsize]) for i in xrange(0,len(tempslope), systemsize)]
                # print "#/exp: ", systemsize, "tot #: ", len(temparray),
                # print "Sig: ", sigma, "Average: ", np.average(temparray), "m/s Std: ", np.std(temparray), "m/s"
                print >>NTranFileHandle, sigma, Nmin, np.average(temparray), np.std(temparray)
        NTranFileHandle.close()
    
def ntranssys1():    
    NTranFileHandle = open(OutputMonteCarloDir + "ntranodd.txt",'w')
    print >>NTranFileHandle, "#sigma, Nmin, average, std"
    NTranFileHandle.close()
    Nminlist = range(3,17,2)
    SigmaList = [150]
    addruns = 100000
    for sigma in SigmaList:
        for Nmin in Nminlist:
            print "\n",datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S"), sigma
            NTranFileHandle = open(OutputMonteCarloDir + "ntranodd.txt",'a')
            tempstr = "g." + str(int(sigma))
            gaussian_setup(sigma)
            add_N_trans_runs(tempstr,addruns,Nmin,1)
            tempslope = np.concatenate(transslopevalues[tempstr])
            print Nmin,
            temparray = tempslope
            print >>NTranFileHandle, sigma, Nmin, np.average(temparray), np.std(temparray)
            NTranFileHandle.close()
    
