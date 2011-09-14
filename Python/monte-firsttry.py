#!/usr/bin/env python
# encoding: utf-8
"""
monte-firsttry.py

Created by Jonathan Whitmore on 2010-06-21.
"""

from config_calibration import *

c_light = 299792458.

# Create dictionaries for wavelength, calibration, and error
wav = {}
cal = {}
err = {}
pix = {}

for core_exposure in exposures_to_analyze:
    for chip in chips:
        tempwav = []
        tempcal = []
        temperr = []
        tempstring = core_exposure + "." + chip
        CalibrationFileList = glob.glob(OutputCalibrationDir + astro_object + "." + tempstring + "*")
        for CalibrationFile in CalibrationFileList:
            for row in csv.reader(open(CalibrationFile),delimiter=' '):
                tempwav.append(float(row[0]))
                tempcal.append(float(row[1]))
                temperr.append(float(row[2]))
        
        # For exposure names assign the wav 
        wav[tempstring] = np.array(tempwav)
        cal[tempstring] = np.array(tempcal)
        err[tempstring] = np.array(temperr)
        

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


def add_monte(calibration,addruns,numlines,systems):
    """add monte carlo runs on fly keeping memory"""
    starting_length = len(calibration[2])
    while len(calibration[2]) < addruns + starting_length:
        delta_alpha = []
        x_val = []
        z_record = []
        slope = []
        # Murphy had 143 systems
        while len(delta_alpha) < systems:
            average_vshift = []
            save_index = []
            # Take list of lab wavelength values, mult by random redshift
            z = np.random.uniform(0.2,3.7)
            # print z
            redshifted = lwav_ang * (1 + z)
            # Only have calibration data in 5000 - 6200 range
            for i in np.where(redshifted[np.where(redshifted > 5000)] < 6200.)[0]:
                for j in range(len(calibration[3])):
                    if(np.min(np.abs(calibration[3][j] - redshifted[i])) < tolerance):
                        average_vshift.append(calibration[4][j][np.argmin(np.abs(calibration[3][j] - redshifted[i]))])
                        save_index.append(i)
                        break
            vshift = np.array(average_vshift)
            if(len(vshift) > numlines):
                z_record.append(z)
                delta_alpha.append( vshift ) 
                helpy = alpha_x_axis[save_index]
                x_val.append(helpy)
                slope.append(linregress(helpy,vshift)[0])
        calibration[0].append(x_val)
        calibration[1].append(delta_alpha)
        calibration[5].append(slope)
        calibration[2].append(z_record)
    pass 

xvalues = {}
alphavalues = {}
zvalues = {}
slopevalues = {}

transxvalues = {}
transalphavalues = {}
transzvalues = {}
transslopevalues = {}

for key in wav:
    xvalues[key] = []
    alphavalues[key] = []
    zvalues[key] = []
    slopevalues[key] = []
    transxvalues[key] = []
    transalphavalues[key] = []
    transzvalues[key] = []
    transslopevalues[key] = []

def add_N_min_runs(key,addruns,numlines,systems):
    """add monte carlo runs on fly keeping memory
    numlines is the minimum number of transitions
    that is accepted. numlines = 3
    3 transitions and more can be selected
    """
    # print datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    starting_length = len(zvalues[key])
    while len(zvalues[key]) < addruns + starting_length:
        delta_alpha = []
        x_val = []
        z_record = []
        slope = []
        while len(delta_alpha) < systems:
            average_vshift = []
            save_index = []
            # Take list of lab wavelength values, mult by random redshift
            z = np.random.uniform(0.2,3.7)
            redshiftarray = lwav_ang * (1.0 + z)
            for redshift in np.where(redshiftarray[np.where(redshiftarray > np.min(wav[key]))] < np.max(wav[key]))[0]:
                if len(np.where(np.abs(wav[key] - redshiftarray[redshift]) < tolerance)[0]) > 0:
                    wavelengthindex = ra.choice(np.where(np.abs(wav[key] - redshiftarray[redshift]) < tolerance)[0])
                    average_vshift.append(cal[key][wavelengthindex])
                    save_index.append(redshift)
            vshift = np.array(average_vshift)
            # print vshift, len(vshift)
            if(len(vshift) >= numlines):
                z_record.append(z)
                delta_alpha.append(vshift) 
                helpy = alpha_x_axis[save_index]
                x_val.append(helpy)
                slope.append(linregress(helpy,vshift)[0])
        xvalues[key].append(x_val)
        alphavalues[key].append(delta_alpha)
        slopevalues[key].append(slope)
        zvalues[key].append(z_record)
    # print datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")    
    pass 

def add_N_trans_runs(key,addruns,numlines,systems):
    """add monte carlo runs on fly keeping memory
    numlines is the exact number of transitions
    that is accepted. numlines = 3
    3 transitions will be selected
    """
    # print datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    starting_length = len(transzvalues[key])
    while len(transzvalues[key]) < addruns + starting_length:
        delta_alpha = []
        x_val = []
        z_record = []
        slope = []
        while len(delta_alpha) < systems:
            average_vshift = []
            save_index = []
            # Take list of lab wavelength values, mult by random redshift
            z = np.random.uniform(0.2,3.7)
            redshiftarray = lwav_ang * (1.0 + z)
            for redshift in np.where(redshiftarray[np.where(redshiftarray > np.min(wav[key]))] < np.max(wav[key]))[0]:
                if len(np.where(np.abs(wav[key] - redshiftarray[redshift]) < tolerance)[0]) > 0:
                    wavelengthindex = ra.choice(np.where(np.abs(wav[key] - redshiftarray[redshift]) < tolerance)[0])
                    average_vshift.append(cal[key][wavelengthindex])
                    save_index.append(redshift)
            vshift = np.array(average_vshift)
            # print vshift, len(vshift)
            if(len(vshift) >= numlines):
                chosenindex = ra.sample(range(len(vshift)),numlines)
                z_record.append(z)
                delta_alpha.append(vshift[chosenindex]) 
                helpy = alpha_x_axis[save_index][chosenindex]
                x_val.append(helpy)
                slope.append(linregress(helpy,vshift[chosenindex])[0])
        transxvalues[key].append(x_val)
        transalphavalues[key].append(delta_alpha)
        transslopevalues[key].append(slope)
        transzvalues[key].append(z_record)
    # print datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
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
    
