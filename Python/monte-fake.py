#!/usr/bin/env python
# encoding: utf-8
"""
monte.py

Created by Jonathan Whitmore on 2009-11-16
Updated: 2009-12-14
"""
# Needs load.py 

import time
import datetime
import random
from load import *

rancal1 = cal1.copy()
rancal2 = cal2.copy()
rancal3 = cal3.copy()
ra.shuffle(rancal1)
ra.shuffle(rancal2)
ra.shuffle(rancal3)

random1 = (rancal1 - cal1)
random2 = (rancal2 - cal2)
random3 = (rancal3 - cal3)

q_values_list = 'qvalues.txt'
q_reader = csv.reader(open(q_values_list), delimiter=' ')
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

alpha_x_axis = 2.0 * c * lwav_cm * qval
tolerance = 0.2

lower_bound = 3000.0 # Lower limit of VLT

fakewav1 = []
for i in expo_wav1:
    fakewav1.append( i / ( expo_wav1[0][0] / lower_bound) )
for counter in range(5):
    boo = len(fakewav1)
    print counter, fakewav1[boo-1][-1]
    for i in expo_wav1:
        fakewav1.append( i / ( expo_wav1[0][0] / fakewav1[boo-1][-1] ) )

fakewav2 = []
for i in expo_wav2:
    fakewav2.append( i / ( expo_wav2[0][0] / lower_bound) )
for counter in range(5):
    boo = len(fakewav2)
    print counter, fakewav2[boo-1][-1]
    for i in expo_wav2:
        fakewav2.append( i / ( expo_wav2[0][0] / fakewav2[boo-1][-1] ) )

fakewav3 = []
for i in expo_wav3:
    fakewav3.append( i / ( expo_wav3[0][0] / lower_bound) )
for counter in range(5):
    boo = len(fakewav3)
    print counter, fakewav3[boo-1][-1]
    for i in expo_wav3:
        fakewav3.append( i / ( expo_wav3[0][0] / fakewav3[boo-1][-1] ) )
        # print fakewav3[counter-1][-1], fakewav3[counter][-1], fakewav3[counter+1][-1]

fakecal1 = expo_cal1 * 6
fakecal2 = expo_cal2 * 6
fakecal3 = expo_cal3 * 6

# wf1 = []
# wf2 = []
# wf3 = []
# for i in xrange(len(wav)):
#     wf1.append(wav[i] / ( wav1[0] / lower_bound ) )
# wf2 = []
# for i in xrange(len(wav)):
#     wf2.append(wav[i] / ( wav1[0] / wf1[-1][-1] ) )
# wf3 = []
# for i in xrange(len(wav)):
#     wf3.append(wav[i] / ( wav1[0] / wf2[-1][-1] ) )
# wf4 = []
# for i in xrange(len(wav)):
#     wf4.append(wav[i] / ( wav1[0] / wf3[-1][-1] ) )
# wf5 = []
# for i in xrange(len(wav)):
#     wf5.append(wav[i] / ( wav1[0] / wf4[-1][-1] ) )
# wf6 = []
# for i in xrange(len(wav)):
#     wf6.append(wav[i] / ( wav1[0] / wf5[-1][-1] ) )
# 
# falsewav1 = wf1 + wf2 + wf3 + wf4 + wf5 + wf6
# falsecal1 = cal + cal + cal + cal + cal + cal

expo_wtest1 = []
expo_wtest2 = []
expo_wtest3 = []
expo_ctest1 = []
expo_ctest2 = []
expo_ctest3 = []
expo_fasttest = []

for i in range(len(fakewav1)):
    expo_wtest1.append(fakewav1[i].copy())
    expo_wtest2.append(fakewav2[i].copy())
    expo_wtest3.append(fakewav3[i].copy())
    expo_ctest1.append(fakecal1[i].copy())
    expo_ctest2.append(fakecal2[i].copy())
    expo_ctest3.append(fakecal3[i].copy())
    
    # expo_test2.append(expo_wav1[i].copy())
    # expo_test3.append(expo_wav1[i].copy())
    # expo_test4.append(expo_wav1[i].copy())
            
runs = 200
systems = 143

monte_sin_1 = []
monte_sin_2 = []
monte_test_1 = []
monte_test_2 = []
monte_test_3 = []
monte_test_4 = []
monte_expo_1 = []
monte_expo_2 = []
monte_expo_3 = []
monte_fast_test1 = []
monte_fast_test2 = []
monte_fast_test3 = []

for i in range(6):
    # 0 # 2 q c lambda
    # 1 # velocity
    # 2 # z's used
    # 3 # wavelength
    # 4 # calibration
    # 5 # slopes
    monte_test_1.append([])
    monte_test_2.append([])
    monte_test_3.append([])
    # monte_test_4.append([])
    monte_expo_1.append([])
    monte_expo_2.append([])
    monte_expo_3.append([])
    monte_fast_test1.append([])
    monte_fast_test2.append([])
    monte_fast_test3.append([])

# monte_sin_1[3] = expo_wav1
# monte_sin_1[4] = expo_sin1
# monte_sin_2[3] = expo_wav1
# monte_sin_2[4] = expo_sin2
monte_test_1[3] = expo_wtest1
monte_test_2[3] = expo_wtest2
monte_test_3[3] = expo_wtest3
monte_test_1[4] = expo_ctest1
monte_test_2[4] = expo_ctest2
monte_test_3[4] = expo_ctest3
monte_expo_1[3] = expo_wav1
monte_expo_1[4] = expo_cal1
monte_expo_2[3] = expo_wav2
monte_expo_2[4] = expo_cal2
monte_expo_3[3] = expo_wav3
monte_expo_3[4] = expo_cal3

# Sin
monte_fast_test1[3] = np.arange(2900,10000,0.2)
monte_fast_test1[4] = np.sin(monte_fast_test1[3] * np.pi / 150) * 250.0

# Smaller Sin
monte_fast_test2[3] = np.arange(2900,10000,0.2)
monte_fast_test2[4] = np.sin(monte_fast_test2[3] * np.pi / 150) * 84.2 * 2 

# Gaussian 
monte_fast_test3[3] = np.arange(2900,10000,0.2)
monte_fast_test3[4] = np.random.standard_normal(len(monte_fast_test3[3])) * 84.2 * 2 


def clear_monte(calibration,addruns,numlines,systems):
    """add monte carlo runs on fly clearing memory"""
    calibration[0] = []
    calibration[1] = []
    calibration[2] = []
    calibration[5] = []
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
            redshifted = lwav_ang * (1 + z)
            # Only have calibration data in 5000 - 6200 range
            for i in np.where(redshifted[np.where(redshifted > calibration[3][0][0])] < calibration[3][-1][-1])[0]:
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
            for i in np.where(redshifted[np.where(redshifted > calibration[3][0][0])] < calibration[3][-1][-1])[0]:
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


def clear_test_monte(calibration,addruns,numlines,systems):
    """add monte carlo runs on fly clearing memory"""
    calibration[0] = []
    calibration[1] = []
    calibration[2] = []
    calibration[5] = []
    starting_length = len(calibration[2])
    while len(calibration[2]) < addruns + starting_length:
        delta_alpha = []
        x_val = []
        z_record = []
        slope = []
        while len(delta_alpha) < systems:
            average_vshift = []
            save_index = []
            z = np.random.uniform(0.2,3.7)
            redshifted = lwav_ang * (1 + z)
            for i in np.where(redshifted[np.where(redshifted > calibration[3][0])] < calibration[3][-1])[0]:
                if(np.min(np.abs(calibration[3] - redshifted[i])) < tolerance):
                    average_vshift.append(calibration[4][np.argmin(np.abs(calibration[3] - redshifted[i]))])
                    save_index.append(i)
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


def add_test_monte(calibration,addruns,numlines,systems):
    """add monte carlo runs on fly clearing memory"""
    starting_length = len(calibration[2])
    while len(calibration[2]) < addruns + starting_length:
        delta_alpha = []
        x_val = []
        z_record = []
        slope = []
        while len(delta_alpha) < systems:
            average_vshift = []
            save_index = []
            z = np.random.uniform(0.2,3.7)
            redshifted = lwav_ang * (1 + z)
            for i in np.where(redshifted[np.where(redshifted> calibration[3][0])] < calibration[3][-1])[0]:
                if(np.min(np.abs(calibration[3] - redshifted[i])) < tolerance):
                    average_vshift.append(calibration[4][np.argmin(np.abs(calibration[3] - redshifted[i]))])
                    save_index.append(i)
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


# for i in xrange(10):
#     print i
#     add_monte(monte_test_1,10,1,10)
for i in xrange(10):
    print i
    add_test_monte(monte_fast_test1,10,1,200)
    add_test_monte(monte_fast_test2,10,1,200)
    add_test_monte(monte_fast_test3,10,1,200)


for i in xrange(10):
    print i
    add_monte(monte_test_2,10,1,200)
    

for i in xrange(10):
    print i
    add_monte(monte_test_3,10,1,200)
    
# def full_output_test():
#     """record z, n transitions, da/a, everything in monte carlo"""
    
    
    
# 
# 
# slope_test1 = [[],[],[]]
# slope_expo1 = [[],[],[]]
# slope_expo2 = [[],[],[]]
# slope_expo3 = [[],[],[]]
# zval_test1 = [[],[],[]]
# zval_expo1 = [[],[],[]]
# zval_expo2 = [[],[],[]]
# zval_expo3 = [[],[],[]]
# 
# # slope_test1.append(np.concatenate(monte_fast_test1[5][:]))
# # slope_expo1.append(np.concatenate(monte_expo_1[5][:]))
# # slope_expo2.append(np.concatenate(monte_expo_2[5][:]))
# # slope_expo3.append(np.concatenate(monte_expo_3[5][:]))
# 
# addruns = 100
# systems = 100
# 
# # First
# for blah in range(50):
#     print blah
#     for i in [0,1,2]:
#         print "expo1", datetime.datetime.now()
#         clear_monte(monte_expo_1,addruns,i*2 + 1,systems)
#         slope_expo1[i].append(np.concatenate(monte_expo_1[5][:]))
#         zval_expo1[i].append(np.concatenate(monte_expo_1[2][:]))
#         
#         
#         print "expo2", datetime.datetime.now()
#         clear_monte(monte_expo_2,addruns,i*2 + 1,systems)
#         slope_expo2[i].append(np.concatenate(monte_expo_2[5][:]))
#         zval_expo2[i].append(np.concatenate(monte_expo_2[2][:]))
# 
# 
#         print "expo3", datetime.datetime.now()
#         clear_monte(monte_expo_3,addruns,i*2 + 1,systems)
#         slope_expo3[i].append(np.concatenate(monte_expo_3[5][:]))
#         zval_expo3[i].append(np.concatenate(monte_expo_3[2][:]))
# 
# 
#         print "test1", datetime.datetime.now()
#         clear_test_monte(monte_fast_test1,addruns,i*2 + 1,systems)
#         slope_test1[i].append(np.concatenate(monte_fast_test1[5][:]))
#         zval_test1[i].append(np.concatenate(monte_fast_test1[2][:]))
# 
# 
# 
# # Second
# slope_expo1[0] = np.concatenate(slope_expo1[0][:])
# slope_expo2[0] = np.concatenate(slope_expo2[0][:])
# slope_expo3[0] = np.concatenate(slope_expo3[0][:])
# slope_test1[0] = np.concatenate(slope_test1[0][:])
# slope_expo1[1] = np.concatenate(slope_expo1[1][:])
# slope_expo2[1] = np.concatenate(slope_expo2[1][:])
# slope_expo3[1] = np.concatenate(slope_expo3[1][:])
# slope_test1[1] = np.concatenate(slope_test1[1][:])
# slope_expo1[2] = np.concatenate(slope_expo1[2][:])
# slope_expo2[2] = np.concatenate(slope_expo2[2][:])
# slope_expo3[2] = np.concatenate(slope_expo3[2][:])
# slope_test1[2] = np.concatenate(slope_test1[2][:])
# 
# np.savetxt('ex1slope0.ascii',slope_expo1[0])
# np.savetxt('ex1slope1.ascii',slope_expo1[1])
# np.savetxt('ex1slope2.ascii',slope_expo1[2])  
# np.savetxt('ex2slope0.ascii',slope_expo2[0])
# np.savetxt('ex2slope1.ascii',slope_expo2[1])
# np.savetxt('ex2slope2.ascii',slope_expo2[2])                                   
# np.savetxt('ex3slope0.ascii',slope_expo3[0])
# np.savetxt('ex3slope1.ascii',slope_expo3[1])
# np.savetxt('ex3slope2.ascii',slope_expo3[2])
# np.savetxt('test1slope0.ascii',slope_test1[0])
# np.savetxt('test1slope1.ascii',slope_test1[1])
# np.savetxt('test1slope2.ascii',slope_test1[2])
# 
# zval_expo1[0] = np.concatenate(zval_expo1[0][:])
# zval_expo2[0] = np.concatenate(zval_expo2[0][:])
# zval_expo3[0] = np.concatenate(zval_expo3[0][:])
# zval_test1[0] = np.concatenate(zval_test1[0][:])
# zval_expo1[1] = np.concatenate(zval_expo1[1][:])
# zval_expo2[1] = np.concatenate(zval_expo2[1][:])
# zval_expo3[1] = np.concatenate(zval_expo3[1][:])
# zval_test1[1] = np.concatenate(zval_test1[1][:])
# zval_expo1[2] = np.concatenate(zval_expo1[2][:])
# zval_expo2[2] = np.concatenate(zval_expo2[2][:])
# zval_expo3[2] = np.concatenate(zval_expo3[2][:])
# zval_test1[2] = np.concatenate(zval_test1[2][:])
# 
# np.savetxt('ex1zval0.ascii',zval_expo1[0])
# np.savetxt('ex1zval1.ascii',zval_expo1[1])
# np.savetxt('ex1zval2.ascii',zval_expo1[2])  
# np.savetxt('ex2zval0.ascii',zval_expo2[0])
# np.savetxt('ex2zval1.ascii',zval_expo2[1])
# np.savetxt('ex2zval2.ascii',zval_expo2[2])                                   
# np.savetxt('ex3zval0.ascii',zval_expo3[0])
# np.savetxt('ex3zval1.ascii',zval_expo3[1])
# np.savetxt('ex3zval2.ascii',zval_expo3[2])
# np.savetxt('test1zval0.ascii',zval_test1[0])
# np.savetxt('test1zval1.ascii',zval_test1[1])
# np.savetxt('test1zval2.ascii',zval_test1[2])
# 
# 
# expo1_by_n = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
# expo1_by_z = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
# expo2_by_n = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
# expo2_by_z = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
# expo3_by_n = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
# expo3_by_z = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
# test1_by_n = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
# test1_by_z = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
# for i in range(len(monte_expo_1[0])):
#     for j in range(len(monte_expo_1[0][i])):
#         expo1_by_n[len(monte_expo_1[0][i][j])].append(monte_expo_1[5][i][j])
#         expo1_by_z[len(monte_expo_1[0][i][j])].append(monte_expo_1[2][i][j])
# for i in range(len(monte_expo_2[0])):
#     for j in range(len(monte_expo_2[0][i])):
#         expo2_by_n[len(monte_expo_2[0][i][j])].append(monte_expo_2[5][i][j])
#         expo2_by_z[len(monte_expo_2[0][i][j])].append(monte_expo_2[2][i][j])
# for i in range(len(monte_expo_3[0])):
#     for j in range(len(monte_expo_3[0][i])):
#         expo3_by_n[len(monte_expo_3[0][i][j])].append(monte_expo_3[5][i][j])
#         expo3_by_z[len(monte_expo_3[0][i][j])].append(monte_expo_3[2][i][j])
# for i in range(len(monte_fast_test1[0])):
#     for j in range(len(monte_fast_test1[0][i])):
#         test1_by_n[len(monte_fast_test1[0][i][j])].append(monte_fast_test1[5][i][j])
#         test1_by_z[len(monte_fast_test1[0][i][j])].append(monte_fast_test1[2][i][j])
# 
# 



# fast
counter = 2
for i in [two_slope,three_slope,four_slope,five_slope,six_slope]:
    print counter, np.average(i), np.std(i) / 2
    counter = counter + 1


# regular
counter = 2
for i in [two_slope,three_slope,four_slope,five_slope,six_slope]:
    print counter, np.average(i), np.std(i) 
    counter = counter + 1

for i in [two_slope2,three_slope2,four_slope2,five_slope2,six_slope2]:
    print "2", np.average(i), np.std(i)

for i in [two_slope3,three_slope3,four_slope3,five_slope3,six_slope3]:
    print "3", np.average(i), np.std(i)



# Expo1
# Take minimum number of transitions
slope = []
for mintrans in range(2,10):
    slope.append([])
    for i in range(len(monte_expo_1[0])):
        for j in range(len(monte_expo_1[0][i])):
            if(len(monte_expo_1[0][i][j]) >= mintrans):
                indices = random.sample(range(len(monte_expo_1[0][i][j])), mintrans)
                x_temp = [monte_expo_1[0][i][j][x] for x in indices]
                y_temp = [monte_expo_1[1][i][j][x] for x in indices]
                slope[mintrans-2].append(linregress(x_temp,y_temp)[0])
    print mintrans, np.average(slope[mintrans-2]), np.std(slope[mintrans-2]) 


# Repeated Expo1
# Take minimum number of transitions
slope = []
for mintrans in range(2,24):
    slope.append([])
    for i in range(len(monte_test_1[0])):
        for j in range(len(monte_test_1[0][i])):
            if(len(monte_test_1[0][i][j]) >= mintrans):
                indices = random.sample(range(len(monte_test_1[0][i][j])), mintrans)
                x_temp = [monte_test_1[0][i][j][x] for x in indices]
                y_temp = [monte_test_1[1][i][j][x] for x in indices]
                slope[mintrans-2].append(linregress(x_temp,y_temp)[0])
    print mintrans, np.average(slope[mintrans-2]), np.std(slope[mintrans-2]) 

# Repeated Expo2
# Take minimum number of transitions
slope = []
for mintrans in range(2,24):
    slope.append([])
    for i in range(len(monte_test_2[0])):
        for j in range(len(monte_test_2[0][i])):
            if(len(monte_test_2[0][i][j]) >= mintrans):
                indices = random.sample(range(len(monte_test_2[0][i][j])), mintrans)
                x_temp = [monte_test_2[0][i][j][x] for x in indices]
                y_temp = [monte_test_2[1][i][j][x] for x in indices]
                slope[mintrans-2].append(linregress(x_temp,y_temp)[0])
    print mintrans, np.average(slope[mintrans-2]), np.std(slope[mintrans-2]) 


# Repeated Expo3
# Take minimum number of transitions
slope = []
for mintrans in range(2,24):
    slope.append([])
    for i in range(len(monte_test_3[0])):
        for j in range(len(monte_test_3[0][i])):
            if(len(monte_test_3[0][i][j]) >= mintrans):
                indices = random.sample(range(len(monte_test_3[0][i][j])), mintrans)
                x_temp = [monte_test_3[0][i][j][x] for x in indices]
                y_temp = [monte_test_3[1][i][j][x] for x in indices]
                slope[mintrans-2].append(linregress(x_temp,y_temp)[0])
    print mintrans, np.average(slope[mintrans-2]), np.std(slope[mintrans-2]) 



# Scale errors
# Fast test 1 
# Take minimum number of transitions
slope1 = []
for mintrans in range(2,24):
    slope1.append([])
    for i in range(len(monte_fast_test1[0])):
        for j in range(len(monte_fast_test1[0][i])):
            if(len(monte_fast_test1[0][i][j]) >= mintrans):
                indices = random.sample(range(len(monte_fast_test1[0][i][j])), mintrans)
                x_temp = [monte_fast_test1[0][i][j][x] for x in indices]
                y_temp = [monte_fast_test1[1][i][j][x] for x in indices]
                slope1[mintrans-2].append(linregress(x_temp,y_temp)[0])
    print mintrans, np.average(slope1[mintrans-2]), np.std(slope1[mintrans-2]) 

# Test 2 -- errors
# Fast test 2
# Take minimum number of transitions
slope2 = []
for mintrans in range(2,24):
    slope2.append([])
    for i in range(len(monte_fast_test2[0])):
        for j in range(len(monte_fast_test2[0][i])):
            if(len(monte_fast_test2[0][i][j]) >= mintrans):
                indices = random.sample(range(len(monte_fast_test2[0][i][j])), mintrans)
                x_temp = [monte_fast_test2[0][i][j][x] for x in indices]
                y_temp = [monte_fast_test2[1][i][j][x] for x in indices]
                slope2[mintrans-2].append(linregress(x_temp,y_temp)[0])
    print mintrans, np.average(slope2[mintrans-2]), np.std(slope2[mintrans-2]) / 2.0 * ( np.pi / 2. )

# Test 3 -- errors
# Fast test 3
slope3 = []
for mintrans in range(2,24):
    slope3.append([])
    for i in range(len(monte_fast_test3[0])):
        for j in range(len(monte_fast_test3[0][i])):
            if(len(monte_fast_test3[0][i][j]) >= mintrans):
                indices = random.sample(range(len(monte_fast_test3[0][i][j])), mintrans)
                x_temp = [monte_fast_test3[0][i][j][x] for x in indices]
                y_temp = [monte_fast_test3[1][i][j][x] for x in indices]
                slope3[mintrans-2].append(linregress(x_temp,y_temp)[0])
    print mintrans, np.average(slope3[mintrans-2]), np.std(slope3[mintrans-2]) / 2.0 




# for i in [two_slope,three_slope,four_slope,five_slope,six_slope]:
#     print len(i)

two_tally = 0
three_tally = 0 
four_tally = 0
five_tally = 0
six_tally = 0
seven_tally = 0
eight_tally = 0
nine_tally = 0
ten_tally = 0
eleven_tally = 0 
twelve_tally = 0 
thirteen_tally = 0 
fourteen_tally = 0 
fifteen_tally = 0 
sixteen_tally = 0 
seventeen_tally = 0 
eighteen_tally = 0 
nineteen_tally = 0 
twenty_tally = 0 
twentyone_tally = 0 
twentytwo_tally = 0 
twentythree_tally = 0 
for i in range(len(monte_test_1[0])):
    for j in range(len(monte_test_1[0][i])):
        if(len(monte_test_1[0][i][j]) == 2 ):
            two_tally = two_tally + 1
        if(len(monte_test_1[0][i][j]) == 3 ):
            three_tally = three_tally + 1
        if(len(monte_test_1[0][i][j]) == 4 ):
            four_tally = four_tally + 1
        if(len(monte_test_1[0][i][j]) == 5 ):
            five_tally = five_tally + 1
        if(len(monte_test_1[0][i][j]) == 6 ):
            six_tally = six_tally + 1
        if(len(monte_test_1[0][i][j]) == 7 ):
            seven_tally = seven_tally + 1
        if(len(monte_test_1[0][i][j]) == 8 ):
            eight_tally = eight_tally + 1
        if(len(monte_test_1[0][i][j]) == 9 ):
            nine_tally = nine_tally + 1
        if(len(monte_test_1[0][i][j]) == 10 ):
            ten_tally = ten_tally + 1
        if(len(monte_test_1[0][i][j]) == 11 ):
            eleven_tally = eleven_tally + 1
        if(len(monte_test_1[0][i][j]) == 12 ):
            twelve_tally = twelve_tally + 1
        if(len(monte_test_1[0][i][j]) == 13 ):
            thirteen_tally = thirteen_tally + 1 
        if(len(monte_test_1[0][i][j]) == 14 ):
            fourteen_tally = fourteen_tally + 1
        if(len(monte_test_1[0][i][j]) == 15 ):
            fifteen_tally = fifteen_tally + 1
        if(len(monte_test_1[0][i][j]) == 16 ):
            sixteen_tally = sixteen_tally + 1
        if(len(monte_test_1[0][i][j]) == 17 ):
            seventeen_tally = seventeen_tally + 1 
        if(len(monte_test_1[0][i][j]) == 18 ):
            eighteen_tally = eighteen_tally + 1
        if(len(monte_test_1[0][i][j]) == 19 ):
            nineteen_tally = nineteen_tally + 1
        if(len(monte_test_1[0][i][j]) == 20 ):
            twenty_tally = twenty_tally + 1
        if(len(monte_test_1[0][i][j]) == 21 ):
            twentyone_tally = twentyone_tally + 1 
        if(len(monte_test_1[0][i][j]) == 22 ):
            twentytwo_tally = twentytwo_tally + 1
        if(len(monte_test_1[0][i][j]) == 23 ):
            twentythree_tally = twentythree_tally + 1 
        

ztally = [[],[]]
hmm = 2
for k in range(23):
    ztally[0].append(hmm)
    ztally[1].append([])
    for i in range(len(monte_test_2[0])):
        for j in range(len(monte_test_2[0][i])):
            if( len(monte_test_2[0][i][j]) == hmm):
                ztally[1][k].append(monte_test_2[2][i][j])
    hmm = hmm + 1
 
pl.plot(ztally[0],[len(ztally[1][x]) for x in xrange(23)])

Ntrans = [[],[]]
hmm = 2
for k in range(23):
    Ntrans[0].append(hmm)
    Ntrans[1].append([])
    for i in range(len(monte_test_3[0])):
        for j in range(len(monte_test_3[0][i])):
            if(len(monte_test_3[0][i][j]) == hmm):
                Ntrans[1][k].append(monte_test_3[5][i][j])
    hmm = hmm + 1

for i in range(22):
    print i+2, np.average(Ntrans[1][i]), np.std(Ntrans[1][i])


# 
# hmm = 2 
# tallyx = []
# tallyy = []       
# for i in [ two_tally, three_tally, four_tally, five_tally, six_tally, seven_tally, eight_tally, nine_tally, ten_tally, eleven_tally, twelve_tally, thirteen_tally, fourteen_tally, fifteen_tally, sixteen_tally, seventeen_tally, eighteen_tally, nineteen_tally, twenty_tally, twentyone_tally, twentytwo_tally, twentythree_tally ]:
#     tallyx.append(hmm)
#     tallyy.append(i)
#     print hmm, i
#     hmm = hmm + 1
# 
# 
# 
# 
# 
# 
# 
# 
# for system_size in [50,143]:
#     k=0
#     for j in [slope_expo1,slope_expo2,slope_expo3,slope_test1]:
#         for i in range(3):
#             if( i == 0 ):
#                 if( k == 0 ):
#                     print "Exposure 1"
#                 if( k == 1 ):
#                     print "Exposure 2"
#                 if( k == 2 ):
#                     print "Exposure 3"
#                 if( k == 3 ):
#                     print "Test Case"
#                 print "N_min,\t mean,\t std,\t System Size, std/sqrt, std_of_mean"
#             print 2 * i + 2, np.average(j[i]), np.std([j[i]]), system_size, np.std([j[i]])/ np.sqrt(system_size), np.std([np.average(j[i][internal:internal+143]) for internal in xrange(0,len(j[i]),system_size)])
#         k = k + 1
# 

