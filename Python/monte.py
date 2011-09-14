#!/usr/bin/env python
# encoding: utf-8
"""
monte.py

Created by Jonathan Whitmore on 2009-11-16
Updated: 2009-12-14
"""
# Needs load.py 

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

# redshifted = lwav * (1 + z) * 1.0e8
tolerance = 0.2
expo_sin1 = []
expo_sin2 = []
expo_test1 = []
expo_test2 = []
expo_test3 = []
expo_test4 = []
expo_fasttest = []

for i in range(len(expo_wav1)):
    expo_sin1.append(expo_wav1[i].copy())
    expo_sin2.append(expo_wav1[i].copy())
    expo_test1.append(expo_wav1[i].copy())
    expo_test2.append(expo_wav1[i].copy())
    expo_test3.append(expo_wav1[i].copy())
    expo_test4.append(expo_wav1[i].copy())
    
for i in range(len(expo_sin1)):
    expo_sin1[i] = np.sin(expo_sin1[i] / 2.0 ) * 100.0 
    expo_sin2[i] = np.sin(expo_sin2[i] / 20.0 ) * 100.0 
    expo_test1[i] = np.sin(expo_test1[i] * np.pi / 150.0 ) * 84.0 * ( np.pi / 2.0 )
    expo_test2[i] = np.random.standard_normal(len(expo_sin1[i])) * 84.0
    expo_test3[i] = np.sin(expo_test3[i] / 150.0 ) * 20.0 
    expo_test4[i] = np.sin(expo_test4[i] / 150.0 ) / 10.0 
        
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
    monte_sin_1.append([])
    monte_sin_2.append([])
    monte_test_1.append([])
    monte_test_2.append([])
    monte_test_3.append([])
    monte_test_4.append([])
    monte_expo_1.append([])
    monte_expo_2.append([])
    monte_expo_3.append([])
    monte_fast_test1.append([])
    monte_fast_test2.append([])
    monte_fast_test3.append([])

monte_sin_1[3] = expo_wav1
monte_sin_1[4] = expo_sin1
monte_sin_2[3] = expo_wav1
monte_sin_2[4] = expo_sin2
monte_test_1[3] = expo_wav1
monte_test_2[3] = expo_wav1
monte_test_3[3] = expo_wav1
monte_test_4[3] = expo_wav1
monte_test_1[4] = expo_test1
monte_test_2[4] = expo_test2
monte_test_3[4] = expo_test3
monte_test_4[4] = expo_test4
monte_expo_1[3] = expo_wav1
monte_expo_1[4] = expo_cal1
monte_expo_2[3] = expo_wav2
monte_expo_2[4] = expo_cal2
monte_expo_3[3] = expo_wav3
monte_expo_3[4] = expo_cal3

monte_fast_test1[3] = np.arange(2900,10000,0.2)
monte_fast_test1[4] = np.sin(monte_fast_test1[3] * np.pi / 150) * 250.0


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
            for i in np.where(redshifted[np.where(redshifted > 5000)] < 6200.)[0]:
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
            for i in np.where(redshifted[np.where(redshifted > 5000)] < 6200.)[0]:
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



slope_test1 = [[],[],[]]
slope_expo1 = [[],[],[]]
slope_expo2 = [[],[],[]]
slope_expo3 = [[],[],[]]
zval_test1 = [[],[],[]]
zval_expo1 = [[],[],[]]
zval_expo2 = [[],[],[]]
zval_expo3 = [[],[],[]]

# slope_test1.append(np.concatenate(monte_fast_test1[5][:]))
# slope_expo1.append(np.concatenate(monte_expo_1[5][:]))
# slope_expo2.append(np.concatenate(monte_expo_2[5][:]))
# slope_expo3.append(np.concatenate(monte_expo_3[5][:]))

addruns = 100
systems = 100


for blah in range(10):
    print blah
    add_monte(monte_test_2,200,1,100)


# First
for blah in range(50):
    print blah
    for i in [0,1,2]:
        print "expo1", datetime.datetime.now()
        clear_monte(monte_expo_1,addruns,i*2 + 1,systems)
        slope_expo1[i].append(np.concatenate(monte_expo_1[5][:]))
        zval_expo1[i].append(np.concatenate(monte_expo_1[2][:]))
        
        
        print "expo2", datetime.datetime.now()
        clear_monte(monte_expo_2,addruns,i*2 + 1,systems)
        slope_expo2[i].append(np.concatenate(monte_expo_2[5][:]))
        zval_expo2[i].append(np.concatenate(monte_expo_2[2][:]))


        print "expo3", datetime.datetime.now()
        clear_monte(monte_expo_3,addruns,i*2 + 1,systems)
        slope_expo3[i].append(np.concatenate(monte_expo_3[5][:]))
        zval_expo3[i].append(np.concatenate(monte_expo_3[2][:]))


        print "test1", datetime.datetime.now()
        clear_test_monte(monte_fast_test1,addruns,i*2 + 1,systems)
        slope_test1[i].append(np.concatenate(monte_fast_test1[5][:]))
        zval_test1[i].append(np.concatenate(monte_fast_test1[2][:]))



# Second
slope_expo1[0] = np.concatenate(slope_expo1[0][:])
slope_expo2[0] = np.concatenate(slope_expo2[0][:])
slope_expo3[0] = np.concatenate(slope_expo3[0][:])
slope_test1[0] = np.concatenate(slope_test1[0][:])
slope_expo1[1] = np.concatenate(slope_expo1[1][:])
slope_expo2[1] = np.concatenate(slope_expo2[1][:])
slope_expo3[1] = np.concatenate(slope_expo3[1][:])
slope_test1[1] = np.concatenate(slope_test1[1][:])
slope_expo1[2] = np.concatenate(slope_expo1[2][:])
slope_expo2[2] = np.concatenate(slope_expo2[2][:])
slope_expo3[2] = np.concatenate(slope_expo3[2][:])
slope_test1[2] = np.concatenate(slope_test1[2][:])

np.savetxt('ex1slope0.ascii',slope_expo1[0])
np.savetxt('ex1slope1.ascii',slope_expo1[1])
np.savetxt('ex1slope2.ascii',slope_expo1[2])  
np.savetxt('ex2slope0.ascii',slope_expo2[0])
np.savetxt('ex2slope1.ascii',slope_expo2[1])
np.savetxt('ex2slope2.ascii',slope_expo2[2])                                   
np.savetxt('ex3slope0.ascii',slope_expo3[0])
np.savetxt('ex3slope1.ascii',slope_expo3[1])
np.savetxt('ex3slope2.ascii',slope_expo3[2])
np.savetxt('test1slope0.ascii',slope_test1[0])
np.savetxt('test1slope1.ascii',slope_test1[1])
np.savetxt('test1slope2.ascii',slope_test1[2])

zval_expo1[0] = np.concatenate(zval_expo1[0][:])
zval_expo2[0] = np.concatenate(zval_expo2[0][:])
zval_expo3[0] = np.concatenate(zval_expo3[0][:])
zval_test1[0] = np.concatenate(zval_test1[0][:])
zval_expo1[1] = np.concatenate(zval_expo1[1][:])
zval_expo2[1] = np.concatenate(zval_expo2[1][:])
zval_expo3[1] = np.concatenate(zval_expo3[1][:])
zval_test1[1] = np.concatenate(zval_test1[1][:])
zval_expo1[2] = np.concatenate(zval_expo1[2][:])
zval_expo2[2] = np.concatenate(zval_expo2[2][:])
zval_expo3[2] = np.concatenate(zval_expo3[2][:])
zval_test1[2] = np.concatenate(zval_test1[2][:])

np.savetxt('ex1zval0.ascii',zval_expo1[0])
np.savetxt('ex1zval1.ascii',zval_expo1[1])
np.savetxt('ex1zval2.ascii',zval_expo1[2])  
np.savetxt('ex2zval0.ascii',zval_expo2[0])
np.savetxt('ex2zval1.ascii',zval_expo2[1])
np.savetxt('ex2zval2.ascii',zval_expo2[2])                                   
np.savetxt('ex3zval0.ascii',zval_expo3[0])
np.savetxt('ex3zval1.ascii',zval_expo3[1])
np.savetxt('ex3zval2.ascii',zval_expo3[2])
np.savetxt('test1zval0.ascii',zval_test1[0])
np.savetxt('test1zval1.ascii',zval_test1[1])
np.savetxt('test1zval2.ascii',zval_test1[2])


expo1_by_n = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
expo1_by_z = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
expo2_by_n = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
expo2_by_z = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
expo3_by_n = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
expo3_by_z = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
test1_by_n = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
test1_by_z = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
for i in range(len(monte_expo_1[0])):
    for j in range(len(monte_expo_1[0][i])):
        expo1_by_n[len(monte_expo_1[0][i][j])].append(monte_expo_1[5][i][j])
        expo1_by_z[len(monte_expo_1[0][i][j])].append(monte_expo_1[2][i][j])
for i in range(len(monte_expo_2[0])):
    for j in range(len(monte_expo_2[0][i])):
        expo2_by_n[len(monte_expo_2[0][i][j])].append(monte_expo_2[5][i][j])
        expo2_by_z[len(monte_expo_2[0][i][j])].append(monte_expo_2[2][i][j])
for i in range(len(monte_expo_3[0])):
    for j in range(len(monte_expo_3[0][i])):
        expo3_by_n[len(monte_expo_3[0][i][j])].append(monte_expo_3[5][i][j])
        expo3_by_z[len(monte_expo_3[0][i][j])].append(monte_expo_3[2][i][j])
for i in range(len(monte_fast_test1[0])):
    for j in range(len(monte_fast_test1[0][i])):
        test1_by_n[len(monte_fast_test1[0][i][j])].append(monte_fast_test1[5][i][j])
        test1_by_z[len(monte_fast_test1[0][i][j])].append(monte_fast_test1[2][i][j])


twox = []
twoy = []
two_slope = []
for i in range(len(monte_fast_test1[0])):
    for j in range(len(monte_fast_test1[0][i])):
        two_indices = random.sample(range(len(monte_fast_test1[0][i][j])), 2)
        twox.append([monte_fast_test1[0][i][j][two_indices[0]],monte_fast_test1[0][i][j][two_indices[1]]])
        twoy.append([monte_fast_test1[1][i][j][two_indices[0]],monte_fast_test1[1][i][j][two_indices[1]]])
        two_slope.append(linregress([monte_fast_test1[0][i][j][two_indices[0]],monte_fast_test1[0][i][j][two_indices[1]]],[monte_fast_test1[1][i][j][two_indices[0]],monte_fast_test1[1][i][j][two_indices[1]]])[0])


threex = []
threey = []
three_slope = []
fourx = []
foury = []
four_slope = []
five_slope = []
six_slope = []
seven_slope = []
eight_slope = []
nine_slope = []
for i in range(len(monte_fast_test1[0])):
    for j in range(len(monte_fast_test1[0][i])):
        if(len(monte_fast_test1[0][i][j]) > 2 ):
            three_indices = random.sample(range(len(monte_fast_test1[0][i][j])), 3)
            x_temp = [monte_fast_test1[0][i][j][three_indices[0]],monte_fast_test1[0][i][j][three_indices[1]],monte_fast_test1[0][i][j][three_indices[2]]]
            y_temp = [monte_fast_test1[1][i][j][three_indices[0]],monte_fast_test1[1][i][j][three_indices[1]],monte_fast_test1[1][i][j][three_indices[2]]]
            threex.append(x_temp)
            threey.append(y_temp)
            three_slope.append(linregress(x_temp,y_temp)[0])
        if(len(monte_fast_test1[0][i][j]) > 3 ):
            four_indices = random.sample(range(len(monte_fast_test1[0][i][j])), 4)
            x_temp = [monte_fast_test1[0][i][j][four_indices[0]],monte_fast_test1[0][i][j][four_indices[1]],monte_fast_test1[0][i][j][four_indices[2]],monte_fast_test1[0][i][j][four_indices[3]]]
            y_temp = [monte_fast_test1[1][i][j][four_indices[0]],monte_fast_test1[1][i][j][four_indices[1]],monte_fast_test1[1][i][j][four_indices[2]],monte_fast_test1[1][i][j][four_indices[3]]]
            fourx.append(x_temp)
            foury.append(y_temp)
            four_slope.append(linregress(x_temp,y_temp)[0])
        if(len(monte_fast_test1[0][i][j]) > 4 ):
            five_indices = random.sample(range(len(monte_fast_test1[0][i][j])), 5)
            x_temp = [monte_fast_test1[0][i][j][five_indices[0]],monte_fast_test1[0][i][j][five_indices[1]],monte_fast_test1[0][i][j][five_indices[2]],monte_fast_test1[0][i][j][five_indices[3]], monte_fast_test1[0][i][j][five_indices[4]]]
            y_temp = [monte_fast_test1[1][i][j][five_indices[0]],monte_fast_test1[1][i][j][five_indices[1]],monte_fast_test1[1][i][j][five_indices[2]],monte_fast_test1[1][i][j][five_indices[3]],monte_fast_test1[1][i][j][five_indices[4]] ]
            five_slope.append(linregress(x_temp,y_temp)[0])
        if(len(monte_fast_test1[0][i][j]) > 5 ):
            six_indices = random.sample(range(len(monte_fast_test1[0][i][j])), 6)
            x_temp = [monte_fast_test1[0][i][j][six_indices[0]],monte_fast_test1[0][i][j][six_indices[1]],monte_fast_test1[0][i][j][six_indices[2]],monte_fast_test1[0][i][j][six_indices[3]], monte_fast_test1[0][i][j][six_indices[4]],monte_fast_test1[0][i][j][six_indices[5]]]
            y_temp = [monte_fast_test1[1][i][j][six_indices[0]],monte_fast_test1[1][i][j][six_indices[1]],monte_fast_test1[1][i][j][six_indices[2]],monte_fast_test1[1][i][j][six_indices[3]],monte_fast_test1[1][i][j][six_indices[4]], monte_fast_test1[1][i][j][six_indices[5]] ]
            six_slope.append(linregress(x_temp,y_temp)[0])
    
# 
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
#
Ntrans = [[],[]]
hmm = 2
for k in range(23):
    Ntrans[0].append(hmm)
    Ntrans[1].append([])
    for i in range(len(monte_test_2[0])):
        for j in range(len(monte_test_2[0][i])):
            if(len(monte_test_2[0][i][j]) == hmm):
                Ntrans[1][k].append(monte_test_2[5][i][j])
    hmm = hmm + 1

for i in range(22):
    print i+2, np.average(Ntrans[1][i]), np.std(Ntrans[1][i])







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

# Third
# testx = []
# testy = []
# for k in range(3):
#     testx.append([])
#     testy.append([])
#     for j in range(20,2000,20):
#         print j, np.std([np.average(slope_test[k][i:i+j]) for i in xrange(0,len(slope_test[k]),j)] )
#         testx[k].append(j)
#         testy[k].append(np.std([np.average(slope_test[k][i:i+j]) for i in xrange(0,len(slope_test[k]),j)] ))
# 
# # Fourth
# expox = []
# expoy = []
# for k in range(3):
#     expox.append([])
#     expoy.append([])
#     for j in range(20,2000,20):
#         print j, np.std([np.average(slope_expo[k][i:i+j]) for i in xrange(0,len(slope_expo[k]),j)] )
#         expox[k].append(j)
#         expoy[k].append(np.std([np.average(slope_expo[k][i:i+j]) for i in xrange(0,len(slope_expo[k]),j)] ))


# pl.scatter(expox[0],expoy[0], color="blue", label="Expo 1")
# pl.scatter(expox[1],expoy[1], color="green", label="Expo 2")
# pl.scatter(expox[2],expoy[2], color="red", label="Expo 3")
# pl.legend()
# 
# pl.scatter(testx[0],testy[0], color="blue", label="small range, 50")
# pl.scatter(testx[1],testy[1], color="green", label="large range, 250")
# pl.scatter(testx[2],testy[2], color="red", label="large range, 250")
# pl.legend()


# np.savetxt("Output/test2_" + runtime,np.average(monte_sin_2[5],1))
# add_monte(monte_expo_1,10,2)
# np.savetxt("Output/raw1_" + runtime,np.average(monte_expo_1[5],1))
    # add_monte(monte_expo_2,100)
    # np.savetxt("Output/raw2_" + runtime,np.average(monte_expo_2[5],1))
    # add_monte(monte_expo_3,100)
    # np.savetxt("Output/raw3_" + runtime,np.average(monte_expo_3[5],1))

# add_monte(monte_sin_1,90,2)
# add_monte(monte_sin_2,90,2)
# 
# for i in range(2,10):
#     clear_monte(monte_sin_1,50,i)
#     print "Test 1: ", i, np.average(monte_sin_1[5]), np.std(monte_sin_1[5])
#     clear_monte(monte_sin_2,50,i)
#     print "Test 2: ", i, np.average(monte_sin_2[5]), np.std(monte_sin_2[5])
# 
# for i in range(2,10):
#     clear_monte(monte_expo_1,50,i)
#     print "Expo 1: ", i, np.average(monte_expo_1[5]), np.std(monte_expo_1[5])
#     clear_monte(monte_expo_2,50,i)
#     print "Expo 2: ", i, np.average(monte_expo_2[5]), np.std(monte_expo_2[5])


# np.savetxt('string',np.average(monte_sin_1[5],1))




# for slope in [monte_sin_1,monte_expo_1, monte_expo_2, monte_expo_3]:
#     for i in range(len(slope[0])):
#         slope[5].append([])
#         for j in range(len(slope[0][0])):
#             slope[5][i].append(linregress(slope[0][i][j],slope[1][i][j])[0])
#     print np.average(slope[5]), np.std(slope[5])    
# 
# totalx_sin_1 = []
# totaly_sin_1 = []
# totalx_expo_1 = []
# totalx_expo_2 = []
# totalx_expo_3 = []
# totaly_expo_1 = []
# totaly_expo_2 = []
# totaly_expo_3 = []
# 
# for i in range(len(monte_sin_1[0][:])):
#     totalx_sin_1.append([])
#     totaly_sin_1.append([])
#     totalx_sin_1[i] = np.concatenate(monte_sin_1[0][i][:])
#     totaly_sin_1[i] = np.concatenate(monte_sin_1[1][i][:])
# 
# for i in range(len(monte_expo_1[0][:])):
#     totalx_expo_1.append([])
#     totaly_expo_1.append([])
#     totalx_expo_1[i] = np.concatenate(monte_expo_1[0][i][:])
#     totaly_expo_1[i] = np.concatenate(monte_expo_1[1][i][:])
# 
# for i in range(len(monte_expo_2[0][:])):
#     totalx_expo_2.append([])
#     totaly_expo_2.append([])
#     totalx_expo_2[i] = np.concatenate(monte_expo_2[0][i][:])
#     totaly_expo_2[i] = np.concatenate(monte_expo_2[1][i][:])
# 
# for i in range(len(monte_expo_3[0][:])):
#     totalx_expo_3.append([])
#     totaly_expo_3.append([])
#     totalx_expo_3[i] = np.concatenate(monte_expo_3[0][i][:])
#     totaly_expo_3[i] = np.concatenate(monte_expo_3[1][i][:])
# 
# # totalx_expo_1 = np.concatenate(totalx_expo_1)
# # totalx_expo_2 = np.concatenate(totalx_expo_2)
# # totalx_expo_3 = np.concatenate(totalx_expo_3)
# # totaly_expo_1 = np.concatenate(totaly_expo_1)
# # totaly_expo_2 = np.concatenate(totaly_expo_2)
# # totaly_expo_3 = np.concatenate(totaly_expo_3)
# 
# # np.savetxt('raw1',np.transpose((totalx_expo_1,totaly_expo_1)))
# # np.savetxt('raw2',np.transpose((totalx_expo_2,totaly_expo_2)))
# # np.savetxt('raw3',np.transpose((totalx_expo_3,totaly_expo_3)))
# 
# # print "expo 1:", linregress(totalx_expo_1,totaly_expo_1)
# # print "expo 2:", linregress(totalx_expo_2,totaly_expo_2)
# # print "expo 3:", linregress(totalx_expo_3,totaly_expo_3)
# 
# print "average along axis"
# 
# # add_monte(monte_sin_1,100)
# # print "progress"
# # add_monte(monte_expo_1,100)
# # print "progress"
# # add_monte(monte_expo_2,100)
# # print "progress"
# # add_monte(monte_expo_3,100)
# # print "Done."
# # monte_sin_1[5] = []
# # # monte_sin_2[5] = []
# # monte_expo_1[5] = []
# monte_expo_2[5] = []
# monte_expo_3[5] = []
# for slope in [monte_sin_1,monte_expo_1, monte_expo_2, monte_expo_3]:
#     for i in range(len(slope[0])):
#         slope[5].append([])
#         for j in range(len(slope[0][0])):
#             slope[5][i].append(linregress(slope[0][i][j],slope[1][i][j])[0])
#     print np.average(slope[5]), np.std(slope[5])    
# 
# totals_expo_1 = []
# totals_expo_2 = []
# totals_expo_3 = []
# totals_expo_1 = np.concatenate(monte_expo_1[5][:])
# totals_expo_2 = np.concatenate(monte_expo_2[5][:])
# totals_expo_3 = np.concatenate(monte_expo_3[5][:])
# 
# for i in range(len(monte_expo_2[0][:])):
#     totals_expo_2.append([])
#     totals_expo_2[i] = np.concatenate(monte_expo_2[5][i][:])
# 
# for i in range(len(monte_expo_3[0][:])):
#     totals_expo_3.append([])
#     totals_expo_3[i] = np.concatenate(monte_expo_3[5][i][:])
#     
# 
# totals_expo_1 = np.concatenate(totals_expo_1)
# print np.average(totals_expo_1)
# totals_expo_2 = np.concatenate(totals_expo_2)
# print np.average(totals_expo_2)
# totals_expo_3 = np.concatenate(totals_expo_3)
# print np.average(totals_expo_3)

# totaly_expo_1 = np.concatenate(totaly_expo_1)
# totaly_expo_2 = np.concatenate(totaly_expo_2)
# totaly_expo_3 = np.concatenate(totaly_expo_3)
# np.savetxt('raw1',np.transpose((totalx_expo_1,totaly_expo_1)))
# np.savetxt('raw2',np.transpose((totalx_expo_2,totaly_expo_2)))
# np.savetxt('raw3',np.transpose((totalx_expo_3,totaly_expo_3)))
# np.savetxt('raw1',totals_expo_1)
# np.savetxt('raw2',totals_expo_2)
# np.savetxt('raw3',totals_expo_3)
# 
# 
# 
# while len(monte_sin_1_z) < runs:
#     delta_alpha = []
#     x_val = []
#     z_record = []
#     # Murphy had 143 systems
#     while len(delta_alpha) < systems:
#         average_vshift = []
#         save_index = []
#         # Take list of lab wavelength values, mult by random redshift
#         z = np.random.uniform(0.2,3.7)
#         # print z
#         redshifted = lwav * (1 + z * 1.0e8)
#         # Only have calibration data in 5000 - 6200 range
#         for i in np.where(redshifted[np.where(redshifted > 5000)] < 6200.)[0]:
#             for j in range(len(expo_wav1)):
#                 if(np.min(np.abs(expo_wav1[j] - redshifted[i])) < tolerance):
#                     average_vshift.append(expo_sin1[j][np.argmin(np.abs(expo_wav1[j] - redshifted[i]))])
#                     save_index.append(i)
#                     break
#         vshift = np.array(average_vshift)
#         if(len(vshift) > 1):
#             z_record.append(z)
#             delta_alpha.append( vshift ) 
#             x_val.append(2.0 * c * qval[save_index] * lwav[save_index])
#     if(len(monte_sin_1_z) % 50 == 0):
#          print "still working... "
#     monte_sin_1_v.append(delta_alpha)
#     monte_sin_1_x.append(x_val)
#     monte_sin_1_z.append(z_record)
# print "Sin 1 Summary: ", np.average(monte_sin_1), np.std(monte_sin_1)
# 
# 
# slope = []
# for i in range(len(monte_sin_1_x)):
#     slope.append([])
#     for j in range(len(monte_sin_1_x[0])):
#         slope[i].append(linregress(monte_sin_1_x[i][j],monte_sin_1_v[i][j])[0])
# 

# while len(monte_sin_1) < runs:
#     delta_alpha = []
#     # Murphy had 143 systems
#     while len(delta_alpha) < systems:
#         average_vshift = []
#         save_index = []
#         # Take list of lab wavelength values, mult by random redshift
#         z = np.random.uniform(0.2,3.7)
#         # print z
#         redshifted = lwav * (1 + z * 1.0e8)
#         # Only have calibration data in 5000 - 6200 range
#         for i in np.where(redshifted[np.where(redshifted > 5000)] < 6200.)[0]:
#             for j in range(len(expo_wav1)):
#                 if(np.min(np.abs(expo_wav1[j] - redshifted[i])) < tolerance):
#                     average_vshift.append(expo_sin1[j][np.argmin(np.abs(expo_wav1[j] - redshifted[i]))])
#                     save_index.append(i)
#                     break
#         vshift = average_vshift - np.average(average_vshift)
#         if(len(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])) > 1):
#             delta_alpha.append(np.average(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])))
#     if(len(monte_sin_1) % 50 == 0):
#         print "check sin: ", len(delta_alpha), len(monte_sin_1), np.average(delta_alpha)
#     monte_sin_1.append(delta_alpha)
# print "Sin 1 Summary: ", np.average(monte_sin_1), np.std(monte_sin_1)



# 
# 
# monte_exposure_1 = []
# while len(monte_exposure_1) < runs:
#     delta_alpha = []
#     # print "safe"
#     # Murphy had 143 systems
#     while len(delta_alpha) < systems:
#         average_vshift = []
#         save_index = []
#         # Take list of lab wavelength values, mult by random redshift
#         z = np.random.uniform(0.2,3.7)
#         # print z
#         redshifted = lwav * (1 + z * 1.0e8)
#         # Only have calibration data in 5000 - 6200 range
#         for i in np.where(redshifted[np.where(redshifted > 5000)] < 6200.)[0]:
#             for j in range(len(expo_wav1)):
#                 if(np.min(np.abs(expo_wav1[j] - redshifted[i])) < tolerance):
#                     average_vshift.append(expo_cal1[j][np.argmin(np.abs(expo_wav1[j] - redshifted[i]))])
#                     save_index.append(i)
#                     break
#         vshift = average_vshift - np.average(average_vshift)
#         if(len(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])) > 1):
#             delta_alpha.append(np.average(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])))
#     if(len(monte_exposure_1) % 50 == 0):
#         print "check exposure 1: ", len(delta_alpha), len(monte_exposure_1), np.average(delta_alpha)
#     monte_exposure_1.append(delta_alpha)
# print "Exposure 1 Summary: ", np.average(monte_exposure_1), np.std(monte_exposure_1)
# 
# monte_exposure_2 = []
# while len(monte_exposure_2) < runs:
#     delta_alpha = []
#     # print "safe"
#     # Murphy had 143 systems
#     while len(delta_alpha) < systems:
#         average_vshift = []
#         save_index = []
#         # Take list of lab wavelength values, mult by random redshift
#         z = np.random.uniform(0.2,3.7)
#         # print z
#         redshifted = lwav * (1 + z * 1.0e8)
#         # Only have calibration data in 5000 - 6200 range
#         for i in np.where(redshifted[np.where(redshifted > 5000)] < 6200.)[0]:
#             for j in range(len(expo_wav2)):
#                 if(np.min(np.abs(expo_wav2[j] - redshifted[i])) < tolerance):
#                     average_vshift.append(expo_cal2[j][np.argmin(np.abs(expo_wav2[j] - redshifted[i]))])
#                     save_index.append(i)
#                     break
#         vshift = average_vshift - np.average(average_vshift)
#         if(len(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])) > 1):
#             delta_alpha.append(np.average(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])))
#     if(len(monte_exposure_2) % 50 == 0):
#         print "check exposure 2: ", len(delta_alpha), len(monte_exposure_2), np.average(delta_alpha)
#     monte_exposure_2.append(delta_alpha)
# print "Exposure 2 Summary: ", np.average(monte_exposure_2), np.std(monte_exposure_2)
# 
# monte_exposure_3 = []
# while len(monte_exposure_3) < runs:
#     delta_alpha = []
#     # print "safe"
#     # Murphy had 143 systems
#     while len(delta_alpha) < systems:
#         average_vshift = []
#         save_index = []
#         # Take list of lab wavelength values, mult by random redshift
#         z = np.random.uniform(0.2,3.7)
#         # print z
#         redshifted = lwav * (1 + z * 1.0e8)
#         # Only have calibration data in 5000 - 6200 range
#         for i in np.where(redshifted[np.where(redshifted > 5000)] < 6200.)[0]:
#             for j in range(len(expo_wav3)):
#                 if(np.min(np.abs(expo_wav3[j] - redshifted[i])) < tolerance):
#                     average_vshift.append(expo_cal3[j][np.argmin(np.abs(expo_wav3[j] - redshifted[i]))])
#                     save_index.append(i)
#                     break
#         vshift = average_vshift - np.average(average_vshift)
#         if(len(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])) > 1):
#             delta_alpha.append(np.average(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])))
#     if(len(monte_exposure_3) % 50 == 0):
#         print "check exposure 3: ", len(delta_alpha), len(monte_exposure_3), np.average(delta_alpha)
#     monte_exposure_3.append(delta_alpha)
# print "Exposure 3 Summary: ", np.average(monte_exposure_3), np.std(monte_exposure_3)
# 
# 


# runs = 5000
# systems = 143
# monte_exposure_1 = []
# while len(monte_exposure_1) < runs:
#     delta_alpha = []
#     # Murphy had 143 systems
#     while len(delta_alpha) < systems:
#         average_vshift = []
#         save_index = []
#         # Take list of lab wavelength values, mult by random redshift
#         z = np.random.uniform(0.2,3.7)
#         # print z
#         redshifted = lwav * (1 + z * 1.0e8)
#         # Only have calibration data in 5000 - 6200 range
#         for i in np.where(redshifted[np.where(redshifted > 5000)] < 6200.)[0]:
#             for j in range(len(expo_wav1)):
#                 if(np.min(np.abs(expo_wav1[j] - redshifted[i])) < tolerance):
#                     average_vshift.append(expo_cal1[j][np.argmin(np.abs(expo_wav1[j] - redshifted[i]))])
#                     save_index.append(i)
#                     break
#         vshift = average_vshift - np.average(average_vshift)
#         # print - vshift / (2.0 * c * qval[save_index] * lwav[save_index])
#         if(len(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])) > 0):
#             delta_alpha.append(np.average(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])))
#     if len(monte_exposure_1) % 25 == 0:
#         print "Expo 1", len(delta_alpha), len(monte_exposure_1), np.average(delta_alpha)
#     monte_exposure_1.append(delta_alpha)
# print "Expo 1 Summary: ", np.average(monte_exposure_1)
# 
# monte_exposure_2 = []
# while len(monte_exposure_2) < runs:
#     delta_alpha = []
#     # Murphy had 143 systems
#     while len(delta_alpha) < systems:
#         average_vshift = []
#         save_index = []
#         # Take list of lab wavelength values, mult by random redshift
#         z = np.random.uniform(0.2,3.7)
#         # print z
#         redshifted = lwav * (1 + z * 1.0e8)
#         # Only have calibration data in 5000 - 6200 range
#         for i in np.where(redshifted[np.where(redshifted > 5000)] < 6200.)[0]:
#             for j in range(len(expo_wav2)):
#                 if(np.min(np.abs(expo_wav2[j] - redshifted[i])) < tolerance):
#                     average_vshift.append(expo_cal2[j][np.argmin(np.abs(expo_wav2[j] - redshifted[i]))])
#                     save_index.append(i)
#                     break
#         vshift = average_vshift - np.average(average_vshift)
#         if(len(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])) > 0):
#             delta_alpha.append(np.average(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])))
#     if len(monte_exposure_2) % 25 == 0:
#         print "Expo 2", len(delta_alpha), len(monte_exposure_2), np.average(delta_alpha)
#     monte_exposure_2.append(delta_alpha)
# print "Expo 2 Summary: ", np.average(monte_exposure_2)
# 
# monte_exposure_3 = []
# while len(monte_exposure_3) < runs:
#     delta_alpha = []
#     # Murphy had 143 systems
#     while len(delta_alpha) < systems:
#         average_vshift = []
#         save_index = []
#         # Take list of lab wavelength values, mult by random redshift
#         z = np.random.uniform(0.2,3.7)
#         # print z
#         redshifted = lwav * (1 + z * 1.0e8)
#         # Only have calibration data in 5000 - 6200 range
#         for i in np.where(redshifted[np.where(redshifted > 5000)] < 6200.)[0]:
#             for j in range(len(expo_wav3)):
#                 if(np.min(np.abs(expo_wav3[j] - redshifted[i])) < tolerance):
#                     average_vshift.append(expo_cal3[j][np.argmin(np.abs(expo_wav3[j] - redshifted[i]))])
#                     save_index.append(i)
#                     break
#         vshift = average_vshift - np.average(average_vshift)
#         if(len(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])) > 0):
#             delta_alpha.append(np.average(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])))
#     if len(monte_exposure_3) % 25 == 0:
#         print "Expo 3", len(delta_alpha), len(monte_exposure_3), np.average(delta_alpha)
#     monte_exposure_3.append(delta_alpha)
# print "Expo 3 Summary: ", np.average(monte_exposure_3)
# 
# 
# print "Expo 1 Summary: ", np.average(monte_exposure_1)
# print "Expo 2 Summary: ", np.average(monte_exposure_2)
# print "Expo 3 Summary: ", np.average(monte_exposure_3)





# 
# carlo2 = []
# for monte in range(runs):
#     z_list = np.random.uniform(0.8,3.0,130)
#     hold_alpha = []
#     for zcount in range(len(z_list)):
#         average_vshift = []
#         save_index = []
#         redshifted = lwav * (1 + z_list[zcount]) * 1.0e8
#         for i in np.where(redshifted[np.where(redshifted > 5033.)] < 6205.)[0]:
#             for j in range(len(expo_wav2)):
#                 if(np.min(np.abs(expo_wav2[j] - redshifted[i])) < tolerance):
#                     average_vshift.append(expo_cal2[j][np.argmin(np.abs(expo_wav2[j] - redshifted[i]))])
#                     save_index.append(i)
#                     break
#         vshift = average_vshift - np.average(average_vshift)
#         hold_alpha.append(np.array(save_index) / (2.0 * c * qval[save_index] * lwav[save_index]))
#     delta_alpha = []
#     for i in range(len(z_list)):
#         if(len(hold_alpha[i]) > 0):
#             delta_alpha.append(np.average(hold_alpha[i]))
#     carlo2.append(np.average(delta_alpha))
#     print "Expo 2", len(delta_alpha), np.average(delta_alpha), monte
# 
# carlo3 = []
# for monte in range(runs):
#     z_list = np.random.uniform(0.8,3.0,130)
#     hold_alpha = []
#     for zcount in range(len(z_list)):
#         average_vshift = []
#         save_index = []
#         redshifted = lwav * (1 + z_list[zcount]) * 1.0e8
#         for i in np.where(redshifted[np.where(redshifted > 5000)] < 6200.)[0]:
#             for j in range(len(expo_wav3)):
#                 if(np.min(np.abs(expo_wav3[j] - redshifted[i])) < tolerance):
#                     average_vshift.append(expo_cal3[j][np.argmin(np.abs(expo_wav3[j] - redshifted[i]))])
#                     save_index.append(i)
#                     break
#         vshift = average_vshift - np.average(average_vshift)
#         hold_alpha.append(np.array(save_index) / (2.0 * c * qval[save_index] * lwav[save_index]))
#     delta_alpha = []
#     for i in range(len(z_list)):
#         if(len(hold_alpha[i]) > 0):
#             delta_alpha.append(np.average(hold_alpha[i]))
#     carlo3.append(np.average(delta_alpha))
#     print "Expo 3", len(delta_alpha), np.average(delta_alpha), monte