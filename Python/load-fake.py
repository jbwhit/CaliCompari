#!/usr/bin/env python
# encoding: utf-8
"""
load.py

Created by Jonathan Whitmore on 2010-02-09
"""

import csv
import numpy as np
import scipy as sp
from scipy import arange, optimize, special, interpolate
import scipy.interpolate as si
import scipy.signal as ss 
from scipy.ndimage import *
import minuit as mi # remove from ubuntu
import matplotlib.pylab as pl # remove from ubuntu  
import sys
import os
import glob
import random as ra
from scipy.stats import linregress


c = 299792458.0

bin_file_list = glob.glob('0*.ascii')
bin_file_list.sort()

wav = []
cal = []
err = []
orwav = []
orcal = []
orerr = []

for filename in bin_file_list:
    
    tempwav, tempcal, temperr, temporwav, temporcal, temporerr = np.loadtxt(filename,unpack=True)
    wav.append(tempwav)
    cal.append(tempcal) 
    err.append(temperr)
    orwav.append(temporwav)
    orcal.append(temporcal)
    orerr.append(temporerr)

# Proper one 
print bin_file_list[:23]
print bin_file_list[23:46]
print bin_file_list[46:]

expo_wav1 = wav[:23]
expo_wav2 = wav[23:46]
expo_wav3 = wav[46:]
expo_cal1 = cal[:23]
expo_cal2 = cal[23:46]
expo_cal3 = cal[46:]

wav1 = np.concatenate(wav[:23])
cal1 = np.concatenate(cal[:23])
err1 = np.concatenate(err[:23])
orwav1 = np.concatenate(orwav[:23])
orcal1 = np.concatenate(orcal[:23])
orerr1 = np.concatenate(orerr[:23])

wav2 = np.concatenate(wav[23:46])
cal2 = np.concatenate(cal[23:46])
err2 = np.concatenate(err[23:46])
orwav2 = np.concatenate(orwav[23:46])
orcal2 = np.concatenate(orcal[23:46])
orerr2 = np.concatenate(orerr[23:46])

wav3 = np.concatenate(wav[46:])
cal3 = np.concatenate(cal[46:])
err3 = np.concatenate(err[46:])
orwav3 = np.concatenate(orwav[46:])
orcal3 = np.concatenate(orcal[46:])
orerr3 = np.concatenate(orerr[46:])

all_wav = np.concatenate([wav1,wav2,wav3])
all_cal = np.concatenate([cal1,cal2,cal3])
all_err = np.concatenate([err1,err2,err3])
all_orwav = np.concatenate([orwav1,orwav2,orwav3])
all_orcal = np.concatenate([orcal1,orcal2,orcal3])
all_orerr = np.concatenate([orerr1,orerr2,orerr3])


lmid1 = []
lmid2 = []
lmid3 = []
lmid = []

for i in range(len(bin_file_list)/3):
    i1 = i
    i2 = i + 23
    i3 = i + 46
    lmid1.append(cal[i1][len(cal[i1])/2] - orcal[i1][len(cal[i1])/2])
    lmid.append(cal[i1][len(cal[i1])/2] - orcal[i1][len(cal[i1])/2])
    lmid2.append(cal[i2][len(cal[i2])/2] - orcal[i2][len(cal[i2])/2])
    lmid.append(cal[i2][len(cal[i2])/2] - orcal[i2][len(cal[i2])/2])
    lmid3.append(cal[i3][len(cal[i3])/2] - orcal[i3][len(cal[i3])/2])
    lmid.append(cal[i3][len(cal[i3])/2] - orcal[i3][len(cal[i3])/2])

diff = []
diff1 = []
diff2 = []
diff3 = []
for j in np.arange(0,9):
    diff.append([])
    diff1.append([])
    diff2.append([])
    diff3.append([])
    for i in range(len(bin_file_list)/3):
        i1 = i
        i2 = i + 23
        i3 = i + 46
        diff1[j].append(cal[i1][(j+1) * len(cal[i1])/10] - orcal[i1][(j+1) * len(cal[i1])/10])
        diff[j].append(cal[i1][(j+1) * len(cal[i1])/10] - orcal[i1][(j+1) * len(cal[i1])/10])
        diff2[j].append(cal[i2][(j+1) * len(cal[i2])/10] - orcal[i2][(j+1) * len(cal[i2])/10])
        diff[j].append(cal[i2][(j+1) * len(cal[i2])/10] - orcal[i2][(j+1) * len(cal[i2])/10])
        diff3[j].append(cal[i3][(j+1) * len(cal[i3])/10] - orcal[i3][(j+1) * len(cal[i3])/10])
        diff[j].append(cal[i3][(j+1) * len(cal[i3])/10] - orcal[i3][(j+1) * len(cal[i3])/10])
    

# 
# rancal1 = cal1.copy()
# rancal2 = cal2.copy()
# rancal3 = cal3.copy()
# ra.shuffle(rancal1)
# ra.shuffle(rancal2)
# ra.shuffle(rancal3)
# 
# 
# q1608 = -1030.; q1611 = 1560.
# l16080 = 1608.45085e-8; l16110 = 1611.20034e-8
# w16080 = 1./ (l16080); w16110 = 1. / l16110
# 
# q1741 = -1400.; q1751 = -700.
# l17410 = 1741.5531e-8; l17510 = 1751.9157e-8
# w17410 = 1. / l17410; w17510 = 1. / l17510
# 
# t1608_1611 = (2 * c * (q1611*l16110 - q1608*l16080) )
# t1741_1751 = (2 * c * (q1751*l17510 - q1741*l17410) )
# # 
# random1 = (rancal1 - cal1)
# random2 = (rancal2 - cal2)
# random3 = (rancal3 - cal3)
# 
# 
# q_values_list = 'qvalues.txt'
# q_reader = csv.reader(open(q_values_list), delimiter=' ')
# trans = []
# lwav = []
# fval = []
# gamma = []
# mass = []
# qval = []
# 
# for row in q_reader:
#     trans.append(row[0])
#     lwav.append(float(row[1]))
#     fval.append(float(row[2]))
#     gamma.append(float(row[3]))
#     mass.append(float(row[4]))
#     qval.append(float(row[5]))
# 
# qval = np.array(qval)    
# lwav = np.array(lwav) * 1e-8
# 
# # redshifted = lwav * (1 + z) * 1.0e8
# tolerance = 0.2
# expo_sin1 = []
# for i in range(len(expo_wav1)):
#     expo_sin1.append(expo_wav1[i].copy())
# 
# for i in range(len(expo_sin1)):
#     expo_sin1[i] = np.sin(expo_sin1[i]) * 100. + 120.3
# # print "safe"
# runs = 5000
# systems = 143
# monte_sin_1 = []
# while len(monte_sin_1) < runs:
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
#         # print "safe"
#         # Only have calibration data in 5000 - 6200 range
#         for i in np.where(redshifted[np.where(redshifted > 5000)] < 6200.)[0]:
#             for j in range(len(expo_wav1)):
#                 if(np.min(np.abs(expo_wav1[j] - redshifted[i])) < tolerance):
#                     average_vshift.append(expo_sin1[j][np.argmin(np.abs(expo_wav1[j] - redshifted[i]))])
#                     save_index.append(i)
#                     break
#         vshift = average_vshift - np.average(average_vshift)
#         # print "safe"
#         # print - vshift / (2.0 * c * qval[save_index] * lwav[save_index])
#         if(len(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])) > 0):
#             delta_alpha.append(np.average(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])))
#     if(len(monte_sin_1) % 50 == 0):
#         print "Sin-check 1", len(delta_alpha), len(monte_sin_1), np.average(delta_alpha)
#     monte_sin_1.append(delta_alpha)
# print "Sin 1 Summary: ", np.average(monte_sin_1)
# 
# 
# 
# # runs = 5000
# # systems = 143
# # monte_exposure_1 = []
# # while len(monte_exposure_1) < runs:
# #     delta_alpha = []
# #     # Murphy had 143 systems
# #     while len(delta_alpha) < systems:
# #         average_vshift = []
# #         save_index = []
# #         # Take list of lab wavelength values, mult by random redshift
# #         z = np.random.uniform(0.2,3.7)
# #         # print z
# #         redshifted = lwav * (1 + z * 1.0e8)
# #         # Only have calibration data in 5000 - 6200 range
# #         for i in np.where(redshifted[np.where(redshifted > 5000)] < 6200.)[0]:
# #             for j in range(len(expo_wav1)):
# #                 if(np.min(np.abs(expo_wav1[j] - redshifted[i])) < tolerance):
# #                     average_vshift.append(expo_cal1[j][np.argmin(np.abs(expo_wav1[j] - redshifted[i]))])
# #                     save_index.append(i)
# #                     break
# #         vshift = average_vshift - np.average(average_vshift)
# #         # print - vshift / (2.0 * c * qval[save_index] * lwav[save_index])
# #         if(len(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])) > 0):
# #             delta_alpha.append(np.average(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])))
# #     if len(monte_exposure_1) % 25 == 0:
# #         print "Expo 1", len(delta_alpha), len(monte_exposure_1), np.average(delta_alpha)
# #     monte_exposure_1.append(delta_alpha)
# # print "Expo 1 Summary: ", np.average(monte_exposure_1)
# # 
# # monte_exposure_2 = []
# # while len(monte_exposure_2) < runs:
# #     delta_alpha = []
# #     # Murphy had 143 systems
# #     while len(delta_alpha) < systems:
# #         average_vshift = []
# #         save_index = []
# #         # Take list of lab wavelength values, mult by random redshift
# #         z = np.random.uniform(0.2,3.7)
# #         # print z
# #         redshifted = lwav * (1 + z * 1.0e8)
# #         # Only have calibration data in 5000 - 6200 range
# #         for i in np.where(redshifted[np.where(redshifted > 5000)] < 6200.)[0]:
# #             for j in range(len(expo_wav2)):
# #                 if(np.min(np.abs(expo_wav2[j] - redshifted[i])) < tolerance):
# #                     average_vshift.append(expo_cal2[j][np.argmin(np.abs(expo_wav2[j] - redshifted[i]))])
# #                     save_index.append(i)
# #                     break
# #         vshift = average_vshift - np.average(average_vshift)
# #         if(len(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])) > 0):
# #             delta_alpha.append(np.average(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])))
# #     if len(monte_exposure_2) % 25 == 0:
# #         print "Expo 2", len(delta_alpha), len(monte_exposure_2), np.average(delta_alpha)
# #     monte_exposure_2.append(delta_alpha)
# # print "Expo 2 Summary: ", np.average(monte_exposure_2)
# # 
# # monte_exposure_3 = []
# # while len(monte_exposure_3) < runs:
# #     delta_alpha = []
# #     # Murphy had 143 systems
# #     while len(delta_alpha) < systems:
# #         average_vshift = []
# #         save_index = []
# #         # Take list of lab wavelength values, mult by random redshift
# #         z = np.random.uniform(0.2,3.7)
# #         # print z
# #         redshifted = lwav * (1 + z * 1.0e8)
# #         # Only have calibration data in 5000 - 6200 range
# #         for i in np.where(redshifted[np.where(redshifted > 5000)] < 6200.)[0]:
# #             for j in range(len(expo_wav3)):
# #                 if(np.min(np.abs(expo_wav3[j] - redshifted[i])) < tolerance):
# #                     average_vshift.append(expo_cal3[j][np.argmin(np.abs(expo_wav3[j] - redshifted[i]))])
# #                     save_index.append(i)
# #                     break
# #         vshift = average_vshift - np.average(average_vshift)
# #         if(len(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])) > 0):
# #             delta_alpha.append(np.average(- vshift / (2.0 * c * qval[save_index] * lwav[save_index])))
# #     if len(monte_exposure_3) % 25 == 0:
# #         print "Expo 3", len(delta_alpha), len(monte_exposure_3), np.average(delta_alpha)
# #     monte_exposure_3.append(delta_alpha)
# # print "Expo 3 Summary: ", np.average(monte_exposure_3)
# # 
# # 
# # print "Expo 1 Summary: ", np.average(monte_exposure_1)
# # print "Expo 2 Summary: ", np.average(monte_exposure_2)
# # print "Expo 3 Summary: ", np.average(monte_exposure_3)
# 
# 
# 
# 
# 
# # 
# # carlo2 = []
# # for monte in range(runs):
# #     z_list = np.random.uniform(0.8,3.0,130)
# #     hold_alpha = []
# #     for zcount in range(len(z_list)):
# #         average_vshift = []
# #         save_index = []
# #         redshifted = lwav * (1 + z_list[zcount]) * 1.0e8
# #         for i in np.where(redshifted[np.where(redshifted > 5033.)] < 6205.)[0]:
# #             for j in range(len(expo_wav2)):
# #                 if(np.min(np.abs(expo_wav2[j] - redshifted[i])) < tolerance):
# #                     average_vshift.append(expo_cal2[j][np.argmin(np.abs(expo_wav2[j] - redshifted[i]))])
# #                     save_index.append(i)
# #                     break
# #         vshift = average_vshift - np.average(average_vshift)
# #         hold_alpha.append(np.array(save_index) / (2.0 * c * qval[save_index] * lwav[save_index]))
# #     delta_alpha = []
# #     for i in range(len(z_list)):
# #         if(len(hold_alpha[i]) > 0):
# #             delta_alpha.append(np.average(hold_alpha[i]))
# #     carlo2.append(np.average(delta_alpha))
# #     print "Expo 2", len(delta_alpha), np.average(delta_alpha), monte
# # 
# # carlo3 = []
# # for monte in range(runs):
# #     z_list = np.random.uniform(0.8,3.0,130)
# #     hold_alpha = []
# #     for zcount in range(len(z_list)):
# #         average_vshift = []
# #         save_index = []
# #         redshifted = lwav * (1 + z_list[zcount]) * 1.0e8
# #         for i in np.where(redshifted[np.where(redshifted > 5000)] < 6200.)[0]:
# #             for j in range(len(expo_wav3)):
# #                 if(np.min(np.abs(expo_wav3[j] - redshifted[i])) < tolerance):
# #                     average_vshift.append(expo_cal3[j][np.argmin(np.abs(expo_wav3[j] - redshifted[i]))])
# #                     save_index.append(i)
# #                     break
# #         vshift = average_vshift - np.average(average_vshift)
# #         hold_alpha.append(np.array(save_index) / (2.0 * c * qval[save_index] * lwav[save_index]))
# #     delta_alpha = []
# #     for i in range(len(z_list)):
# #         if(len(hold_alpha[i]) > 0):
# #             delta_alpha.append(np.average(hold_alpha[i]))
# #     carlo3.append(np.average(delta_alpha))
# #     print "Expo 3", len(delta_alpha), np.average(delta_alpha), monte