#!/usr/bin/env python
# encoding: utf-8
"""
load.py

Created by Jonathan Whitmore on 2009-11-16
Updated: 2009-12-14
"""

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
# # print np.average(rancal1 - cal1) / (2 * c * (q1611*l16110 - q1608*l16080) ), \
# # np.std(rancal1 - cal1) / (2 * c * (q1611*l16110 - q1608*l16080) )
# # print np.average(rancal2 - cal2) / (2 * c * (q1611*l16110 - q1608*l16080) ), \
# # np.std(rancal2 - cal2) / (2 * c * (q1611*l16110 - q1608*l16080) )
# # print np.average(rancal3 - cal3) / (2 * c * (q1611*l16110 - q1608*l16080) ), \
# # np.std(rancal3 - cal3) / (2 * c * (q1611*l16110 - q1608*l16080) )
# # 
# # print np.average(random1 / t1608_1611), np.std(random1 / t1608_1611)
# # print np.average(random2 / t1608_1611), np.std(random2 / t1608_1611)
# # print np.average(random3 / t1608_1611), np.std(random3 / t1608_1611)
# # 
# # print np.average(random1 / t1741_1751), np.std(random1 / t1741_1751)
# # print np.average(random2 / t1741_1751), np.std(random2 / t1741_1751)
# # print np.average(random3 / t1741_1751), np.std(random3 / t1741_1751)
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
# 
# 
#         # 138: pl.hist(all_cal,range=[-500,850],bins=200, color="black", alpha=.25)
#         # 139: pl.hist(cal1,range=[-500,850],bins=200, color="blue", alpha=.5,label="Exposure 1")
#         # 140: pl.hist(cal2,range=[-500,850],bins=200, color="green", alpha=.5,label="Exposure 2")
#         # 141: pl.hist(cal3,range=[-500,850],bins=200, color="red", alpha=.5,label="Exposure 3")
#         # 142: pl.legend()
#         # 143: pl.xlabel("velocity offset (m/s)")
#         # 144: pl.ylabel("counts")
#         # 145: np.average(cal1)
#         # 146: pl.average(cal1)
#         # 147: print pl.average(cal1), pl.std(cal1)
#         # 148: print pl.average(cal2), pl.std(cal2)
#         # 149: print pl.average(cal3), pl.std(cal3)
#         # 150: pl.average(all_cal)
#         # 151: pl.average(orcal1)
#         # 152: pl.average(orcal2)
#         # 153: pl.average(orcal3)
#         # 155: pl.title("Histogram of Offsets by Exposure")
#         # 156: pl.savefig("histogram.eps")
#         # 157: pl.savefig("histogram1.pdf")
# 
# 
# pl.hist(np.average(monte_total,1),range=[-0.000025,0.000015],bins=50,color="black",alpha=.4,label="All")
# pl.hist(np.average(monte_exposure_1,1),range=[-0.000025,0.000015],bins=50,color="blue",alpha=.4,label="Exposure 1")
# pl.hist(np.average(monte_exposure_2,1),range=[-0.000025,0.000015],bins=50,color="green",alpha=.4,label="Exposure 2")
# pl.hist(np.average(monte_exposure_3,1),range=[-0.000025,0.000015],bins=50,color="red",alpha=.4,label="Exposure 3")
# pl.vlines(np.average(monte_total),0,450,color="blue",label='Average')
# pl.vlines(0,0,1200,color="black",label='Zero')
# pl.title('Monte Carlo 143 systems')
# pl.legend()
# pl.savefig('monte-2000.png')
# 
# pl.hist(np.average(monte_sin_1,1),range=[-0.000025,0.000015],bins=50,color="blue",alpha=.4,label="Sin Offset")
# pl.vlines(np.average(monte_sin_1),0,200,color="blue",label='Average')
# pl.vlines(0,0,200,color="black",label='Zero')
# pl.title('Monte Carlo Sin Test')
# pl.legend()
# pl.savefig('monte-sin-test.png')
# 
# pl.hist(np.average(steph,1),range=[-0.000025,0.000015],bins=50,color="black",alpha=.2)
# pl.hist(np.average(monte_exposure_1,1),range=[-0.000025,0.000015],bins=50,color="blue",alpha=.2)
# pl.hist(np.average(monte_exposure_2,1),range=[-0.000025,0.000015],bins=50,color="green",alpha=.2)
# pl.hist(np.average(monte_exposure_3,1),range=[-0.000025,0.000015],bins=50,color="red",alpha=.2)
# pl.hist(np.average(monte_exposure_1,1),range=[-0.000025,0.000015],bins=50,color="blue",alpha=.2)
# pl.hist(np.average(monte_exposure_2,1),range=[-0.000025,0.000015],bins=50,color="green",alpha=.2)
# pl.hist(np.average(monte_exposure_3,1),range=[-0.000025,0.000015],bins=50,color="red",alpha=.2)
# pl.vlines(0,0,1200,color="black")
# pl.vlines(np.average(steph),0,1200,color="blue")
# 
# pl.hist(np.average(monte_sin_1,1),bins=25)
# pl.vlines(0,0,700)
# pl.title('Sin and Offset')

pl.hist(cal1,bins=100,range=[-500,1000],color="blue",alpha=0.5,label="Expo 1")
pl.hist(cal2,bins=100,range=[-500,1000],color="green",alpha=0.5,label="Expo 2")
pl.hist(cal3,bins=100,range=[-500,1000],color="red",alpha=0.5,label="Expo 3")
pl.xlabel("Miscalibration in m/s")
pl.ylabel("Counts")
pl.title("Histogram of Offsets by Exposure")
pl.legend()
pl.savefig('Histogram1.eps')
pl.close()
