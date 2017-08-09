#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 15:17:57 2016
@author: david


plt.clf()

plt.xlim((0,1)
plt.yscale('log')

plt.legend(loc='center right')
plt.legend(loc='lower right')
plt.legend(loc='upper right')
plt.legend(loc='lower left')

plt.xlabel(u"guess [%]")  
plt.ylabel(u"DC")

plt.ylabel(u"warp [mm]")
plt.xlabel(u"slice")

plt.tight_layout()

# print entire array:
np.set_printoptions(threshold='nan')

"""

import FunITK as fun
from FunITK import Volume
import datetime
import numpy as np
import matplotlib.pyplot as plt
import os

idxSlice = 200
ph2_CT = Volume(path="../data/phantom2/CT_x1", method="CT", resample=1, ref=idxSlice)
ph2_CT_x4 = Volume(path="../data/phantom2/CT_x4", method="CT", resample=4, ref=idxSlice)
ph2_CT_x9 = Volume(path="../data/phantom2/CT_x9", method="CT", resample=9, ref=idxSlice)
ph2_CT_x25 = Volume(path="../data/phantom2/CT_x25", method="CT", resample=25, ref=idxSlice)
ph2_CT_x100 = Volume(path="../data/phantom2/CT_x100", method="CT", resample=100, ref=idxSlice)

ph2_MR = Volume(path="../data/phantom2/MR_x1", method="MR", resample=1, ref=idxSlice)
ph2_MR_x4 = Volume(path="../data/phantom2/MR_x4", method="MR", resample=4, ref=idxSlice)
ph2_MR_x9 = Volume(path="../data/phantom2/MR_x9", method="MR", resample=9, ref=idxSlice)
ph2_MR_x25 = Volume(path="../data/phantom2/MR_x25", method="MR", resample=25, ref=idxSlice)
ph2_MR_x100 = Volume(path="../data/phantom2/MR_x100", method="MR", resample=100, ref=idxSlice)

vol_list = [[ph2_CT, ph2_CT_x4, ph2_CT_x9, ph2_CT_x25, ph2_CT_x100],[ph2_MR, ph2_MR_x4, ph2_MR_x9, ph2_MR_x25, ph2_MR_x100]]
modality, sets = np.shape(vol_list)

length = ph2_CT.zSize
spacing = ph2_CT.zSpace
sliceNumbers = np.arange(length, dtype=int)

# for data centered around iso-centre, this is real x-axis:
iso = 183
dist = ( (sliceNumbers - iso ) ).round(2)

warp = np.zeros((sets, length, 2))
warpMagnitude = np.zeros((sets, length, 1))
    
# 2 DC for CT, 4 DC for MR (2 using MR.centroid, 2 using CT.centroid!)
DC_CT = np.zeros((sets, length, 2))
DC_CT_average = np.zeros((sets, 2))
DC_MR = np.zeros((sets, length, 4))
DC_MR_average = np.zeros((sets, 4))
iterate = 51

for i in range(sets):
    vol_list[0][i].getCentroid()
    vol_list[1][i].getCentroid()
    # this calculates the coordinate difference of MR.centroid relative to CT.centroid
    warp[i] = fun.sitk_coordShift(vol_list[0][i].centroid, vol_list[1][i].centroid)
    # this calculates the norm (=absolute distance) between the centroids in each slice
    warpMagnitude[i] = fun.sitk_coordDist(warp[i])
    
    

#fig = plt.figure()
#plt.ylim(ymin=0.8, ymax=1)
#plt.xlim(xmin=(3.5-.1), xmax=(4.5+.1))
#calculates DC for CT
for i in range(sets):
    vol_list[0][i].getCentroid()
    a = vol_list[0][i].getDice()
    aa = vol_list[0][i].diceAverage
    b = vol_list[0][i].getDice(centroid=vol_list[0][i].centroid, iterations=iterate,
                               # plot=True,
                               # save='{}_x{}-{}iter'.format(vol_list[1][i].method, vol_list[1][i].resample, iterate)
                               )
    bb = vol_list[0][i].diceAverage
    DC_CT[i] = np.column_stack((a,b))
    DC_CT_average[i] = aa,bb
    

#fig = plt.figure()
#plt.ylim(ymin=0.5, ymax=.95)
#plt.xlim(xmin=(1.8-.1), xmax=(2.8+.1))
#calculates DC for MR, using first CT COM, then its own COM
#this way self.bestRadius is still set to the radius yielding the best DC
#independently of the CT COM
for i in range(sets):
    vol_list[1][i].getCentroid()
    c = vol_list[1][i].getDice(centroid=vol_list[0][i].centroid)
    cc = vol_list[1][i].diceAverage
    d = vol_list[1][i].getDice(centroid=vol_list[0][i].centroid, iterations=iterate)
    dd = vol_list[1][i].diceAverage
    a = vol_list[1][i].getDice()
    aa = vol_list[1][i].diceAverage
    b = vol_list[1][i].getDice(centroid=vol_list[1][i].centroid, iterations=iterate,
                               # plot=True,
                               # save='{}_x{}-{}iter'.format(vol_list[1][i].method, vol_list[1][i].resample, iterate)
                               )
    bb = vol_list[1][i].diceAverage

    DC_MR[i] = np.column_stack((a,b,c,d))
    DC_MR_average[i] = aa,bb,cc,dd

'''
# x and y warp
fig = plt.figure()
#plt.ylim(ymin=-2.1, ymax=.5)
plt.xlim(xmin=dist[0], xmax=dist[-1])
plt.plot(dist, warp[i])
plt.legend(('x-shift', 'y-shift'),loc=2)
plt.ylabel(u"warp [mm]")
plt.xlabel(u"z-axis [mm]")
#plt.title('Economic Cost over Time')
#plt.show()


# warpMagnitude
fig = plt.figure()
plt.ylim(ymin=0, ymax=warpMagnitude[i].max())
plt.xlim(xmin=dist[0], xmax=dist[-1])
plt.plot(dist, warpMagnitude[i])
#plt.legend(('warpMagnitude'),loc=0)
plt.ylabel(u"warp [mm]")
plt.xlabel(u"z-axis [mm]")

# DC for CT and MRI and MRI (CT COM)
fig = plt.figure()
plt.ylim(ymin=0, ymax=1)
plt.xlim(xmin=dist[0], xmax=dist[-1])
plt.plot(dist, DC_CT[i,:,1])
plt.plot(dist, DC_MR[i,:,1])
plt.plot(dist, DC_MR[i,:,3])
plt.legend(('CT', 'MR', 'MR (CT COM)'),loc=0)
plt.ylabel(u"DC")
plt.xlabel(u"z-axis [mm]")

# warpDC
fig = plt.figure()
plt.ylim(ymin=0, ymax=warpDC_opti[i].max())
plt.xlim(xmin=dist[0], xmax=dist[-1])
plt.plot(dist, warpDC[i])
plt.plot(dist, warpDC_opti[i])
plt.legend(('DC', 'DC optimised'),loc=0)
plt.ylabel(u"warpDC [mm]")
plt.xlabel(u"z-axis [mm]")
'''


'''
# http://stackoverflow.com/questions/16621351/how-to-use-python-numpy-savetxt-to-write-strings-and-float-number-to-an-ascii-fi
now = datetime.datetime.now()

COLUMNS  = 'sliceNo dist warp_x  warp_y  warpMagnitude  DC_CT  DC_CT_opti  DC_MR  DC_MR_opti  DC_MR_CT-COM  DC_MR_opti_CT-COM'
#COLUMNS  = 'sliceNo & dist & $warp_x$  & $warp_y$  & $warp$ & $DC_{CT}$  & $DC^*_{CT}$ & $DC_{MR}$  & $DC^*_{MR}$ & $DC_{MR(CT-COM)}$ & $DC^*_{MR(CT-COM)}$'
for i in range(sets):
    DATA = np.column_stack((sliceNumbers.astype(str),
                            dist.astype(str),
                            warp[i].round(4).astype(str),
                            warpMagnitude[i].round(4).astype(str),
                            DC_CT[i,:,0].round(4).astype(str),
                            DC_CT[i,:,1].round(4).astype(str),
                            DC_MR[i,:,0].round(4).astype(str),
                            DC_MR[i,:,1].round(4).astype(str),
                            DC_MR[i,:,2].round(4).astype(str),
                            DC_MR[i,:,3].round(4).astype(str)))
 #   text = np.row_stack((NAMES, DATA))
    head0 = "{}_x{}\n path: {}\n thresholds: {}, {}\n DC-average: {} (radius = 4)\n DC-average (opti): {} (bestRadius: {})\n".format(vol_list[0][i].method,
    vol_list[0][i].resample, vol_list[0][i].path, vol_list[0][i].lower,
    vol_list[0][i].upper, DC_CT_average[i][0].round(4),
    DC_CT_average[i][1].round(4), vol_list[0][i].bestRadius)
    
    head1 = "{}_x{}\n path: {}\n thresholds: {}, {}\n DC-average: {} (radius = 2)\n DC-average (opti): {} (bestRadius: {})\n DC-average (CT-COM): {}\n DC-average (CT-COM, opti): {}\n".format(vol_list[1][i].method,
    vol_list[1][i].resample, vol_list[1][i].path, vol_list[1][i].lower,
    vol_list[1][i].upper, DC_MR_average[i][0].round(4),
    DC_MR_average[i][1].round(4), vol_list[1][i].bestRadius,
    DC_MR_average[i][2].round(4), DC_MR_average[i][3].round(4))
    
    head = str(now) + '\n'+ head0 + head1 + '\n' + COLUMNS
    np.savetxt('../data/output_txt/phantom2_out_txt/CT-MR_x{}_{}_{}.txt'.format(vol_list[0][i].resample, 
               now.date(), now.time()), DATA, delimiter="   &  ", header=head,
               comments="# ", fmt='%3s')
'''


'''
# creates mask (pixel values either 0 or 1)
for i in range(sets):
    vol_list[0][i].getMask()
# creates CT.masked using CT.mask,
# but assigns each slice the centroid distance*1000*spacing as pixel value
    vol_list[0][i].applyMask(replaceArray=warpMagnitude[i])
# exports 3D image as .mha file
    fun.sitk_write(vol_list[0][i].masked, "../data_final_export_-1/", "{}_x{}_warpMagnitude.mha".format(vol_list[0][i].method, vol_list[0][i].resample))
    
    vol_list[0][i].applyMask(replaceArray=DC_MR[i,:,1])
    fun.sitk_write(vol_list[0][i].masked, "../data_final_export_-1/", "{}_x{}_DC-MR-opti.mha".format(vol_list[0][i].method, vol_list[0][i].resample))
    
    vol_list[0][i].applyMask(replaceArray=DC_MR[i,:,3])
    fun.sitk_write(vol_list[0][i].masked, "../data_final_export_-1/", "{}_x{}_DC-MR-opti_CT-COM.mha".format(vol_list[0][i].method, vol_list[0][i].resample))


# instead of opening the created file manually, you can use this lines in
# the IPython console to start 3D Slicer and open it automatically:
# %env SITK_SHOW_COMMAND /home/davidblacher/Downloads/Slicer-4.5.0-1-linux-amd64/Slicer
# %env SITK_SHOW_COMMAND /home/david/Downloads/Slicer-4.5.0-1-linux-amd64/Slicer
# sitk.Show(CT.masked)
'''