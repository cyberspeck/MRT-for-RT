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

"""

import FunITK as fun
from FunITK import Volume
import datetime
import numpy as np
# https://docs.scipy.org/doc/numpy/reference/generated/numpy.set_printoptions.html
# np.set_printoptions(threshold=np.nan)

# CT images not usable @173-192
# MR images has airbubble @306

pathCT = "../data_final/CT_x1/"
pathMR = "../data_final/MR_x1/"
pathCT_x4 = "../data_final/CT_x4/"
pathMR_x4 = "../data_final/MR_x4/"
pathCT_x9 = "../data_final/CT_x9/"
pathMR_x9 = "../data_final/MR_x9/"
#pathCT_x16 = "../data_final/CT_x16/"
#pathMR_x16 = "../data_final/MR_x16/"
pathCT_x25 = "../data_final/CT_x25/"
pathMR_x25 = "../data_final/MR_x25/"
pathCT_x100 = "../data_final/CT_x100/"
pathMR_x100 = "../data_final/MR_x100/"
idxSlice = 10

CT = Volume(path=pathCT, method="CT", resample=1, ref=idxSlice)
CT_x4 = Volume(path=pathCT_x4, method="CT", resample=4, ref=idxSlice)
CT_x9 = Volume(path=pathCT_x9, method="CT", resample=9, ref=idxSlice)
#CT_x16 = Volume(path=pathCT_x16, method="CT", resample=16, ref=idxSlice)
CT_x25 = Volume(path=pathCT_x25, method="CT", resample=25, ref=idxSlice)
CT_x100 = Volume(path=pathCT_x100, method="CT", resample=100, ref=idxSlice)

MR = Volume(path=pathMR, method="MR", resample=1, ref=idxSlice)
MR_x4 = Volume(path=pathMR_x4, method="MR", resample=4, ref=idxSlice)
MR_x9 = Volume(path=pathMR_x9, method="MR", resample=9, ref=idxSlice)
#MR_x16 = Volume(path=pathMR_x16, method="MR", resample=16, ref=idxSlice)
MR_x25 = Volume(path=pathMR_x25, method="MR", resample=25, ref=idxSlice)
MR_x100 = Volume(path=pathMR_x100, method="MR", resample=100, ref=idxSlice)

vol_list = [[CT, CT_x4, CT_x9, CT_x25, CT_x100],[MR, MR_x4, MR_x9, MR_x25, MR_x100]]
modality, sets = np.shape(vol_list)

sliceNumbers = np.arange(CT.zSize, dtype=int)
warp = np.zeros((sets, CT.zSize, 2))
warpMagnitude = np.zeros((sets, CT.zSize, 1))

# one calculated with DC and one with DC opti
warpDC = np.zeros((sets, CT.zSize, 1))
warpDC_opti = np.zeros((sets, CT.zSize, 1))
    
# 2 DC for CT, 4 DC for MR (2 using MR.centroid, 2 using CT.centroid!)
DC_CT = np.zeros((sets, CT.zSize, 2))
DC_CT_average = np.zeros((sets, 2))
DC_MR = np.zeros((sets, CT.zSize, 4))
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
    warpDC[i] = warpMagnitude[i] * (1 - a)
    warpDC_opti[i] = warpMagnitude[i] * (1 - b)


'''
fig = plt.figure()
plt.ylim(ymin=-2.1, ymax=.5)
plt.xlim(xmin=0, xmax=CT.zSize)
plt.plot(warp[i])
plt.ylabel(u"warp [mm]")
plt.xlabel(u"slice")

i=4
fig = plt.figure()
plt.ylim(ymin=.3, ymax=1)
plt.xlim(xmin=0, xmax=CT.zSize)
plt.plot(DC_CT[i,:,1])
plt.plot(DC_MR[i,:,1])
plt.plot(DC_MR[i,:,3])
plt.ylabel(u"DC")

fig = plt.figure()
plt.ylim(ymin=0, ymax=0.21)
plt.xlim(xmin=0, xmax=CT.zSize)
plt.plot(warpDC[i])
plt.plot(warpDC_opti[i])
plt.ylabel(u"warp*DC [mm]")
plt.xlabel(u"slice")
'''

# http://stackoverflow.com/questions/16621351/how-to-use-python-numpy-savetxt-to-write-strings-and-float-number-to-an-ascii-fi
now = datetime.datetime.now()
COLUMNS  = 'sliceNr  warp_X  warp_Y  warpMagnitude  DC_CT  DC_CT_opti  DC_MR  DC_MR_opti  DC_MR_CT-COM  DC_MR_opti_CT-COM  warpDC  warpDC_opti'
for i in range(sets):
    DATA = np.column_stack((sliceNumbers.astype(str),
                            warp[i].round(4).astype(str),
                            warpMagnitude[i].round(4).astype(str),
                            DC_CT[i,:,0].round(4).astype(str),
                            DC_CT[i,:,1].round(4).astype(str),
                            DC_MR[i,:,0].round(4).astype(str),
                            DC_MR[i,:,1].round(4).astype(str),
                            DC_MR[i,:,2].round(4).astype(str),
                            DC_MR[i,:,3].round(4).astype(str),
                            warpDC[i,:].round(4).astype(str),
                            warpDC_opti[i,:].round(4).astype(str)))
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
    np.savetxt('CT-MR_x{}_{}_{}.txt'.format(vol_list[0][i].resample, 
               now.date(), now.time()), DATA, delimiter="     ", header=head,
               comments="# ", fmt='%3s')

'''
# creates mask (pixel values either 0 or 1)
for i in range(sets):
    vol_list[0][i].getMask()
# creates CT.masked using CT.mask,
# but assigns each slice the centroid distance*1000*spacing as pixel value
    vol_list[0][i].applyMask(replaceArray=warpMagnitude[i])
# exports 3D image as .mha file
    fun.sitk_write(vol_list[0][i].masked, "../data_final/export/", "{}_x{}_warpMagnitude.mha".format(vol_list[0][i].method, vol_list[0][i].resample))
    
    vol_list[0][i].applyMask(replaceArray=DC_MR[i,:,1])
    fun.sitk_write(vol_list[0][i].masked, "../data_final/export/", "{}_x{}_DC-MR-opti.mha".format(vol_list[0][i].method, vol_list[0][i].resample))
    
    vol_list[0][i].applyMask(replaceArray=DC_MR[i,:,3])
    fun.sitk_write(vol_list[0][i].masked, "../data_final/export/", "{}_x{}_DC-MR-opti_CT-COM.mha".format(vol_list[0][i].method, vol_list[0][i].resample))
    
    vol_list[0][i].applyMask(replaceArray=warpDC_opti[i], scale=10000)
    fun.sitk_write(vol_list[0][i].masked, "../data_final/export/", "{}_x{}_DC-warpDC_opti.mha".format(vol_list[0][i].method, vol_list[0][i].resample))


# instead of opening the created file manually, you can use this lines in
# the IPython console to start 3D Slicer and open it automatically:
# %env SITK_SHOW_COMMAND /home/davidblacher/Downloads/Slicer-4.5.0-1-linux-amd64/Slicer
# %env SITK_SHOW_COMMAND /home/david/Downloads/Slicer-4.5.0-1-linux-amd64/Slicer
# sitk.Show(CT.masked)
'''



