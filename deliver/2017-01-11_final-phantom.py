
# coding: utf-8

# In[1]:

import FunctionsCustom as fun
from FunctionsCustom import Volume
import datetime
import numpy as np


# In[2]:

pathCT = "../data_final/CT_x1/"
pathMR = "../data_final/MR_x1/"
pathCT_x4 = "../data_final/CT_x4/"
pathMR_x4 = "../data_final/MR_x4/"
pathCT_x9 = "../data_final/CT_x9/"
pathMR_x9 = "../data_final/MR_x9/"
pathCT_x16 = "../data_final/CT_x16/"
pathMR_x16 = "../data_final/MR_x16/"
pathCT_x25 = "../data_final/CT_x25/"
pathMR_x25 = "../data_final/MR_x25/"
pathCT_x100 = "../data_final/CT_x100/"
pathMR_x100 = "../data_final/MR_x100/"
idxSlice = 10


# In[3]:

CT = Volume(path=pathCT, method="CT", info="x1", ref=idxSlice)
CT_x4 = Volume(path=pathCT_x4, method="CT", info="x4", ref=idxSlice)
CT_x9 = Volume(path=pathCT_x9, method="CT", info="x9", ref=idxSlice)
CT_x16 = Volume(path=pathCT_x16, method="CT", info="x16", ref=idxSlice)
CT_x25 = Volume(path=pathCT_x25, method="CT", info="x25", ref=idxSlice)
CT_x100 = Volume(path=pathCT_x100, method="CT", info="x100", ref=idxSlice)

MR = Volume(path=pathMR, method="MR", info="x1", ref=idxSlice)
MR_x4 = Volume(path=pathMR_x4, method="MR", info="x4", ref=idxSlice)
MR_x9 = Volume(path=pathMR_x9, method="MR", info="x9", ref=idxSlice)
MR_x16 = Volume(path=pathMR_x16, method="MR", info="x16", ref=idxSlice)
MR_x25 = Volume(path=pathMR_x25, method="MR", info="x25", ref=idxSlice)
MR_x100 = Volume(path=pathMR_x100, method="MR", info="x100", ref=idxSlice)


# In[6]:

vol_list = [[CT, CT_x4, CT_x9, CT_x16, CT_x25, CT_x100],[MR, MR_x4, MR_x9, MR_x16, MR_x25, MR_x100]]

for volumes in vol_list:
    for vol in volumes:
        #vol.getCentroid(percentLimit='auto', iterations=10)
        vol.getCentroid(threshold='auto')
        vol.showCentroid()
        vol.getMask()
        vol.getDice()


# In[8]:

sliceNumbers = np.arange(CT.zSize, dtype=int)
(methods, sets) = np.shape(vol_list)
distortion = np.zeros((sets, CT.zSize, 2))
distortionNorm = np.zeros((sets, CT.zSize, 1))
dice_CT_MR = np.zeros((sets, CT.zSize, 2))
for index in range(sets):
    # this calculates the coordinate difference of MR.centroid relative to CT.centroid
    distortion[index] = fun.coordShift(vol_list[0][index].centroid, vol_list[1][index].centroid)
    # this calculates the norm (=absolute distance) between the centroids in each slice
    distortionNorm[index] = fun.coordDist(distortion[index])
    # collects dice-score of CT and MRI data
    dice_CT_MR[index] = np.column_stack((vol_list[0][index].dice, vol_list[1][index].dice))


# In[9]:

now = datetime.datetime.now()
NAMES  = ['sliceNumber', 'distortionX', 'distortionY', 'distortionNorm', 'dice_CT', 'dice_MR']
for index in range(sets):
    DATA = np.column_stack((sliceNumbers.astype(str), distortion[index].astype(str), distortionNorm[index].astype(str), dice_CT_MR[index,:,0].astype(str), dice_CT_MR[index,:,1].astype(str)))
    text = np.row_stack((NAMES, DATA))
    head0 = "{}_{}\n path: {}\n thresholds: {}, {}\n".format(vol_list[0][index].method, vol_list[0][index].info, vol_list[0][index].path, vol_list[0][index].lower, vol_list[0][index].upper)
    head1 = "{}_{}\n path: {}\n thresholds: {}, {}\n".format(vol_list[1][index].method, vol_list[1][index].info, vol_list[1][index].path, vol_list[1][index].lower, vol_list[1][index].upper)
    head = str(now) + '\n'+ head0 + head1
    np.savetxt('CT-MR_{}_{}_{}.txt'.format(vol_list[0][index].info, now.date(), now.time()), text, delimiter="   ", header=head, comments="# ", fmt='%3s')

