# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 15:17:57 2016

For future developer: Segmentation and automatic registration of separate rods
could be done similar to this: https://blancosilva.wordpress.com/2010/12/15/image-processing-with-numpy-scipy-and-matplotlibs-in-sage/
@author: david

Separate Matplotlib figure windows pop up for Spyder 2.3/3.0:
Tools > Preferences > Ipython Console > Graphics > Graphics Backend > Backend: “automatic”


plt.clf()

plt.xlim((0,1)
plt.yscale('log')
plt.legend(loc='center right')
plt.xlabel(u"guess [%]")  
plt.ylabel(u"dc")
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
pathCT_x16 = "../data_final/CT_x16/"
pathMR_x16 = "../data_final/MR_x16/"
pathCT_x25 = "../data_final/CT_x25/"
pathMR_x25 = "../data_final/MR_x25/"
pathCT_x100 = "../data_final/CT_x100/"
pathMR_x100 = "../data_final/MR_x100/"
idxSlice = 10

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

vol_list = [[CT, CT_x4, CT_x9, CT_x16, CT_x25, CT_x100],[MR, MR_x4, MR_x9, MR_x16, MR_x25, MR_x100]]

for volumes in vol_list:
    for vol in volumes:
        #vol.getCentroid(percentLimit='auto', iterations=10)
        vol.getCentroid(threshold='auto')
        vol.showCentroid()
        vol.getMask()
        vol.getDice()

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


# http://stackoverflow.com/questions/16621351/how-to-use-python-numpy-savetxt-to-write-strings-and-float-number-to-an-ascii-fi
now = datetime.datetime.now()
NAMES  = ['sliceNumber', 'distortionX', 'distortionY', 'distortionNorm', 'dice_CT', 'dice_MR']
for index in range(sets):
    DATA = np.column_stack((sliceNumbers.astype(str), distortion[index].astype(str), distortionNorm[index].astype(str), dice_CT_MR[index,:,0].astype(str), dice_CT_MR[index,:,1].astype(str)))
    text = np.row_stack((NAMES, DATA))
    head0 = "{}_{}\n path: {}\n thresholds: {}, {}\n".format(vol_list[0][index].method, vol_list[0][index].info, vol_list[0][index].path, vol_list[0][index].lower, vol_list[0][index].upper)
    head1 = "{}_{}\n path: {}\n thresholds: {}, {}\n".format(vol_list[1][index].method, vol_list[1][index].info, vol_list[1][index].path, vol_list[1][index].lower, vol_list[1][index].upper)
    head = str(now) + '\n'+ head0 + head1
    np.savetxt('CT-MR_{}_{}_{}.txt'.format(vol_list[0][index].info, now.date(), now.time()), text, delimiter="   ", header=head, comments="# ", fmt='%3s')


'''
print("\n")
print("CT, centroid[:,:,{}]: {}".format(idxSlice, CT.centroid[idxSlice]))
print("MR, centroid[:,:,{}]: {}".format(idxSlice, MR.centroid[idxSlice]))
print("distrotion[:,:,{}]: {}".format(idxSlice, distortion[idxSlice]))


print("distrotionNorm[:,:,{}]: {}".format(idxSlice, distortionNorm[idxSlice]))

# creates mask (pixel values either 0 or 1)
CT.getMask()
# creates CT.masked using CT.mask,
# but assigns each slice the centroid distance*1000*spacing as pixel value
CT.applyMask(replaceArray=distortionNorm, spacing=CT.img.GetSpacing()[0])

# exports 3D image as .mha file
fun.sitk_write(CT.masked, "../data/", "CT_distortionNorm.mha")
'''
# instead of opening the created file manually, you can use this lines in
# the IPython console to start 3D Slicer and open it automatically:
# %env SITK_SHOW_COMMAND /home/davidblacher/Downloads/Slicer-4.5.0-1-linux-amd64/Slicer
# %env SITK_SHOW_COMMAND /home/david/Downloads/Slicer-4.5.0-1-linux-amd64/Slicer
# sitk.Show(CT.masked)
