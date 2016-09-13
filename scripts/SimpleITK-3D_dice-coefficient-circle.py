# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 13:20:39 2016

@author: david
"""


import FunITK as fun
from FunITK import Volume

pathCT = "../data/cropped_CT/"
pathMR = "../data/cropped_MR/"
pathCT_04 = "../data/cropped_CT_resample_04/"
pathMR_04 = "../data/cropped_MR_resample_04/"
pathCT_012 = "../data/cropped_CT_resample_012/"
pathMR_012 = "../data/cropped_MR_resample_012/"
idxSlice = 10


CT = Volume(path=pathCT, method="CT", ref=idxSlice, seeds=[(8, 6, idxSlice)])
#CT.getThresholds()
#CT.getMask()
#CT.applyMask()
#CT.getCentroid(threshold='auto', scale = 1)
#CT.showCentroid()
'''
limits = np.linspace(0.1, 1.9, num=20)
centroids = np.zeros((len(limits), CT.zSize, 2))
cts = np.zeros((len(limits), 2))
for index, p in enumerate(limits, start=0):
    centroids[index] = CT.getCentroid(threshold='auto', scale = p)
    print(p)
    cts[index] = centroids[index,idxSlice]
# plot differences in centroids compared to centroid calculated with
# threshold='auto' using coordDist()!
'''

'''
MR = Volume(path=pathMR, method="MR", ref = idxSlice, seeds=[(8, 6, idxSlice)])
MR.getCentroid()
MR.showCentroid()
MR.getMask()
MR.showMask()


CT_012 = Volume(path=pathCT_012, method="CT", ref=idxSlice, seeds=[(80, 60, idxSlice)])
#CT_012.show()
#CT_012.getThresholds()
CT_012.getCentroid()
CT_012.showCentroid()
CT_012.getMask()
CT_012.showMask()
'''
#CT_04 = Volume(path=pathCT_04, method="CT", ref=idxSlice, seeds=[(25, 20, idxSlice)])
#CT_04.getThreshold()
#CT_04.showCentroid()
#CT_04.showMask()

#MR_04 = Volume(path=pathMR_04, method="MR", ref = idxSlice, seeds=[(25, 20, idxSlice)])
#MR_04.getThreshold()


'''
MR = Volume(path=pathMR, method="MR", ref=idxSlice, seeds=[(6, 8, idxSlice)])
MR.getCentroid()
MR.showCentroid()
MR.getMask()
MR.showMask()

rMR=MR.getDice(show=idxSlice)
'''

# CT.applyMask(replaceArray=CT.dice)
# fun.sitk_write(CT.masked)

# print("Dice coefficient for each slice of MR-Volume (radius = 2.3): \n{}".format(fun.dice_circle(MR.mask, MR.centroid, radius=rMR)))


# instead of opening the created file manually, you can use this lines in
# the IPython console to start 3D Slicer and open it automatically:
# %env SITK_SHOW_COMMAND /home/david/Downloads/Slicer-4.5.0-1-linux-amd64/Slicer
# sitk.Show(CT.mask)