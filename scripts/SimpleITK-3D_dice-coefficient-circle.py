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


CT = Volume(path=pathCT, denoise=True, method="CT", ref=idxSlice, seeds=[(8, 6, idxSlice)])
CT.show()
'''
CT.getCentroid()
CT.showCentroid()
CT.getMask()
CT.showMask()
rCT=CT.getDice()

reader = sitk.ImageSeriesReader()
filenames = reader.GetGDCMSeriesFileNames(pathCT_012)
reader.SetFileNames(filenames)
img = reader.Execute()
fun.sitk_show(img)
centroid=fun.sitk_centroid(img)
fun.centroid_show(img, centroid)
blub=fun.dice_circle(img, radius=30, centroid=centroid, show=1)


CT_4 = Volume(path=pathCT_04)
CT_4.show(ref=10)
fun.sitk_show(CT_4.img, ref=10)

CT_04.getCentroid()
CT_04.showCentroid()
CT_04.getMask()
CT_04.showMask()
rCT_04=CT_04.getDice()

CT_012 = Volume(path=pathCT_012, method="CT", ref=idxSlice, seeds=[(6, 8, idxSlice)])
CT_012.getCentroid()
CT_012.showCentroid()
CT_012.getMask()
CT_012.showMask()
rCT_012=CT_012.getDice()


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