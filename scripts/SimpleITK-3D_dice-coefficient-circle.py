# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 13:20:39 2016

@author: david
"""


import FunITK as fun
from FunITK import Volume

pathCT = "../data/cropped_CT-a/"
pathMR = "../data/cropped_MR-d-a/"
idxSlice = 10

CT = Volume(path=pathCT, method="CT", ref=idxSlice, seeds=[(6, 8, idxSlice)])
CT.getCentroid()
CT.getMask()

MR = Volume(path=pathMR, method="MR", ref=idxSlice, seeds=[(6, 8, idxSlice)])
MR.getCentroid(new=True)
MR.showCentroid()
MR.getMask(new=True)
MR.showMask()

rCT=CT.getDice()
rMR=MR.getDice()


CT.applyMask(replaceArray=CT.dice)
# fun.sitk_write(CT.masked)

# print("Dice coefficient for each slice of MR-Volume (radius = 2.3): \n{}".format(
# fun.dice(MR.mask, MR.centroid, radius=2.3)))


# instead of opening the created file manually, you can use this lines in
# the IPython console to start 3D Slicer and open it automatically:
# %env SITK_SHOW_COMMAND /home/david/Downloads/Slicer-4.5.0-1-linux-amd64/Slicer
# sitk.Show(CT.mask)