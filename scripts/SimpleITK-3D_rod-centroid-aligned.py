# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 15:17:57 2016

For future developer: Segmentation and automatic registration of separate rods
could be done similar to this:https://blancosilva.wordpress.com/2010/12/15/image-processing-with-numpy-scipy-and-matplotlibs-in-sage/

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
'''
CT = Volume(path=pathCT, method="CT", ref=idxSlice, seeds=[(7, 8, idxSlice)])
a = CT.getCentroid(percentLimit='auto', iterations=10)
CT.showCentroid()
b = CT.getCentroid(threshold='auto')
CT.getMask()
CT.getDice(b, CT.mask)
CT.showCentroid()

CT_04 = Volume(path=pathCT_04, method="CT", ref=idxSlice, seeds=[(20, 24, idxSlice)])
c = CT_04.getCentroid(percentLimit='auto', iterations=10)
CT_04.showCentroid()
d = CT_04.getCentroid(threshold='auto')
CT_04.getMask()
CT_04.getDice(d, CT_04.mask)
CT_04.showCentroid()

CT_012 = Volume(path=pathCT_012, method="CT", ref=idxSlice, seeds=[(65, 80, idxSlice)])
e = CT_012.getCentroid(percentLimit='auto', iterations=5)
CT_012.showCentroid()
f = CT_012.getCentroid(threshold='auto')
CT_012.showCentroid()
CT_012.getMask()
CT_012.getDice(f, CT_012.mask)

'''
MR = Volume(path=pathMR, method="MR", ref=idxSlice, seeds=[(8, 6, idxSlice)])
MR.showSeed(pixel=True)
g = MR.getCentroid(percentLimit='auto')
MR.showCentroid()
h = MR.getCentroid(threshold='auto')
MR.getMask()
MR.getDice(h, MR.mask)
'''
MR_04 = Volume(path=pathMR_04, method="MR", ref=idxSlice, seeds=[(25, 18, idxSlice)])
MR_04.showSeed(pixel=True)
i = MR_04.getCentroid(percentLimit='auto', iterations=5)
MR_04.showCentroid()
#j = MR_04.getCentroid(threshold='auto')
#MR_04.getMask()
#MR_04.getDice(i, MR_04.mask)

MR_012 = Volume(path=pathMR_012, method="MR", ref=idxSlice, seeds=[(83, 60, idxSlice)])
k = MR_012.getCentroid(percentLimit='auto', iterations=10, halfShift=0.1)
MR_012.showCentroid()
l = MR_012.getCentroid(threshold='auto')
MR_012.getMask()
MR_012.showMask()
MR_012.getDice(l, MR_012.mask)

# this calculates the coordinate difference of MR.centroid relative to CT.centroid
distortion = fun.coordShift(CT.centroid, MR.centroid)
print("\n")
print("CT, centroid[:,:,{}]: {}".format(idxSlice, CT.centroid[idxSlice]))
print("MR, centroid[:,:,{}]: {}".format(idxSlice, MR.centroid[idxSlice]))
print("distrotion[:,:,{}]: {}".format(idxSlice, distortion[idxSlice]))

# this calculates the norm (=absolute distance) between the centroids in each slice
distortionNorm = fun.coordDist(distortion)
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
# %env SITK_SHOW_COMMAND /home/david/Downloads/Slicer-4.5.0-1-linux-amd64/Slicer
# sitk.Show(CT.masked)