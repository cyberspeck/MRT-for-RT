# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 15:17:57 2016

For future developer: Segmentation and automatic registration of separate rods
could be done similar to this:https://blancosilva.wordpress.com/2010/12/15/image-processing-with-numpy-scipy-and-matplotlibs-in-sage/

@author: david
"""

import FunITK as fun
from FunITK import Volume

#CT images not usable 176-192
pathCT = "../data_final/CT_x1/"
pathMR = "../data_final/MR_x1/"
pathCT_x4 = "../data_final/CT_x4/"
pathMR_x4 = "../data_final/MR_x4/"
pathCT_x9 = "../data_final/CT_x9/"
pathMR_x9 = "../data_final/MR_x9/"
idxSlice = 10

CT = Volume(path=pathCT, method="CT", info="x1", ref=idxSlice, seeds=[(4, 20, idxSlice)])
#a = CT.getCentroid(percentLimit='auto', iterations=10)
#CT.showCentroid()
b = CT.getCentroid(threshold='auto')
CT.showCentroid()
#CT.getMask()
#CT.getDice(b, CT.mask)

CT_x4 = Volume(path=pathCT_x4, method="CT", info="x4", ref=idxSlice, seeds=[(8, 40, idxSlice)])
#c = CT_x4.getCentroid(percentLimit='auto', iterations=10)
#CT_x4.showCentroid()
d = CT_x4.getCentroid(threshold='auto')
CT_x4.showCentroid()
#CT_x4.getMask()
#CT_x4.getDice(d, CT_x4.mask)

CT_x9 = Volume(path=pathCT_x9, method="CT", info="x9", ref=idxSlice, seeds=[(12, 60, idxSlice)])
#e = CT_x9.getCentroid(percentLimit='auto', iterations=5)
#CT_x9.showCentroid()
f = CT_x9.getCentroid(threshold='auto')
CT_x9.showCentroid()
#CT_x9.getMask()
#CT_x9.getDice(f, CT_x9.mask)


MR = Volume(path=pathMR, method="MR", info="x1", ref=idxSlice, seeds=[(7, 21, idxSlice)])
#MR.showSeed(pixel=True)
#g = MR.getCentroid(percentLimit='auto')
#MR.showCentroid()
h = MR.getCentroid(threshold='auto')
MR.showCentroid()
#MR.getMask()
#MR.getDice(h, MR.mask)

MR_x4 = Volume(path=pathMR_x4, method="MR", info="x4", ref=idxSlice, seeds=[(14, 42, idxSlice)])
#MR_x4.showSeed(pixel=True)
#i = MR_x4.getCentroid(percentLimit='auto', iterations=5)
#MR_x4.showCentroid()
j = MR_x4.getCentroid(threshold='auto')
MR_x4.showCentroid()
#MR_x4.getMask()
#MR_x4.getDice(i, MR_x4.mask)

MR_x9 = Volume(path=pathMR_x9, method="MR", info="x9", ref=idxSlice, seeds=[(21, 63, idxSlice)])
#k = MR_x9.getCentroid(percentLimit='auto', iterations=10, halfShift=0.1)
#MR_x9.showCentroid()
l = MR_x9.getCentroid(threshold='auto')
MR_x9.showCentroid()
#MR_x9.getMask()
#MR_x9.getDice(l, MR_x9.mask)
'''
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