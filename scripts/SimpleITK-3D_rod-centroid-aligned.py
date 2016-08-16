# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 15:17:57 2016

For future developer: Segmentation and automatic registration of separate rods
could be done similar to this:https://blancosilva.wordpress.com/2010/12/15/image-processing-with-numpy-scipy-and-matplotlibs-in-sage/

@author: david
"""

import FunITK as fun
from FunITK import Volume

pathCT = "../data/cropped_CT-a/"
pathMR = "../data/cropped_MR-d-a/"
idxSlice = 10

CT = Volume(path=pathCT, method="CT", ref=idxSlice, seeds=[(6, 8, idxSlice)])
CT.getCentroid()
CT.showCentroid()

MR = Volume(path=pathMR, method="MR", ref=idxSlice, seeds=[(6, 8, idxSlice)])
MR.getCentroid()
MR.showCentroid()

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

# instead of opening the created file manually, you can use this lines in
# the IPython console to start 3D Slicer and open it automatically:
# %env SITK_SHOW_COMMAND /home/david/Downloads/Slicer-4.5.0-1-linux-amd64/Slicer
# sitk.Show(CT.masked)