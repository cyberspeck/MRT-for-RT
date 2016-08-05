# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 15:17:57 2016

@author: david
"""

import FunITK as fun
from FunITK import Volume

pathCT = "../data/cropped_CT-a/"
pathMR = "../data/cropped_MR-d-a/"
idxSlice = 10

CT = Volume(path=pathCT, method="CT", ref=idxSlice, seeds=[(6,8,idxSlice)])
# CT.show()
# CT.showSeed()
CT.getMask()
#CT.applyMask()
#CT.showMask()
#CT.showMasked()
CT.getCentroid()
#CT.showCentroid()

MR = Volume(path=pathMR, method="MR", ref=idxSlice, seeds=[(6,8,idxSlice)])
MR.getCentroid()
#MR.showCentroid()

distortion = coordShift(CT.centroid, MR.centroid)
print("CT, centroid[:,:,10]: ", CT.centroid[10])
print("MR, centroid[:,:,10]: ", MR.centroid[10])
print("distrotion[:,:,10]: ", distortion[10])

distortionNorm = np.zeros((CT.zSize, 1))

for slice in range(CT.zSize):
    distortionNorm[slice,:] = np.linalg.norm(distortion[slice,:])

print("distrotionNorm[:,:,10]: ", distortionNorm[10])

CT.applyMask(replaceArray=distortionNorm)

#CT.mask[6,7,CT.zSize-1]
#CT.showMask()
#CT.showMasked()

#%env SITK_SHOW_COMMAND /home/david/Downloads/Slicer-4.5.0-1-linux-amd64/Slicer
#sitk.Show(MR.img)