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

ct = Volume(path=pathCT, method="CT", ref=idxSlice)
ct.centroid(show=idxSlice)

mr = Volume(path=pathMR, method="MR", ref=idxSlice)
mr.centroid(show=idxSlice)


distortion = fun.coordShift(ct.centroid, mr.centroid)
print(distortion)
