
# coding: utf-8

# In[1]:

import os
import numpy as np
from scipy import ndimage
import SimpleITK as sitk
import matplotlib.pyplot as plt
get_ipython().magic('pylab inline')

import FunctionsCustom as fun


# In[2]:

pathCT = "../data/cropped_CT/"
pathMR = "../data/cropped_MR-d/"

imgOriginalCT = fun.sitk_read(pathCT)
imgOriginalMR = fun.sitk_read(pathMR)

idxSlice = 20

labelRod = 1

# imgOriginal_slice = imgOriginal[:,:,idxSlice]

fun.sitk_show(imgOriginalCT[:,:,idxSlice], title="CT, original")
fun.sitk_show(imgOriginalMR[:,:,idxSlice], title="MR, original")


# In[4]:

imgSmoothCT = sitk.CurvatureFlow(image1=imgOriginalCT,
                                    timeStep=0.125,
                                    numberOfIterations=5)

imgSmoothMR = sitk.CurvatureFlow(image1=imgOriginalMR,
                                    timeStep=0.125,
                                    numberOfIterations=5)
                                    
fun.sitk_show(imgSmoothCT[:,:,idxSlice], title="CT, denoised")
fun.sitk_show(imgSmoothMR[:,:,idxSlice], title="MR, denoised")

