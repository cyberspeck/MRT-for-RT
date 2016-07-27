# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 15:17:57 2016

@author: david
based on:
https://pyscience.wordpress.com/2014/10/19/image-segmentation-with-python-and-sitk/
http://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/03_Image_Details.html
"""

import FunITK as fun
import matplotlib.pyplot as plt

pathCT = "../data/cropped_CT-a/"
pathMR = "../data/cropped_MR-d-a/"

imgOriginalCT = fun.sitk_read(pathCT)
imgOriginalMR = fun.sitk_read(pathMR)

xSpace, ySpace, zSpace = imgOriginalMR.GetSpacing()
idxSlice = 10

fun.sitk_show(imgOriginalCT[:,:,idxSlice], title="CT, original")
centroidCT = fun.sitk_centroid(imgOriginalCT, show = 1)

fun.sitk_show(imgOriginalMR[:,:,idxSlice], title="MR, original")
centroidMR = fun.sitk_centroid(imgOriginalMR, show = 1)