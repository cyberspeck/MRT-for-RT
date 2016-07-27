# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 10:22:52 2016

@author: david
based on:
https://pyscience.wordpress.com/2014/10/19/image-segmentation-with-python-and-sitk/
http://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/03_Image_Details.html
"""

import os
import numpy as np
from scipy import ndimage
import SimpleITK as sitk
import matplotlib.pyplot as plt

import FunITK as fun

pathCT = "../data/cropped_CT/"
pathMR = "../data/cropped_MR-d/"

imgOriginalCT = fun.sitk_read(pathCT)
imgOriginalMR = fun.sitk_read(pathMR)

idxSlice = 20

labelPlastic = 2
labelFilling = 1

# imgOriginal_slice = imgOriginal[:,:,idxSlice]

fun.sitk_show(imgOriginalCT[:,:,idxSlice], title="CT, original")
fun.sitk_show(imgOriginalMR[:,:,idxSlice], title="MR, original")

imgSmoothCT = sitk.CurvatureFlow(image1=imgOriginalCT,
                                    timeStep=0.125,
                                    numberOfIterations=5)

imgSmoothMR = sitk.CurvatureFlow(image1=imgOriginalMR,
                                    timeStep=0.125,
                                    numberOfIterations=5)
                                    
fun.sitk_show(imgSmoothCT[:,:,idxSlice], title="CT, denoised")
fun.sitk_show(imgSmoothMR[:,:,idxSlice], title="MR, denoised")



seedFillingCT = [(14,14,idxSlice)]

# https://itk.org/SimpleITKDoxygen/html/classitk_1_1simple_1_1ConnectedThresholdImageFilter.html
'''
imgFillingCT = sitk.ConnectedThreshold(image1=imgSmoothCT, 
                                              seedList=seedFillingCT, 
                                              lower=00, 
                                              upper=110,
                                              replaceValue=labelFilling)
'''

# https://itk.org/SimpleITKDoxygen/html/classitk_1_1simple_1_1NeighborhoodConnectedImageFilter.html
maskFillingCT = sitk.ConfidenceConnected(image1=imgSmoothCT, 
                                              seedList=seedFillingCT,
                                              multiplier=2.5,
                                              numberOfIterations=10,
                                              initialNeighborhoodRadius=2,
                                              #lower=00, 
                                              #upper=110,
                                              replaceValue=labelFilling)

fun.sitk_show(imgSmoothCT[:,:,idxSlice], title="CT, denoised")
fun.sitk_show(maskFillingCT[:,:,idxSlice], title="CT, denoised, filling mask")

imgFillingCT = fun.sitk_mask(imgSmoothCT, maskFillingCT)

fun.sitk_show(imgFillingCT[:,:,idxSlice], title="CT, denoised, filling")

# for shift in range(-10,10):
#    sitk_show(imgFillingCT[:,:,idxSlice+shift])

# to view in 3D Slicer, type this in IPython console or in jupyter notebook:
# %env SITK_SHOW_COMMAND /home/david/Downloads/Slicer-4.5.0-1-linux-amd64/Slicer
# sitk.Show(imgFillingCT)
