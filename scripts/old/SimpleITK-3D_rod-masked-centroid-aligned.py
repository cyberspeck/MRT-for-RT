# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 15:17:57 2016

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

pathCT = "../data/cropped_CT-a/"
pathMR = "../data/cropped_MR-d-a/"

imgOriginalCT = fun.sitk_read(pathCT)
imgOriginalMR = fun.sitk_read(pathMR)

idxSlice = 10

labelRod = 1

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



seedFillingCT = [(8,8,idxSlice)]
seedFillingMR = [(8,8,idxSlice)]

# https://itk.org/SimpleITKDoxygen/html/classitk_1_1simple_1_1ConnectedThresholdImageFilter.html

maskRodCT = sitk.ConnectedThreshold(image1=imgSmoothCT, 
                                              seedList=seedFillingCT, 
                                              lower=00, 
                                              upper=300,
                                              replaceValue=labelRod)

maskRodMR = sitk.ConnectedThreshold(image1=imgSmoothMR, 
                                              seedList=seedFillingMR, 
                                              lower=80, 
                                              upper=1500,
                                              replaceValue=labelRod)

fun.sitk_show(maskRodCT[:,:,idxSlice], title="CT, denoised, rod mask")
fun.sitk_show(maskRodMR[:,:,idxSlice], title="MR, denoised, rod mask")

imgRodCT = fun.sitk_mask(imgSmoothCT, maskRodCT)
imgRodMR = fun.sitk_mask(imgSmoothMR, maskRodMR)

fun.sitk_show(imgRodCT[:,:,idxSlice], title="CT, denoised, rod")
fun.sitk_show(imgRodMR[:,:,idxSlice], title="MR, denoised, rod")

# important to remember:
# sitk.Image saves Volume like this (x,y,z)
# array returned by sitk.GetArrayFromImage(Image)
# is transposed: (z,y,x)

rodCT = sitk.GetArrayFromImage(imgRodCT)
zCT, yCT, xCT = np.shape(rodCT)

rodMR = sitk.GetArrayFromImage(imgRodMR)
zMR, yMR, xMR = np.shape(rodMR)

# create array with centroid of rod in each slice of CT and MRT
centroidCT = np.zeros((zCT, 2))
centroidMR = np.zeros((zMR, 2))

for slice in range(zCT):
    centroidCT[slice,:] = np.array(ndimage.measurements.center_of_mass(rodCT[slice,:,:]))

for slice in range(zMR):
    centroidMR[slice,:] = np.array(ndimage.measurements.center_of_mass(rodMR[slice,:,:]))

centroidCT_int = centroidCT.astype(int)

imgCentroidCT = np.zeros((zCT, yCT, xCT))
for slice in range(zCT):
    imgCentroidCT[slice, centroidCT_int[slice,0], centroidCT_int[slice,1]] = 1


sitkCentroidCT = sitk.GetImageFromArray(imgCentroidCT)
fun.sitk_show(sitkCentroidCT[:,:,idxSlice])

if zCT == zMR:
    centroidDiff = np.zeros((zCT, 2))
    for slice in range(zCT):
        centroidDiff[slice,0] = centroidCT[slice,0] - centroidMR[slice, 0]
        centroidDiff[slice,1] = centroidCT[slice,1] - centroidMR[slice, 1]

# and now without denoising and masking:

rodCT_org = sitk.GetArrayFromImage(imgOriginalCT)
rodMR_org = sitk.GetArrayFromImage(imgOriginalMR)

centroidCT_org = np.zeros((zCT, 2))
centroidMR_org = np.zeros((zMR, 2))

for slice in range(zCT):
    centroidCT_org[slice,:] = np.array(ndimage.measurements.center_of_mass(rodCT_org[slice,:,:]))

for slice in range(zMR):
    centroidMR_org[slice,:] = np.array(ndimage.measurements.center_of_mass(rodMR_org[slice,:,:]))

if zCT == zMR:
    centroidDiff_org = np.zeros((zCT, 2))
    for slice in range(zCT):
        centroidDiff_org[slice,0] = centroidCT_org[slice,0] - centroidMR_org[slice, 0]
        centroidDiff_org[slice,1] = centroidCT_org[slice,1] - centroidMR_org[slice, 1]

# is denoising and masking realy necessary?

centroidDiff_noise = np.zeros((zCT, 2))
for slice in range(zCT):
    centroidDiff_noise[slice,0] = centroidDiff_org[slice,0] - centroidDiff[slice, 0]
    centroidDiff_noise[slice,1] = centroidDiff_org[slice,1] - centroidDiff[slice, 1]



print("denoised centroid ct")
print(centroidCT[:2,:])

print("original centroid ct")
print(centroidCT_org[:2,:])

print("denoised")
print(centroidDiff[:2,:])
print("orginal")
print(centroidDiff_org[:2,:])
print("diff_diff")
print(centroidDiff_noise[:2,:])