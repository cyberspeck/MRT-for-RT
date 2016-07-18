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

pathCT = "../data/cropped_CT/"
pathMR = "../data/cropped_MR-d/"

imgOriginalCT = fun.sitk_read(pathCT)
imgOriginalMR = fun.sitk_read(pathMR)

idxSlice = 20

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



seedFillingCT = [(14,14,idxSlice)]
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
    
fun.array_show(imgCentroidCT[0,:,:])
sitkCentroidCT = sitk.GetImageFromArray(imgCentroidCT)
fun.sitk_show(sitkCentroidCT[:,:,idxSlice])


'''
imgSmoothCT_int = sitk.Cast(sitk.RescaleIntensity(imgSmoothCT[:,:,idxSlice]), sitkCentroidCT[:,:,idxSlice].GetPixelID())
sitkCentroidCT_int = sitk.Cast(sitk.RescaleIntensity(sitkCentroidCT[:,:,idxSlice]), sitkCentroidCT[:,:,idxSlice].GetPixelID())
sitk_show(sitk.LabelOverlay(imgSmoothCT_int, sitkCentroidCT))
'''



'''
#read the images
fixed_image =  sitk.ReadImage('../data/cropped_CT/001-10w-1ct0001.dcm', sitk.sitkFloat32)
moving_image = sitk.ReadImage('../data/cropped_MR-d/001-10w-1mr-d0001.dcm', sitk.sitkFloat32) 
 
#initial alignment of the two volumes
transform = sitk.CenteredTransformInitializer(fixed_image,
                                              moving_image,
                                              sitk.Euler3DTransform(),
                                              sitk.CenteredTransformInitializerFilter.GEOMETRY)

#multi-resolution rigid registration using Mutual Information
registration_method = sitk.ImageRegistrationMethod()
registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
registration_method.SetMetricSamplingPercentage(0.01)
registration_method.SetInterpolator(sitk.sitkLinear)
registration_method.SetOptimizerAsGradientDescent(learningRate=1.0,
                                                  numberOfIterations=100,
                                                  convergenceMinimumValue=1e-6,
                                                  convergenceWindowSize=10)

registration_method.SetOptimizerScalesFromPhysicalShift()
registration_method.SetShrinkFactorsPerLevel(shrinkFactors = [4,2,1])
registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas=[2,1,0])
registration_method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
registration_method.SetInitialTransform(transform)
registration_method.Execute(fixed_image, moving_image)

sitk.WriteTransform(transform, 'ct2mrT1.tfm')
'''