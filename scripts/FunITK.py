# -*- coding: utf-8 -*-
"""
Created on Sat Jul 16 22:45:11 2016

useful FUNcitons to using SimpleITK for working with DICOM images

@author: david
based on:
https://pyscience.wordpress.com/2014/10/19/image-segmentation-with-python-and-SimpleITK/
http://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/03_Image_Details.html
"""

# important to remember:
# sitk.Image saves Volume like this (x,y,z)
# array returned by sitk.GetArrayFromImage(Image)
# is transposed: (z,y,x)

import numpy as np
from scipy import ndimage
import SimpleITK as sitk
import matplotlib.pyplot as plt


def sitk_read(directory, denoise=False):
    '''
    returns DICOM files as SimpleITK image
    '''
    reader = sitk.ImageSeriesReader()
    filenames = reader.GetGDCMSeriesFileNames(directory)
    reader.SetFileNames(filenames)
    if denoise:
        return reader.Execute()
    else:
        imgOriginal = reader.Execute()
        return sitk.CurvatureFlow(image1=imgOriginal,
                                  timeStep=0.125,
                                  numberOfIterations=5)


def sitk_show(img, title=None, margin=0.05, dpi=40, scale=3,
              interpolation='nearest'):
    """
    scale is a scaling factor for the shown image
    """
    arr = sitk.GetArrayFromImage(img)

#    spacing = img.GetSpacing()
#    print(title, " Spacing: ", spacing)
#    figsize = (scale + margin) * arr.shape[0] / dpi, (scale + margin) * arr.shape[1] / dpi
#    extent = (0, arr.shape[1]*spacing[1], arr.shape[0]*spacing[0], 0)
#    fig = plt.figure(figsize=figsize, dpi=dpi)
#    ax = fig.add_axes([margin, margin, scale - 2*margin, scale - 2*margin])
    plt.set_cmap("gray")
#    ax.imshow(arr, extent=extent, interpolation=interpolation)

    if title:
        plt.title(title)
    
    plt.imshow(arr, origin="lower", interpolation=interpolation)
    plt.show()


def sitk_centroid(img, show=False, percentLimit=0.9, interpolation='nearest', title=None):
    '''
    returns array with y&x coordinate of centroid for every slice of img
    centroid[slice, y&x-coordinate]
    '''
#    xSpace, ySpace, zSpace = img.GetSpacing()
    arr = sitk.GetArrayFromImage(img)
    z, y, x = np.shape(arr)
    # create array with centroid coordinates of rod in each slice
    com = np.zeros((z, 2))
    for slice in range(z):
        hist, bins = np.histogram(arr[slice,:,:].ravel(), density=True, bins=100)
        threshold = bins[np.cumsum(hist) * (bins[1] - bins[0]) > percentLimit][0]
        marr = np.ma.masked_less(arr[slice,:,:],threshold)
        com[slice,:] = ndimage.measurements.center_of_mass(marr)

    if show:
        plt.set_cmap("gray")
        if title:
            plt.title(title)

        plt.imshow(arr[show,:,:], origin="lower", interpolation=interpolation)
        plt.scatter(*com[slice,::-1])
        plt.show()
 #       centroid_int = centroid.astype(int)
 #       arrCentroid = np.zeros((z, y, x))
 #       for slice in range(z):
 #           arrCentroid[slice, centroid_int[slice, 0], centroid_int[slice, 1]] = 1
 #       imgCentroid = sitk.GetImageFromArray(arrCentroid)
 #       imgCentroid.SetSpacing(img.GetSpacing())
 #       sitk_show(imgCentroid[:,:,show], title="centroid")

    return com

def coordShift(first, second):
    if np.shape(first) == np.shape(second):
        z, xy = np.shape(first)
        diff = np.zeros((z, 2))
        for slice in range(z):
            diff[slice, 0] = first[slice, 0] - second[slice, 0]
            diff[slice, 1] = first[slice, 1] - second[slice, 1]

    return diff


def sitk_mask(img, mask):
    arr = sitk.GetArrayFromImage(img)
    maskA = sitk.GetArrayFromImage(mask)

    imgMaskedA= arr*maskA

    return sitk.GetImageFromArray(imgMaskedA)


# to view in 3D Slicer, type this in IPython console or in jupyter notebook:
# %env SITK_SHOW_COMMAND /home/david/Downloads/Slicer-4.5.0-1-linux-amd64/Slicer
# sitk.Show(imgFillingCT)
