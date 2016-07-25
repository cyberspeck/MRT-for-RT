# -*- coding: utf-8 -*-
"""
Created on Sat Jul 16 22:45:11 2016

useful FUNcitons to using SimpleITK for working with DICOM images

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


def sitk_show(img, title=None, margin=0.05, dpi=40, scale=3,
              interpolation='nearest'):
    """
    scale is a scaling factor for the shown image
    """
    nda = sitk.GetArrayFromImage(img)
    spacing = img.GetSpacing()
    figsize = (scale + margin) * nda.shape[0] / dpi, (scale + margin) * nda.shape[1] / dpi
    extent = (0, nda.shape[1]*spacing[1], nda.shape[0]*spacing[0], 0)
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_axes([margin, margin, scale - 2*margin, scale - 2*margin])

    plt.set_cmap("gray")
    ax.imshow(nda,extent=extent,interpolation=interpolation)
    
    if title:
        plt.title(title)
    
    plt.show()

def array_show(array):
    plt.figure(dpi=40)
    plt.axes().set_aspect('equal', 'datalim')
    plt.set_cmap(plt.gray())
    plt.pcolor(np.flipud(array))

def sitk_mask(img0, mask):
    img0A = sitk.GetArrayFromImage(img0)
    maskA = sitk.GetArrayFromImage(mask)

    imgMaskedA= img0A*maskA
    
    return sitk.GetImageFromArray(imgMaskedA)

def sitk_read(directory):
    '''
    returns Volume as SimpleITK image
    '''
    reader = sitk.ImageSeriesReader()
    filenames = reader.GetGDCMSeriesFileNames(directory)
    reader.SetFileNames(filenames)
    return reader.Execute()
    
# to view in 3D Slicer, type this in IPython console or in jupyter notebook:
# %env SITK_SHOW_COMMAND /home/david/Downloads/Slicer-4.5.0-1-linux-amd64/Slicer
# sitk.Show(imgFillingCT)
