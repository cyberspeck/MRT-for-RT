# -*- coding: utf-8 -*-
"""
Created on Sat Jul 16 22:45:11 2016

Volume class
custom FUNcitons using SimpleITK

@author: david
based on:
https://pyscience.wordpress.com/2014/10/19/image-segmentation-with-python-and-SimpleITK/
http://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/03_Image_Details.html
http://stackoverflow.com/questions/18435003/ndimages-center-of-mass-to-calculate-the-position-of-a-gaussian-peak
http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.ndimage.measurements.center_of_mass.html

install SimpleITK
https://itk.org/Wiki/SimpleITK/GettingStarted#Generic_Distribution
"""

# important to remember:
# sitk.Image saves Volume like this (x,y,z)
# array returned by sitk.GetArrayFromImage(Image)
# is transposed: (z,y,x)

import numpy as np
from scipy import ndimage
import SimpleITK as sitk
import matplotlib.pyplot as plt


class Volume:
    '''
    SimpleITK.Image with convenient properties and functions
    '''
    def __init__(self, path=None, method=None, denoise=False, ref=1,
                 info=False, seeds=None):
        if(path is None):
            print("Error: no path given!")
        else:
            self.path = path
            self.method = method
            self.denoise = denoise
            self.ref = ref
            self.info = info
            self.centroid = False
            self.mask = False
            self.masked = False
            if seeds:
                self.seeds=seeds

            print("Import DICOM Files from: ", path)
            self.img = sitk_read(path, denoise)

            if (self.img is True and self.denoise is True):
                print("\n...denoising...")

            self.xSpace, self.ySpace, self.zSpace = self.img.GetSpacing()
            self.xSize, self.ySize, self.zSize = self.img.GetSize()
            self.title = method

            if denoise is True:
                a = self.title
                self.title = a + " denoised"

            if info is True:
                a = self.title
                self.title = a + ", " + info

    def show(self, interpolation=None, ref=None):
        if ref is None:
            ref = self.ref

        if interpolation is None:
            a = 'nearest'

        sitk_show(img=self.img, ref=ref, title=self.title, interpolation=a)
        
    def showSeed(self, title=None, interpolation='nearest'):
        if self.seeds is False:
            print("Volume has no seeds yet.")
            return None

        x,y,z = self.seeds[0]
        arr = sitk.GetArrayFromImage(self.img)
        plt.set_cmap("gray")
        if title is None:
            plt.title(self.title + ", seed")

        plt.imshow(arr[z, :, :], interpolation=interpolation)
        plt.scatter(x,y)
        plt.show()

    def getCentroid(self, show=False, percentLimit=0.9,
                      title=None, method=None, new=False):
        if (self.centroid is not False and new is False):
            return self.centroid

        if title is None:
            title = self.title

        if method is None:
            method = self.method

        self.centroid = sitk_centroid(self.img, show=show,
                                      percentLimit=percentLimit,
                                      title=title, method=method)
        return self.centroid        
        
    def showCentroid(self, title=None, interpolation='nearest', ref=None):
        if self.centroid is False:
            print("Volume has no centroid yet. use Volume.getCentroid() first!")
            return None

        if title is None:
            title = self.title
        if ref is None:
            ref = self.ref
        centroid_show(img=self.img, com=self.centroid, title=title,
                      interpolation=interpolation, ref=ref)
        

    def getMask(self, lower=None, upper=None, replaceValue=1, new=False):
        if self.mask and new == False:
            return self.mask

        if self.method is "CT":
            if lower is None:
                lower = 0
            if upper is None:
                upper = 300
                
        if self.method is "MR":
            if lower is None:
                lower = 80
            if upper is None:
                upper = 1500            
 
        self.mask = sitk_getMask(img=self.img, seedList=self.seeds,
                                 lower=lower, upper=upper,
                                 replaceValue=replaceValue)
        return self.mask

    def applyMask(self, mask=None, replaceArray=False):
        if (mask is None):
            if self.mask:
                mask = self.mask
            else:
                print("Volume has no mask yet. use Volume.getMask() first!")
                return None
                        
        self.masked = sitk_applyMask(self.img, mask, replaceArray=replaceArray)
        
        return self.masked
        
    def showMask(self, interpolation=None, ref=None):
        if self.mask is False:
            print("Volume has no mask yet. use Volume.getMask() first!")
            return None

        if ref is None:
            ref = self.ref

        if interpolation is None:
            interpolation = 'nearest'

        title = self.title + ", mask"

        sitk_show(img=self.mask, ref=ref, title=title,
                  interpolation=interpolation)
        
    def showMasked(self, interpolation=None, ref=None):
        if self.masked is False:
            print ("Volume has not been masked yet. use Volume.applyMask() first!")
            return None
        if ref is None:
            ref = self.ref

        if interpolation is None:
            interpolation = 'nearest'

        title = self.title + ", mask"

        sitk_show(img=self.masked, ref=ref, title=title,
                  interpolation=interpolation)

def sitk_read(directory, denoise=False):
    '''
    returns DICOM files as "SimpleITK.Image" data type (3D)
    if denoise is true: uses SimpleITK to denoise data
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


def sitk_show(img, ref=1, title=None, interpolation='nearest'):
    """
    shows plot of img at z=ref
    """
    arr = sitk.GetArrayFromImage(img[:, :, ref])
    plt.set_cmap("gray")

    if title:
        plt.title(title)

    plt.imshow(arr, interpolation=interpolation)
    plt.show()


def sitk_centroid(img, show=False, percentLimit=0.9, interpolation='nearest',
                  title=None, method=None):
    '''
    returns array with y&x coordinate of centroid for every slice of img
    centroid[slice, y&x-coordinate]
    
    if method is None: compute centroid using only brightest 10% of pixels
    if method is CT or MR: using pixels above thershold
        for CT threshold = -900
        for MR threshold = 30
    '''
    arr = sitk.GetArrayFromImage(img)
    z, y, x = np.shape(arr)
    # create array with centroid coordinates of rod in each slice
    com = np.zeros((z, 2))

    if method == "CT":
        threshold = -900 
            # Set Threshold at -900 (air starts at -950)

    if method == "MR":
        threshold = 30
            # Set Threshold at 30 (air level in 3D slicer)

    if (method != "MR" and method != "CT"):
        for slice in range(z):
            hist, bins = np.histogram(arr[slice, :, :].ravel(),
                                      density=True, bins=100)
            threshold = bins[np.cumsum(hist) * (bins[1] - bins[0])
                             > percentLimit][0]
            marr = np.ma.masked_less(arr[slice, :, :], threshold)
            com[slice, ::-1] = ndimage.measurements.center_of_mass(marr)
    else:        
        for slice in range(z):
            marr = np.ma.masked_less(arr[slice, :, :], threshold)
            com[slice, ::-1] = ndimage.measurements.center_of_mass(marr)

    if show:
        centroid_show(img, com=com, title=title,
                      interpolation=interpolation, ref=show)

    return com

def centroid_show(img, com, title=None, interpolation='nearest', ref=1):
        arr = sitk.GetArrayFromImage(img)
        plt.set_cmap("gray")
        if title:
            plt.title(title + ", centroid")

        plt.imshow(arr[ref, :, :], interpolation=interpolation)
        plt.scatter(*com[ref, :])
        plt.show()

def coordShift(first, second):
    '''
    returns array with difference of y&x coordinates for every
    centroid[slice, y&x-coordinate]

    if method is None: compute centroid using only brightest 10% of pixels
    if method is CT or MR: using pixels above thershold
        for CT threshold = -900
        for MR threshold = 30
    '''
    if (np.shape(first) == np.shape(second) and
            np.shape((np.shape(first))) == (2,)):
        z, xy = np.shape(first)
        diff = np.zeros((z, 2))
        for slice in range(z):
            diff[slice, 0] = first[slice, 0] - second[slice, 0]
            diff[slice, 1] = first[slice, 1] - second[slice, 1]
        return diff
    else:
        print("Wrong shape! coordShift returned 'False'")
        return False

def sitk_getMask(img, seedList, upper, lower, replaceValue=1):
    return sitk.ConnectedThreshold(image1=img, seedList=seedList, 
                                   lower=lower, upper=upper,
                                   replaceValue=replaceValue)

def sitk_applyMask(img, mask, replaceArray=False):
    '''
    masks img (SimpleITK.Image) using mask (SimpleITK.Image)
    '''
    if img.GetSize() != mask.GetSize():
        print(mask.GetSize())
        print(img.GetSize())
        
        print("mask and image are not the same size!")
        return False
    
    arr = sitk.GetArrayFromImage(img)
    maskA = sitk.GetArrayFromImage(mask)
    xSize, ySize, zSize = img.GetSize()

    imgMaskedA = arr*maskA
    
    if (np.shape(replaceArray) == (img.GetDepth(), 1) and
     replaceArray is not False):
        for slice in range(zSize):
            for x in range(xSize):
                for y in range(ySize):
                    if maskA[slice, y, x] == 1:
                        #print("a: ", mask[slice, y, x], "type: ", type(mask[slice, y, x]))
                        #print("b: ", np.uint8(1000*replaceArray[slice]), "type: ", type(np.uint8(1000*replaceArray[slice])))
                        imgMaskedA[slice, y, x] = 1000*replaceArray[slice]
                        #print("c: ", mask[slice, y, x]), "type: ", type(mask[slice, y, x])

    return sitk.GetImageFromArray(imgMaskedA)


# to view in 3D Slicer, type this in IPython console or in jupyter notebook:
# %env SITK_SHOW_COMMAND /home/david/Downloads/Slicer-4.5.0-1-linux-amd64/Slicer
# sitk.Show(imgFillingCT)
