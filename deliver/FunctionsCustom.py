
# coding: utf-8

# In[2]:

import numpy as np
from scipy import ndimage
import SimpleITK as sitk
import matplotlib.pyplot as plt


# In[3]:

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

    plt.imshow(arr, origin="lower", interpolation=interpolation)
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
        for slice in range(z):
            # Set Threshold at -900 (air level measured in 3D slicer)
            marr = np.ma.masked_less(arr[slice, :, :], -900)
            com[slice, :] = ndimage.measurements.center_of_mass(marr)

    if method == "MR":
        for slice in range(z):
            # Set Threshold at 30 (air level measured in 3D slicer)
            marr = np.ma.masked_less(arr[slice, :, :], 30)
            com[slice, :] = ndimage.measurements.center_of_mass(marr)

    if (method != "MR" and method != "CT"):
        for slice in range(z):
            hist, bins = np.histogram(arr[slice, :, :].ravel(),
                                      density=True, bins=100)
            threshold = bins[np.cumsum(hist) * (bins[1] - bins[0]) > percentLimit][0]
            marr = np.ma.masked_less(arr[slice, :, :], threshold)
            com[slice, :] = ndimage.measurements.center_of_mass(marr)

    if show:
        plt.set_cmap("gray")
        if title:
            plt.title(title + ", centroid")

        plt.imshow(arr[show, :, :], origin="lower", interpolation=interpolation)
        plt.scatter(*com[slice, ::-1])
        plt.show()

    return com


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


def sitk_mask(img, mask):
    '''
    masks img (SimpleITK.Image) using mask (SimpleITK.Image)
    '''
    arr = sitk.GetArrayFromImage(img)
    maskA = sitk.GetArrayFromImage(mask)

    imgMaskedA = arr*maskA

    return sitk.GetImageFromArray(imgMaskedA)


# to view in 3D Slicer, type this in IPython console or in jupyter notebook:
# %env SITK_SHOW_COMMAND /home/david/Downloads/Slicer-4.5.0-1-linux-amd64/Slicer
# sitk.Show(imgFillingCT)


# In[4]:

class Volume:
    '''
    SimpleITK.Image with convenient properties and functions
    '''
    def __init__(self, path=None, method=None,
                 denoise=False, ref=1, info=False):
        if(path is None):
            print("Error: no path given!")
        else:
            self.path = path
            self.method = method
            self.denoise = denoise
            self.ref = ref
            self.info = info

            print("Import DICOM Files from: ", path)
            self.img = sitk_read(path, denoise)

            if (self.img is True and self.denoise is True):
                print("\n...denoising...")

            self.xSpace, self.ySpace, self.zSpace = self.img.GetSpacing()
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

    def centroid(self, show=False, percentLimit=0.9, title=None, method=None):
        if title is None:
            title = self.title

        if method is None:
            method = self.method

        self.centroid = sitk_centroid(self.img, show=show,
                                      percentLimit=percentLimit,
                                      title=title, method=method)

