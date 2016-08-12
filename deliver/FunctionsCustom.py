
# coding: utf-8

# In[1]:

import numpy as np
from scipy import ndimage
import SimpleITK as sitk
import matplotlib.pyplot as plt
import os
from skimage.draw import circle


# In[2]:

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


def sitk_write(image, output_dir='', filename='3DImage.mha'):
    '''
    saves image as .mha file
    '''
    output_file_name_3D = os.path.join(output_dir, filename)
    sitk.WriteImage(image, output_file_name_3D)


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
            threshold = bins[np.cumsum(hist) * (bins[1] - bins[0]) >
                             percentLimit][0]
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


def coordDist(shift):
    '''
    calculates norm for each entry of array
    returns array with list of calculated values
    '''
    if np.shape(shift)[1] != 2:
        print("shift has wrong shape!")
        return False

    dist = np.zeros((len(shift), 1))
    for slice in range(len(shift)):
        dist[slice, :] = np.linalg.norm(shift[slice, :])
    return dist


def sitk_getMask(img, seedList, upper, lower):

    return sitk.ConnectedThreshold(image1=img, seedList=seedList,
                                   lower=lower, upper=upper,
                                   replaceValue=1)


def sitk_applyMask(img, mask, replaceArray=False, spacing=1):
    '''
    masks img (SimpleITK.Image) using mask (SimpleITK.Image)
    if a replaceArray is given, the spacing*values*1000 of the array will be used
    as pixel intensity for an entire slice each
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
                        imgMaskedA[slice, y, x] = 1000*replaceArray[slice]*spacing

    return sitk.GetImageFromArray(imgMaskedA)


def dice_circle(input_img, centroid, radius=2.1, show=False,
                interpolation='nearest'):
    """
    Dice coefficient, inspired by
     Medpy (http://pythonhosted.org/MedPy/_modules/medpy/metric/binary.html)

    Computes the Dice coefficient (akas Sorensen index) between a binary
    object in an image and a circle.

    The metric is defined as:

        DC=\frac{2|A\cap B|}{|A|+|B|}

    where A is the first and B the second set of samples (here: binary objects)
    sys.argv[2]
    Parameters
    ----------
    input_umg : SimpleITK.Image
        Input data containing objects. Can be any type but will be converted
        into binary: background where 0, object everywhere else.
    centroid : array_like
        array with coordinates for circle centre
    radius : float
        radius for creating reference circles

    Returns
    -------
    dc : array_like
        The Dice coefficient between the object(s) in ```input``` and the
        created cirles. It ranges from 0 (no overlap) to 1 (perfect overlap).
    """

    xSize, ySize, zSize = input_img.GetSize()
    profile = np.zeros((zSize, ySize, xSize), dtype=np.uint8)
    centres = centroid.astype(int)
    for slice in range(zSize):
        rr, cc = circle(centres[slice, 0], centres[slice, 1], radius)
        profile[slice, cc, rr] = 1

    input = sitk.GetArrayFromImage(input_img)

    input = np.atleast_1d(input.astype(np.bool))
    reference = np.atleast_1d(profile.astype(np.bool))

    intersection = np.zeros((zSize, 1))
    size_input = np.zeros((zSize, 1))
    size_reference = np.zeros((zSize, 1))
    dc = np.zeros((zSize, 1))
    for slice in range(zSize):
        intersection[slice] = np.count_nonzero(input[slice, :, :] & reference[slice, :, :])
        size_input[slice] = np.count_nonzero(input[slice, :, :])
        size_reference[slice] = np.count_nonzero(reference[slice, :, :])
#        print("\n intersection[slice, :]: {}".format(intersection[slice, :]))
#        print("size_input[slice]: {}".format(size_input[slice]))
#        print("size_reference[slice]: {}".format(size_reference[slice]))
        try:
            dc[slice] = 2. * intersection[slice] / float(size_input[slice] + size_reference[slice])
        except ZeroDivisionError:
            dc[slice] = 0.0

    if show:
        plt.set_cmap("gray")
        plt.title("profile, radius: {}".format(radius))

        plt.imshow(profile[show, :, :], interpolation=interpolation)
        plt.scatter(*centres[show, :])
        plt.show()

    return dc


# In[3]:

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
                self.seeds = seeds

            print("\n Import DICOM Files from: ", path)
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

        x, y, z = self.seeds[0]
        arr = sitk.GetArrayFromImage(self.img)
        plt.set_cmap("gray")
        if title is None:
            plt.title(self.title + ", seed")

        plt.imshow(arr[z, :, :], interpolation=interpolation)
        plt.scatter(x, y)
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

    def getMask(self, lower=None, upper=None, new=False):
        if self.mask and new is False:
            return self.mask

        if self.method is "CT":
            if lower is None:
                lower = -300
            if upper is None:
                upper = 300

        if self.method is "MR":
            if lower is None:
                lower = 120
            if upper is None:
                upper = 1500

        self.mask = sitk_getMask(img=self.img, seedList=self.seeds,
                                 lower=lower, upper=upper)
        return self.mask

    def applyMask(self, mask=None, replaceArray=False, spacing=1):
        if (mask is None):
            if self.mask:
                mask = self.mask
            else:
                print("Volume has no mask yet. use Volume.getMask() first!")
                return None

        self.masked = sitk_applyMask(self.img, mask, replaceArray=replaceArray,
                                     spacing=spacing)

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
            print("Volume has not been masked yet. use Volume.applyMask() first!")
            return None
        if ref is None:
            ref = self.ref

        if interpolation is None:
            interpolation = 'nearest'

        title = self.title + ", mask"

        sitk_show(img=self.masked, ref=ref, title=title,
                  interpolation=interpolation)

    def getDice(self, show=False, showAll=False):
        '''
        calculates max dice coefficient by trying different radii
        returns radius that leads to best result
        '''
        radii = np.array([2, 2.1, 2.3, 2.9, 3.1, 3.2, 3.7, 4.1, 4.2])
        dcs = np.zeros(len(radii))
        for index, r in enumerate(radii, start=0):
            dcs[index] = np.average(dice_circle(self.mask, self.centroid,
                                    radius=r, show=showAll))

        self.dice = dice_circle(self.mask, self.centroid,
                                radius=radii[dcs.argmax()], show=show)
        print("max dice-coefficient obtained for {} when compared to circle with radius = {}".format(
        self.method, radii[dcs.argmax()]))
        print("max dice-coefficient average for the whole volume is: {}".format(dcs.max()))
        return radii[dcs.argmax()]

