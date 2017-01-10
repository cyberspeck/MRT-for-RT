# -*- coding: utf-8 -*-
"""
Created on Sat Jul 16 22:45:11 2016
@author: david

Volume class
custom FUNcitons using SimpleITK
https://itk.org/Wiki/SimpleITK/GettingStarted#Generic_Distribution

works only with cropped CT and MRI images (showing only one rod),
both Volumes should have the same PixelSpacing,
 and x and y PixelSpacing shoult be equal
sitk_write() creates .mha file with pixel values corresponding
 to distortion in pixel distance * PixelSpacing (mm)

important to remember:
    sitk.Image saves Volume like this (x,y,z)
    array returned by sitk.GetArrayFromImage(Image)
    is transposed: (z,y,x)

based on:
https://pyscience.wordpress.com/2014/10/19/image-segmentation-with-python-and-SimpleITK/
http://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/03_Image_Details.html
http://stackoverflow.com/questions/18435003/ndimages-center-of-mass-to-calculate-the-position-of-a-gaussian-peak
http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.ndimage.measurements.center_of_mass.html

install anaconda
create virtual environment
    conda create -n ENV2 python=2.7 spyder
    source activate ENV2
    pip install SimpleITK
    pip install numpy
    pip install scipy
    pip install matplotlib
    pip install scikit-image
"""

import numpy as np
from scipy import ndimage
import SimpleITK as sitk
import matplotlib.pyplot as plt
import os
from skimage.draw import circle


class Volume:
    '''
    Create a Volume (SimpleITK.Image with convenient properties and functions)
    recommended use:
    create new Volume (optional use denoise=True)
    Volume.getThresholds()

    Parameters
    ----------
    path : string_like
        directory containing DICOM data
    method : string_like, recommended
        either "CT or "MR", used for automatic calculations
    denoise : bool, optional
        If true, the imported data will be denoised using
        SimpleITK.CurvatureFlow(image1=self.img,
                                timeStep=0.125,
                                numberOfIterations=5)
    ref : int, optional
        slice used to make calculations (idealy isocenter) e.g. thresholds
        all plots show this slice
        by default it is set to be in the middle of the image (z-axis)
    info : string, optional
        additional information about imported data, becomes part of title
    seeds : array_like (int,int,int), optional
        coordinates (pixel) of points inside rod, used for segmentation
        by default list of brightest pixel in each slice
    radius: double, optional
        overrides radius value (default CT:4mm, MR:2mm)
    spacing: double, optional
        by default SitpleITK.img.GetSpacing is used to find relation of pixels
        to real length (in mm)
    '''
    def __init__(self, path=None, method=None, denoise=False, ref=None,
                 info=False, seeds='auto', radius=0, spacing=0):
        if(path is None):
            print("Error: no path given!")
        else:
            self.path = path
            self.method = method
            self.denoise = denoise
            self.info = info
            self.centroid = False
            self.mask = False
            self.masked = False
            self.title = method
            self.radius = radius
            self.lower = False
            self.upper = False

            print("\n Import DICOM Files from: ", path)
            self.img = sitk_read(path, self.denoise)

            if (self.img and self.denoise):
                a = self.title
                self.title = a + " denoised"

            if info:
                a = self.title
                self.title = a + ", " + info

            self.xSize, self.ySize, self.zSize = self.img.GetSize()
            if spacing == 0:
                self.xSpace, self.ySpace, self.zSpace = self.img.GetSpacing()

            if type(ref) == int:
                self.ref = ref
            else:
                self.ref = int(self.zSize / 2)

            # niceSlice used to remember which slices show irregularities such
            # as parts of plastic pane (CT)
            # and should therefore not be used to calculate COM, dice, etc.
            self.niceSlice = np.ones((self.zSize, 1), dtype=bool)
            arr = sitk.GetArrayFromImage(self.img)
            average = np.average(arr[ref])
#                print("\nAverage @ ref: ", average)
            for index in range(self.zSize):
                # if average value of slice differs too much -> badSlice
                # difference between ref-Slice and current chosen arbitratry
                if np.absolute(np.average(arr[index]) - average) > 40:
                    print("Irregularities detected in slice {}".format(index))
                    self.niceSlice[index] = False

            if type(seeds) == list:
                self.seeds = seeds
            elif seeds == 'auto':
                self.seeds = []
                for index in range(self.zSize):
                    yMax = int(arr[index].argmax() / self.xSize)
                    xMax = arr[index].argmax() - yMax*self.xSize
                    if self.niceSlice[index] == True:
                        self.seeds.append((xMax, yMax, index))
#                    print("{}: found max at ({},{})".format(index, xMax, yMax))

    def show(self, pixel=False, interpolation=None, ref=None, save=False):
        '''
        plots ref slice of Volume

        Parameters
        ----------
        pixel: bool, optional
            if True, changes axis from mm to pixels
        interpolation: "string", optional, default: 'nearest'
            using build-in interpolation of matplotlib.pyplot.imshow
            Acceptable values are 'none', 'nearest', 'bilinear', 'bicubic',
            'spline16', 'spline36', 'hanning', 'hamming', 'hermite', 'kaiser',
            'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc',
            'lanczos'
        ref: int, optional
            slice to be plotted instead of self.ref (default: 0)
        '''

        if ref is None:
            ref = self.ref

        if interpolation is None:
            a = 'nearest'

        if pixel is False:
            extent = (-self.xSpace/2, self.xSize*self.xSpace - self.xSpace/2, self.ySize*self.ySpace - self.ySpace/2, -self.ySpace/2)
            # The location, in data-coordinates, of the lower-left and upper-right corners
        # (left, right, bottom, top)
        else:
            extent = None

        sitk_show(img=self.img, ref=ref, extent=extent, title=self.title, interpolation=a, save=save)

    def showSeed(self, pixel=False, interpolation='nearest', ref=None, save=False):
        '''
        plots slice containing seed

        Parameters
        ----------
        pixel: bool, optional
            if True, changes axis from mm to pixels
        interpolation: "string", optional, default: 'nearest'
            using build-in interpolation of matplotlib.pyplot.imshow
            Acceptable values are 'none', 'nearest', 'bilinear', 'bicubic',
            'spline16', 'spline36', 'hanning', 'hamming', 'hermite', 'kaiser',
            'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc',
            'lanczos'
        ref: int, optional
            slice of seed to be plotted instead of self.ref (default: zSize/2)
        '''
        if ref is None:
            ref = self.ref
        
        if type(self.seeds[ref]) != tuple:
            print("No seed found @ slice {}".format(ref))
            return None

        extent = None
        if pixel is False:
            extent = (-self.xSpace/2, self.xSize*self.xSpace - self.xSpace/2, self.ySize*self.ySpace - self.ySpace/2, -self.ySpace/2)
            x = (self.seeds[ref][0] * self.xSpace)
            y = (self.seeds[ref][1] * self.xSpace)
        else:
            x, y, z = self.seeds[ref]

        arr = sitk.GetArrayFromImage(self.img)
        fig = plt.figure()
        plt.set_cmap("gray")
        plt.title(self.title + ", seed @ {}".format(self.seeds[ref]))

        plt.imshow(arr[ref, :, :], extent=extent, interpolation=interpolation)
        plt.scatter(x, y)
        plt.show()
        if save != False:
            fig.savefig(str(save) + ".png")

    def getThresholds(self, pixelNumber=0, scale=1):
        '''
        Calculates threshold based on number of pixels representing rod.
        if no pixelNumber is given, self.radius is used to get estimated
        pixelNumber. All calculations based on ref-slice.

        approx. number of pixels being part of rod:
        pn = realRadius^2 * pi / pixelSpacing^2

        Parameters
        ----------
        pixelNumber: int, optional
            if 0, uses self.radis to calculate pixelnumber
            if self.radius also 0, uses self.method instead (CT: 4mm, MR: 2mm)
        scale: double, optional
            factor altering pixelNumber

        Returns
        -------
        Lower and upper threshold value: (double, double)
        '''

        if pixelNumber == 0:
            if self.radius != 0:
                realRadius = self.radius
            else:
                if self.method == "CT":
                    realRadius = 4
                if self.method == "MR":
                    realRadius = 2
                if self.method != "MR" and self.method != "CT":
                    print("method is unknown, please set pixelNumber!")
                    return None
            pixelNumber = np.power(realRadius, 2)*np.pi/np.power(self.xSpace, 2)*scale

        pn = pixelNumber
        arr = sitk.GetArrayFromImage(self.img)
        self.upper = np.double(arr.max())

        hist, bins = np.histogram(arr[self.ref, :, :].ravel(), bins=100)
        self.lower = np.double(bins[np.argmax((np.cumsum(hist[::-1]) < pn)[::-1])])
        print("number of pixels (pn): {}\n lower: {}\n upper: {}".format(pn, self.lower, self.upper))

        return (self.lower, self.upper)

    def getCentroid(self, show=False, percentLimit=False, threshold=False,
                    pixelNumber=0, scale=1, iterations=5, halfShift=0.2,
                    plot=False, save=False):
        '''
        Either used with percentLimit = within (0,1) or "auto"
        or with threshold = number or "auto"
        '''

        if (threshold is False and percentLimit is False) or (threshold is True and percentLimit is True):
            print("Please use either percentLimit or threshold!")
            return None

        if percentLimit == "auto" and threshold is False:
            # value A @ 0.5*(1-halfShift) and B @ 0.5*(1+halfShift) times pn
            # if A > B: next value around 0.25
            # else: around 0.75
            # and so forth
            # calculates 10 centroids with different percentLimits
            # gets dice coefficient for each centroid percentLimit combination
            # returns best result

            arr = sitk.GetArrayFromImage(self.img)
            guess = np.zeros(iterations+1)
            guess[0] = 0.5
            left = 0
            right = 1
            a = 0  # counts how often new iteration has led to smaller score
            thresholdsA = np.zeros((iterations+1,2))
            thresholdsB = np.zeros((iterations,2))
            centroidScoreA = np.zeros(iterations+1)
            centroidScoreB = np.zeros(iterations+1)
            centroidsA = np.zeros((iterations+1, self.zSize, 2))
            centroidsB = np.zeros((iterations+1, self.zSize, 2))
            diceA = np.zeros((iterations+1, self.zSize, 1))
            diceB = np.zeros((iterations+1, self.zSize, 1))
            for index in range(iterations):
                print(" Iteration #{}, A @ {}% = {} pixels".format(index, guess[index]*(1-halfShift)*100, self.xSize*self.ySize*guess[index]*(1-halfShift)))
                thresholdsA[index] = self.getThresholds(pixelNumber=self.xSize*self.ySize*guess[index]*(1-halfShift))
                maskA = sitk.ConnectedThreshold(image1=self.img,
                                                seedList=self.seeds,
                                                lower=self.lower,
                                                upper=self.upper,
                                                replaceValue=1)
                maskedA2 = sitk_applyMask(self.img - arr.min(), maskA)
                maskedA = maskedA2 + arr.min()
                centroidsA[index] = self.xSpace*sitk_centroid(maskedA,
                                                              ref=self.ref,
                                                              threshold=arr.min()+1)
                diceA[index] = self.getDice(centroidsA[index], maskA)
                centroidScoreA[index] = np.average(diceA[index,diceA[index]>-1])
                # sitk_show(maskA, ref= self.ref)
                # self.centroid = centroidsA[index]
                # print("threshold: {} \ncentroid[{}]: {}".format(arr.min()+1, self.ref, self.centroid[self.ref]))
                # self.showCentroid()

                print("\n Iteration #{}, B @ {}% = {} pixels".format(index, guess[index]*(1+halfShift)*100, self.xSize*self.ySize*guess[index]*(1+halfShift)))
                thresholdsB[index] = self.getThresholds(pixelNumber=self.xSize*self.ySize*guess[index]*(1+halfShift))
                maskB = sitk.ConnectedThreshold(image1=self.img,
                                                seedList=self.seeds,
                                                lower=self.lower,
                                                upper=self.upper,
                                                replaceValue=1)
                maskedB2 = sitk_applyMask(self.img - arr.min(), maskB)
                maskedB = maskedB2 + arr.min()
                centroidsB[index] = self.xSpace*sitk_centroid(maskedB,
                                                              ref=self.ref,
                                                              threshold=arr.min()+1)
                diceB[index] = self.getDice(centroidsB[index], maskB)
                centroidScoreB[index] = np.average(diceB[index,diceB[index]>-1])
                # sitk_show(maskB, ref= self.ref)
                # self.centroid = centroidsB[index]
                # print("threshold: {} \ncentroid[{}]: {}".format(arr.min()+1, self.ref, self.centroid[self.ref]))
                # self.showCentroid()
                print("--------------------------")

                if centroidScoreA[index] < centroidScoreB[index]:
                    guess[index+1] = (right + guess[index]) / 2
                    left = guess[index]
                    print("current guess = {}".format((guess[index])))
                elif centroidScoreA[index] >= centroidScoreB[index]:
                    right = guess[index]
                    guess[index+1] = (left + guess[index]) / 2
                    print("current guess = {}".format((guess[index])))
                else:
                    break

                # this checks if new iterations get lower scores: 
                if np.array([centroidScoreA, centroidScoreB]).max() > np.array([centroidScoreA[index], centroidScoreB[index]]).max() and a==0:
                    a = 1
                elif np.array([centroidScoreA, centroidScoreB]).max() > np.array([centroidScoreA[index], centroidScoreB[index]]).max() and a==1:
                    a = 2
                elif np.array([centroidScoreA, centroidScoreB]).max() > np.array([centroidScoreA[index], centroidScoreB[index]]).max() and a==2:
                    # skips to last iteration, and sets guess yielding best
                    # result so far
                    if centroidScoreA.max() >= centroidScoreB.max():
                        guess[iterations] = guess[centroidScoreA.argmax()]
                    else:
                        guess[iterations] = guess[centroidScoreB.argmax()]
                    print("\n\nIteration found region yielding acceptable result, skipping right to")
                    break
                else:
                    a = 0

                print("next guess (#{}) = {} \n \n".format(index+1, guess[index+1]))

            print(" Final iteration (#{}) @ {}% = {} pixels".format(iterations, guess[iterations]*100, self.xSize*self.ySize*guess[iterations]))
            thresholdsA[iterations] = self.getThresholds(pixelNumber=self.xSize*self.ySize*guess[iterations])
            maskA = sitk.ConnectedThreshold(image1=self.img,
                                            seedList=self.seeds,
                                            lower=self.lower,
                                            upper=self.upper,
                                            replaceValue=1)
            maskedA2 = sitk_applyMask(self.img - arr.min(), maskA)
            maskedA = maskedA2 + arr.min()
            centroidsA[iterations] = self.xSpace*sitk_centroid(maskedA,
                                                               ref=self.ref,
                                                               threshold=arr.min()+1)
            centroidScoreA[iterations] = self.getDice(centroidsA[iterations], maskA)
            # sitk_show(maskA, ref= self.ref)
            # self.centroid = centroidsA[iterations]
            # print("threshold: {} \ncentroid[{}]: {}".format(arr.min()+1, self.ref, self.centroid[self.ref]))
            # self.showCentroid()

        
            if centroidScoreA.max() > centroidScoreB.max():
                self.centroid = centroidsA[centroidScoreA.argmax()]
                self.lower, self.upper = thresholdsA[centroidScoreA.argmax()]
                self.dice = diceA[centroidScoreA.argmax()]
                self.diceAverage = centroidScoreA.max()
                print("\nmax dice-coefficient obtained during iteration #{}: {}".format(centroidScoreA.argmax(), centroidScoreA.max()))
            elif (centroidScoreA.max() <= centroidScoreB.max() and centroidScoreB.max() != 0):
                self.centroid = centroidsB[centroidScoreB.argmax()]
                self.lower, self.upper = thresholdsB[centroidScoreB.argmax()]
                self.dice = diceB[centroidScoreB.argmax()]
                self.diceAverage = centroidScoreB.max()
                print("\nmax dice-coefficient obtained during iteration #{}: {}".format(centroidScoreB.argmax(), centroidScoreB.max()))
            else:
                return None
                
            print("\n\n")
            for index in range(np.size(guess)-1):
                print("\niteration #{}".format(index))
                print("A: {}, Score: {}".format(guess[index]*(1-halfShift)*100, centroidScoreA[index]))
                print("B: {}, Score: {}".format(guess[index]*(1+halfShift)*100, centroidScoreB[index]))
            print("\niteration #{}".format(iterations))
            print("A: {}, Score: {}\n".format(guess[iterations]*100, centroidScoreA[iterations]))


            if plot == True:
                percent = np.linspace(0,100, 101, endpoint=True)
                fig = plt.figure()
                for index in range(iterations):
                    if guess[index] > 0 and centroidScoreA[index] > 0:
                        plt.plot(guess[index]*(1-halfShift)*100, centroidScoreA[index], 'bo')
                    if guess[index] > 0 and centroidScoreB[index] > 0:
                        plt.plot(guess[index]*(1+halfShift)*100, centroidScoreB[index], 'go')                
                plt.show()
                if save != False:
                    fig.savefig(str(save) + ".png")
                
        if percentLimit != "auto" and percentLimit is True:
            self.centroid = self.xSpace * sitk_centroid(self.img, ref=self.ref, show=show,
                                                        percentLimit=percentLimit,
                                                        title=self.title)

        if threshold == 'auto':
            self.getThresholds(pixelNumber=pixelNumber, scale=scale)
            self.centroid = self.xSpace * sitk_centroid(self.img, ref=self.ref,
                                                        show=show,
                                                        threshold=self.lower,
                                                        title=self.title)

        if threshold != "auto" and threshold is True:
            self.centroid = self.xSpace * sitk_centroid(self.img, ref=self.ref, show=show,
                                                        threshold=threshold,
                                                        title=self.title)

        for index in range(self.zSize):
            if not self.niceSlice[index]:
                self.centroid[index] = -1, -1
            if self.centroid[index,0] < 0 or self.centroid[index,1] < 0 :
                self.centroid[index] = -1
        return self.centroid


    def showCentroid(self, com2=None, title=None, pixel=False,
                     interpolation='nearest', ref=None, save=False):
        
        if self.centroid is False:
            print("Volume has no centroid yet. use Volume.getCentroid() first!")
            return None

        if title is None:
            title = self.title
        if ref is None:
            ref = self.ref

        if pixel is False:
            extent = (-self.xSpace/2, self.xSize*self.xSpace - self.xSpace/2, self.ySize*self.ySpace - self.ySpace/2, -self.ySpace/2)
            centroid_show(img=self.img, com=self.centroid, com2=com2,
                          extent=extent, save=save, title=title,
                          interpolation=interpolation, ref=ref)
        else:
            centroid_show(img=self.img, com=self.centroid/self.xSpace,
                          com2=com2/self.xSpace, save=save, title=title,
                          interpolation=interpolation, ref=ref)

    def getMask(self, lower=False, upper=False):

        if lower is False and self.lower is not False:
            lower = self.lower
        if upper is False and self.upper is not False:
            upper = self.upper

        if lower is False:
            print("Lower threshold missing!")
            return None
        if upper is False:
            print("Upper threshold missing!")
            return None

        self.mask = sitk_getMask(self.img, self.seeds, upper, lower)
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

    def showMask(self, interpolation=None, ref=None, save=False):
        if self.mask is False:
            print("Volume has no mask yet. use Volume.getMask() first!")
            return None

        if ref is None:
            ref = self.ref

        if interpolation is None:
            interpolation = 'nearest'

        title = self.title + ", mask"

        sitk_show(img=self.mask, ref=ref, title=title,
                  interpolation=interpolation, save=save)

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

    def getDice(self, centroid=None, mask=None, show=False, showAll=False):
        '''
        calculates average dice coefficient ('dc') of the rod by trying
        different radii (neglecting dc of -1)
        returns max obtained average dc after trying different radii

        uses xSpace (should be same as ySpace) to be used with any resolution
        calculating dc for rods with 1.5mm < r < 4mm
        '''
        if centroid is None:
            centroid = self.centroid
        if mask is None:
            mask = self.mask
        
        com = centroid / self.xSpace        

        if self.radius == 0:
            if self.method == "CT":
                radii = np.linspace(2.5, 5.5, num=25)/self.xSpace
            if self.method == "MR":
                radii = np.linspace(0.5, 3.5, num=25)/self.xSpace
            if self.method != "CT" and self.method != "MR":
                # radii = np.linspace(1.5, 4.5, num = 11)*self.xSpace
                print("Unknown method!")
                return None

            dcs = np.zeros(len(radii))
            for index, r in enumerate(radii, start=0):
                dice = dice_circle(img=mask, centroid=com, radius=r, show=showAll)
                dcs[index] = np.average(dice[dice>-1])

            self.dice = dice_circle(img=mask, centroid=com,
                                    radius=radii[dcs.argmax()], show=show)
            print("max dice-coefficient obtained for {} when compared to circle with radius = {}".format(self.method, radii[dcs.argmax()]*self.xSpace))
            print("max dice-coefficient average for the whole volume is: {}".format(dcs.max()))
#            self.radius = radii[dcs.argmax()]
            return np.average(dcs.max())

        else:
            self.dice = dice_circle(img=mask, centroid=com,
                                    radius=self.radius, show=show)
            self.diceAverage = np.average(self.dice)
            print("dice-coefficient average for the whole volume is: {}".format(self.diceAverage))
            return self.dice


def sitk_read(directory, denoise=False):
    '''
    returns DICOM files as "SimpleITK.Image" data type (3D)
    if denoise is True: uses SimpleITK to denoise data
    '''
    reader = sitk.ImageSeriesReader()
    filenames = reader.GetGDCMSeriesFileNames(directory)
    reader.SetFileNames(filenames)
    if denoise:
        print("\n...denoising...")
        imgOriginal = reader.Execute()
        return sitk.CurvatureFlow(image1=imgOriginal,
                                  timeStep=0.125,
                                  numberOfIterations=5)
    else:
        return reader.Execute()


def sitk_write(image, output_dir='', filename='3DImage.mha'):
    '''
    saves image as .mha file
    '''
    output_file_name_3D = os.path.join(output_dir, filename)
    sitk.WriteImage(image, output_file_name_3D)


def sitk_show(img, ref=0, extent=None, title=None, interpolation='nearest', save=False):
    """
    shows plot of img at z=ref
    """
    arr = sitk.GetArrayFromImage(img[:, :, ref])
    fig = plt.figure()
    plt.set_cmap("gray")

    if title:
        plt.title(title)

    plt.imshow(arr, extent=extent, interpolation=interpolation)
    plt.show()
    if save != False:
        fig.savefig(str(save) + ".png")

def sitk_centroid(img, show=False, ref=False, percentLimit=False, save=False,
                  threshold=False, interpolation='nearest', title=None):
    '''
    returns array with y&x coordinate of centroid for every slice of img
    centroid[slice, y&x-coordinate]
    if no pixel has value > threshold:
        centroid x&y-coordinate of that slice = -1,-1
    '''
    if (threshold is False and percentLimit is False) or (threshold is True and percentLimit is True):
        print("Please set either percentLimit or threshold!")
        return None

    arr = sitk.GetArrayFromImage(img)
    z, y, x = np.shape(arr)
    # create array with centroid coordinates of rod in each slice
    com = np.zeros((z, 2))

    if ref is False:
        ref = int(z/2)

    if threshold is False:
        hist, bins = np.histogram(arr[ref, :, :].ravel(),
                                  density=True, bins=100)
        threshold = bins[np.concatenate((np.array([0]), np.cumsum(hist))) *
                         (bins[1] - bins[0]) > percentLimit][0]

    for index in range(z):
        if arr[index].max() > threshold:
            # structuring_element=[[1,1,1],[1,1,1],[1,1,1]]
            segmentation, segments = ndimage.label(arr[index] > threshold)
            # print("segments: {}".format(segments))
            # add ', structuring_element' to label() for recognising
            # diagonal pixels as part of object
            com[index, ::-1] = ndimage.center_of_mass(arr[index, :, :]-threshold,
                                                      segmentation)
            # add ', range(1,segments)' to center_of_mass for list of centroids
            # in each slice (multiple rods!)
        else:
            com[index] = (-1,-1)

    if show:
        if type(show) == bool:
            show == ref
            centroid_show(img, com=com, title=title, save=save,
                          interpolation=interpolation, ref=show)

    return com

def centroid_show(img, com, com2=None, extent=None, title=None, save=False,
                  interpolation='nearest', ref=1):
        arr = sitk.GetArrayFromImage(img)
        fig = plt.figure()
        plt.set_cmap("gray")
        if title:
            plt.title(title + ", centroid")
        x = y = 0
        plt.imshow(arr[ref, :, :], extent=extent, interpolation=interpolation)
        if com2 == None:
            x, y = com[ref]
        else:
            x = [com[ref,0],com2[ref,0]]
            y = [com[ref,1],com2[ref,1]]
        plt.scatter(x, y, c=['b','g'])
        plt.show()
        if save != False:
            fig.savefig(str(save) + ".png")

def coordShift(first, second):
    '''
    returns array with difference of y&x coordinates for every
    centroid[slice, y&x-coordinate]
    '''
    if (np.shape(first) == np.shape(second) and
            np.shape((np.shape(first))) == (2,)):
        z, xy = np.shape(first)
        diff = np.zeros((z, 2))
        for slice in range(z):
            if first[slice,0]==-1 or first[slice,1]==-1 or second[slice,0]==-1 or second[slice,0]==-1:
                diff[slice, 0] = diff[slice, 1] = -1
            else:
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
        if shift[slice,0] == -1 or shift[slice,1] == -1:
            dist[slice,:] =-1
        else:
            dist[slice, :] = np.linalg.norm(shift[slice, :])
    return dist


def sitk_getMask(img, seedList, upper, lower):

    if seedList is False:
        print("no seeds given!")
        return None
        
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
            replaceArray):
        for slice in range(zSize):
            for x in range(xSize):
                for y in range(ySize):
                    if maskA[slice, y, x] == 1:
                        imgMaskedA[slice, y, x] = 1000*replaceArray[slice]*spacing

    return sitk.GetImageFromArray(imgMaskedA)


def dice_circle(img, centroid, radius=2.1, show=False,
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
        created circles. It ranges from 0 (no overlap) to 1 (perfect overlap).
        if centroid coordinates + radius would create circle exceeding image
        size: dc of this slice = -1
        Other errors occuring during the calculation should also result in -1
    """

    xSize, ySize, zSize = img.GetSize()
    profile = np.zeros((zSize, ySize, xSize), dtype=np.uint8)
    centres = centroid.astype(int)    
    dc = np.zeros((zSize, 1))
    for slice in range(zSize):
        if centres[slice,0]+radius < xSize and centres[slice, 1]+radius < ySize and centres[slice,0]-radius > 0 and centres[slice, 1]-radius > 0:

            rr, cc = circle(centres[slice, 0], centres[slice, 1], radius, (xSize,ySize))
            profile[slice, cc, rr] = 1
        else:
            # print("something's fishy!")
            dc[slice]= -1

    input = sitk.GetArrayFromImage(img)

    input = np.atleast_1d(input.astype(np.bool))
    reference = np.atleast_1d(profile.astype(np.bool))

    intersection = np.zeros((zSize, 1))
    size_input = np.zeros((zSize, 1))
    size_reference = np.zeros((zSize, 1))
    for slice in range(zSize):
        intersection[slice] = np.count_nonzero(input[slice, :, :] & reference[slice, :, :])
        size_input[slice] = np.count_nonzero(input[slice, :, :])
        size_reference[slice] = np.count_nonzero(reference[slice, :, :])

        try:
            if (dc[slice] == 0) and (float(size_input[slice] + size_reference[slice]) != 0):
                dc[slice] = 2. * intersection[slice] / float(size_input[slice] + size_reference[slice])

        except ZeroDivisionError:
            dc[slice] = -1

    if show:
        plt.set_cmap("gray")
        plt.title("profile, radius: {}".format(radius))

        plt.imshow(profile[show, :, :], interpolation=interpolation)
        plt.scatter(*centres[show, :])
        plt.show()

    return dc


# to view in 3D Slicer, type this in IPython console or in jupyter notebook:
# %env SITK_SHOW_COMMAND /home/david/Downloads/Slicer-4.5.0-1-linux-amd64/Slicer
# sitk.Show(imgFillingCT)
