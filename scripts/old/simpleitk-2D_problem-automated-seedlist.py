# -*- coding: utf-8 -*-
"""
Created on Thu Jul 8 2016

@author: david
based on:
https://pyscience.wordpress.com/2014/10/19/image-segmentation-with-python-and-simpleitk/
"""

import os
import numpy
import SimpleITK
import matplotlib.pyplot as plt
get_ipython().magic('pylab inline')


def sitk_show(img, title=None, margin=0.05, dpi=40, scale=2, interpolation=None ):
    """
    scale is a scaling factor for the shown image
    """
    nda = SimpleITK.GetArrayFromImage(img)
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

pathDicom = "../data/cropped_CT/"

idxSlice = 50

labelPlastic = 1
labelFilling = 2


reader = SimpleITK.ImageSeriesReader()
filenamesDICOM = reader.GetGDCMSeriesFileNames(pathDicom)
reader.SetFileNames(filenamesDICOM)
imgOriginal = reader.Execute()


imgOriginal_slice = imgOriginal[:,:,idxSlice]


sitk_show(imgOriginal_slice, title="immer noch Supa")


imgSmooth = SimpleITK.CurvatureFlow(image1=imgOriginal_slice,
                                    timeStep=0.125,
                                    numberOfIterations=5)


'''
x= (12,14)
seedFilling = [(14,14),(14,15),x]
imgFilling = SimpleITK.ConnectedThreshold(image1=imgSmooth, 
                                              seedList=seedFilling, 
                                              lower=00, 
                                              upper=110,
                                              replaceValue=labelFilling)
sitk_show(imgFilling)
'''

connectedFilling = SimpleITK.ConnectedThresholdImageFilter()
connectedFilling.SetLower(00)
connectedFilling.SetUpper(110)
connectedFilling.SetReplaceValue(labelFilling)
connectedFilling.AddSeed((14,14))
connectedFilling.AddSeed((14,15))

imgFilling = connectedFilling.Execute(imgSmooth)
sitk_show(imgFilling)


# -get maximum values of plastic for seedList-
# first make array out of denoised img:
arraySmooth = SimpleITK.GetArrayFromImage(imgSmooth)

# all pixels > 125 used as seeds:
# i,j = where(arraySmooth > 125)
# arraySmooth[i,j] = 1000

maxPlastic = np.array(where(arraySmooth > 125)).T

# show seeds in array    
for x,y in maxPlastic:
    arraySmooth[x,y] = 1000
plt.imshow(arraySmooth)

listPlastic = list(map(tuple, maxPlastic))
print("listPlastic: \n type:", type(listPlastic), "\n value:", listPlastic)

# this doesn't work:
'''
imgPlastic = SimpleITK.ConnectedThreshold(image1=imgSmooth, 
                                              seedList=listPlastic, 
                                              lower=110, 
                                              upper=200,
                                              replaceValue=labelPlastic)
'''

# neither does this:
'''
seedPlastic = SimpleITK.VectorUIntList()
for pixel in listPlastic:
    print("pixel:\n type:", type(pixel), "\n value:", pixel,"\n")
    seed = SimpleITK.VectorUInt32(pixel)
    print("seed:\n type:", type(seed), "\n value:", seed)
    seedPlastic.push_back(seed)
'''

# or this:
'''
connectedPlastic = SimpleITK.ConnectedThresholdImageFilter()
connectedPlastic.SetLower(110)
connectedPlastic.SetUpper(200)
connectedPlastic.SetReplaceValue(labelPlastic)
for pixel in listPlastic:
    print("pixel:\n type:", type(pixel), "\n value:", pixel,"\n")
    connectedPlastic.AddSeed(pixel)

imgPlastic = connectedPlastic.Execute(imgSmooth)
sitk_show(imgPlastic)
'''

