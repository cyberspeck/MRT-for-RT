
# coding: utf-8

# In[1]:

# This notebook is based on the very helpfull blog
# https://pyscience.wordpress.com/2014/10/19/image-segmentation-with-python-and-simpleitk/

import os
import numpy
import SimpleITK
import matplotlib.pyplot as plt
get_ipython().magic('pylab inline')


# In[9]:

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


# In[10]:

# Directory where the DICOM files are being stored (in this
# case the 'data/cropped_CT' folder). 
pathDicom = "../data/cropped_CT/"

# Z slice of the DICOM files to process. In the interest of
# simplicity, segmentation will be limited to a single 2D
# image but all processes are entirely applicable to the 3D image
idxSlice = 50

# int labels to assign to the segmented white and gray matter.
# These need to be different integers but their values themselves
# don't matter
labelPlastic = 1
labelFilling = 2


# In[11]:

reader = SimpleITK.ImageSeriesReader()
filenamesDICOM = reader.GetGDCMSeriesFileNames(pathDicom)
reader.SetFileNames(filenamesDICOM)
imgOriginal = reader.Execute()


# In[12]:

# For now we'll only look at a single 2D image
imgOriginal_slice = imgOriginal[:,:,idxSlice]


# In[13]:

# and look at it
sitk_show(imgOriginal_slice, title="immer noch Supa")


# In[14]:

imgSmooth = SimpleITK.CurvatureFlow(image1=imgOriginal_slice,
                                    timeStep=0.125,
                                    numberOfIterations=5)

# imgSmooth[14,14]= 400
# imgSmooth[15,15]= 0
# get maximum values of plastic:
a = SimpleITK.GetArrayFromImage(imgSmooth)
i,j = where(a>120)
maxPlastic = np.array(where(a>120))
# maxPlastic =numpy.nonzero(a> 120)
a[i,j] = 1000
plt.imshow(a)
#seedPlastic = list(map(tuple, maxPlastic))
seedPlastic = list(map(tuple, maxPlastic.T))

sitk_show(imgSmooth, interpolation="nearest")

#seedPlastic = [(12,14), (15,15)]
seedFilling = [(14,14)]

imgFilling = SimpleITK.ConnectedThreshold(image1=imgSmooth, 
                                              seedList=seedFilling, 
                                              lower=00, 
                                              upper=110,
                                              replaceValue=labelFilling)


sitk_show(imgFilling)
#sitk_show(imgFilling, interpolation="nearest")

imgPlastic = SimpleITK.ConnectedThreshold(image1=imgSmooth, 
                                              seedList=seedPlastic, 
                                              lower=110, 
                                              upper=200,
                                              replaceValue=labelPlastic)


sitk_show(imgPlastic)
#sitk_show(imgPlastic, interpolation="nearest")


print("Seed Plastic Value:", imgSmooth[12,14])
print("Seed Filling Value:", imgSmooth[14,14])

# print(plt.hist((imgSmooth), bins = 40))