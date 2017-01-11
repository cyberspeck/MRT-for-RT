
# coding: utf-8

# In[1]:

import FunctionsCustom as fun
from FunctionsCustom import Volume
get_ipython().magic('pylab inline')

pathCT = "../data/cropped_CT-a/"
pathMR = "../data/cropped_MR-d-a/"
idxSlice = 10


# In[2]:

CT = Volume(path=pathCT, method="CT", ref=idxSlice, seeds=[(6, 8, idxSlice)])
MR = Volume(path=pathMR, method="MR", ref=idxSlice, seeds=[(6, 8, idxSlice)])


# In[3]:

CT.show()
MR.show()


# In[4]:

CT.getCentroid()
CT.showCentroid()
MR.getCentroid()
MR.showCentroid()


# In[5]:

# this calculates the coordinate difference of MR.centroid relative to CT.centroid
distortion = fun.coordShift(CT.centroid, MR.centroid)

print("CT, centroid[:,:,{}]: {}".format(idxSlice, CT.centroid[idxSlice]))
print("MR, centroid[:,:,{}]: {}".format(idxSlice, MR.centroid[idxSlice]))
print("distrotion[:,:,{}]: {}".format(idxSlice, distortion[idxSlice]))


# In[6]:

# this calculates the norm (=absolute distance) between the centroids in each slice
distortionNorm = fun.coordDist(distortion)
print("distrotionNorm[:,:,{}]: {}".format(idxSlice, distortionNorm[idxSlice]))


# In[7]:

CT.getMask()
# creates CT.masked using CT.mask,
# but assigns each slice the centroid distance*1000*spacing as pixel value
CT.applyMask(replaceArray=distortionNorm, spacing=CT.img.GetSpacing()[0])

CT.showMasked()
# exports 3D image as .mha file
# fun.sitk_write(CT.masked, "../data/", "CT_distortionNorm.mha")

