
# coding: utf-8

# In[1]:

import FunctionsCustom as fun
from FunctionsCustom import Volume
get_ipython().magic('pylab inline')


# In[2]:

pathCT = "../data/cropped_CT-a/"
pathMR = "../data/cropped_MR-d-a/"
idxSlice = 10

CT = Volume(path=pathCT, method="CT", ref=idxSlice, seeds=[(6, 8, idxSlice)])
CT.getCentroid()
CT.showCentroid()
CT.getMask()
CT.showMask()

rCT=CT.getDice(show=idxSlice)

MR = Volume(path=pathMR, method="MR", ref=idxSlice, seeds=[(6, 8, idxSlice)])
MR.getCentroid()
MR.showCentroid()
MR.getMask()
MR.showMask()

rMR=MR.getDice(show=idxSlice)


# In[3]:

CT.applyMask(replaceArray=CT.dice)
CT.showMasked()

