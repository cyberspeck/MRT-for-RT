
# coding: utf-8

# In[2]:

import FunctionsCustom as fun
from FunctionsCustom import Volume
get_ipython().magic('pylab inline')

pathCT = "../data/cropped_CT-a/"
pathMR = "../data/cropped_MR-d-a/"
idxSlice = 10


# In[3]:

ct = Volume(path=pathCT, method="CT", ref=idxSlice)
mr = Volume(path=pathMR, method="MR", ref=idxSlice)


# In[4]:

ct.show()
mr.show()


# In[5]:

ct.centroid(show=idxSlice)
mr.centroid(show=idxSlice)


# In[6]:

distortion = fun.coordShift(ct.centroid, mr.centroid)
print(distortion)

