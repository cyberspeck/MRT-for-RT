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

#seedPlastic = list(map(tuple, maxPlastic))
listPlastic = list(map(tuple, maxPlastic))

sitk_show(imgSmooth, interpolation="nearest")

# use Seeds as seedList
# seedPlastic = [(12,14), (15,15)]
seedFilling = [(14,14)]

seedPlastic = SimpleITK.VectorUIntList()
for pixel in listPlastic:
    print(type(a))
    print(pixel)
    seed = SimpleITK.VectorUInt32(pixel)
    print(type(seed))
    print(seed)
#    seedPlastic.push_back(seed)

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