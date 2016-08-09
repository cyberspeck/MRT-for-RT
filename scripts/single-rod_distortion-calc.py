# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 18:44:50 2016
inspired by: https://itk.org/SimpleITKDoxygen/html/DicomSeriesReader_8py-example.html

@author: david
"""

import FunITK as fun
from FunITK import Volume
import sys

if len (sys.argv) < 4:
    print("\n Usage: single-rod_distortion-calc.py <input_directory_CT> <input_directory_MR> <output_file> \n \n")
    sys.exit(1)


CT = Volume(path=sys.argv[1], method="CT", seeds=[(6, 8, 10)])
CT.getCentroid()

MR = Volume(path=sys.argv[2], method="MR", seeds=[(6,8,10)])
MR.getCentroid()

# this calculates the coordinate difference of MR.centroid relative to CT.centroid
distortion = fun.coordShift(CT.centroid, MR.centroid)

# this calculates the norm (=absolute distance) between the centroids in each slice
distortionNorm = fun.coordDist(distortion)

# creates mask (pixel values either 0 or 1)
CT.getMask()
# creates CT.masked using CT.mask,
# but assigns each slice the centroid distance as pixel value
CT.applyMask(replaceArray=distortionNorm)

# exports 3D image as .mha file
fun.sitk_write(CT.masked, filename=sys.argv[3])
