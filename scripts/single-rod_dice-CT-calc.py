# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 18:44:50 2016
inspired by:
https://itk.org/SimpleITKDoxygen/html/DicomSeriesReader_8py-example.html

@author: david
"""

import FunITK as fun
from FunITK import Volume
import sys
import numpy as np

if len(sys.argv) < 3:
    print("\n Usage: single-rod_dice-CT-calc.py <input_directory_CT-DICOM> <output_file> \n \n")
    sys.exit(1)


vol = Volume(path=sys.argv[1], method="CT", seeds=[(6, 8, 10)])
vol.getCentroid()
vol.getMask()

vol.getDice()
vol.applyMask(replaceArray=vol.dice)

# exports 3D image as .mha file
fun.sitk_write(vol.masked, filename=sys.argv[2])
