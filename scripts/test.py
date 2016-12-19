# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 19:14:12 2016

@author: david
"""


import FunITK as fun
from FunITK import Volume
import datetime

#CT images not usable 173-192
pathCT = "../data_final/CT_x1/"
pathMR = "../data_final/MR_x1/"
pathCT_x4 = "../data_final/CT_x4/"
pathMR_x4 = "../data_final/MR_x4/"
pathCT_x9 = "../data_final/CT_x9/"
pathMR_x9 = "../data_final/MR_x9/"
pathCT_x16 = "../data_final/CT_x16/"
pathMR_x16 = "../data_final/MR_x16/"
pathCT_x25 = "../data_final/CT_x25/"
pathMR_x25 = "../data_final/MR_x25/"
pathCT_x100 = "../data_final/CT_x100/"
pathMR_x100 = "../data_final/MR_x100/"
idxSlice = 10

CT = Volume(path=pathCT, method="CT", info="x1", ref=idxSlice, seeds=[(4, 20, idxSlice)])
#a = CT.getCentroid(percentLimit='auto', iterations=10)
#CT.showCentroid()
b = CT.getCentroid(threshold='auto')
CT.showCentroid()
CT.getMask()
CT.getDice(b, CT.mask)

CT_x4 = Volume(path=pathCT_x4, method="CT", info="x4", ref=idxSlice, seeds=[(8, 40, idxSlice)])
#c = CT_x4.getCentroid(percentLimit='auto', iterations=10)
#CT_x4.showCentroid()
d = CT_x4.getCentroid(threshold='auto')
CT_x4.showCentroid()
CT_x4.getMask()
CT_x4.getDice(d, CT_x4.mask)

CT_x9 = Volume(path=pathCT_x9, method="CT", info="x9", ref=idxSlice, seeds=[(12, 60, idxSlice)])
#e = CT_x9.getCentroid(percentLimit='auto', iterations=5)
#CT_x9.showCentroid()
f = CT_x9.getCentroid(threshold='auto')
CT_x9.showCentroid()
CT_x9.getMask()
CT_x9.getDice(f, CT_x9.mask)

CT_x16 = Volume(path=pathCT_x16, method="CT", info="x16", ref=idxSlice, seeds=[(16, 80, idxSlice)])
#e = CT_x16.getCentroid(percentLimit='auto', iterations=5)
#CT_x16.showCentroid()
f = CT_x16.getCentroid(threshold='auto')
CT_x16.showCentroid()
CT_x16.getMask()
CT_x16.getDice(f, CT_x16.mask)

CT_x25 = Volume(path=pathCT_x25, method="CT", info="x25", ref=idxSlice, seeds=[(20, 100, idxSlice)])
#CT_x25.getCentroid(percentLimit='auto', iterations=5)
#CT_x25.showCentroid()
CT_x25.getCentroid(threshold='auto')
CT_x25.showCentroid()
CT_x25.getMask()
CT_x25.getDice()

CT_x100 = Volume(path=pathCT_x100, method="CT", info="x100", ref=idxSlice, seeds=[(40, 200, idxSlice)])
#CT_x100.getCentroid(percentLimit='auto', iterations=5)
#CT_x100.showCentroid()
CT_x100.getCentroid(threshold='auto')
CT_x100.showCentroid()
CT_x100.getMask()
CT_x100.getDice()



MR = Volume(path=pathMR, method="MR", info="x1", ref=idxSlice, seeds=[(7, 21, idxSlice)])
#MR.showSeed(pixel=True)
#g = MR.getCentroid(percentLimit='auto')
#MR.showCentroid()
h = MR.getCentroid(threshold='auto')
MR.showCentroid()
MR.getMask()
MR.getDice(h, MR.mask)

MR_x4 = Volume(path=pathMR_x4, method="MR", info="x4", ref=idxSlice, seeds=[(14, 42, idxSlice)])
#MR_x4.showSeed(pixel=True)
#i = MR_x4.getCentroid(percentLimit='auto', iterations=5)
#MR_x4.showCentroid()
j = MR_x4.getCentroid(threshold='auto')
MR_x4.showCentroid()
MR_x4.getMask()
MR_x4.getDice(j, MR_x4.mask)

MR_x9 = Volume(path=pathMR_x9, method="MR", info="x9", ref=idxSlice, seeds=[(21, 63, idxSlice)])
#k = MR_x9.getCentroid(percentLimit='auto', iterations=10, halfShift=0.1)
#MR_x9.showCentroid()
l = MR_x9.getCentroid(threshold='auto')
MR_x9.showCentroid()
MR_x9.getMask()
MR_x9.getDice(l, MR_x9.mask)

MR_x16 = Volume(path=pathMR_x16, method="MR", info="x16", ref=idxSlice, seeds=[(28, 84, idxSlice)])
#k = MR_x16.getCentroid(percentLimit='auto', iterations=10, halfShift=0.1)
#MR_x16.showCentroid()
l = MR_x16.getCentroid(threshold='auto')
MR_x16.showCentroid()
MR_x16.getMask()
MR_x16.getDice(l, MR_x16.mask)

MR_x25 = Volume(path=pathMR_x25, method="MR", info="x25", ref=idxSlice, seeds=[(35, 105, idxSlice)])
#MR_x25.getCentroid(percentLimit='auto', iterations=10, halfShift=0.1)
#MR_x25.showCentroid()
MR_x25.getCentroid(threshold='auto')
MR_x25.showCentroid()
MR_x25.getMask()
MR_x25.getDice()

MR_x100 = Volume(path=pathMR_x100, method="MR", info="x100", ref=idxSlice, seeds=[(70, 210, idxSlice)])
#MR_x100.getCentroid(percentLimit='auto', iterations=10, halfShift=0.1)
#MR_x100.showCentroid()
MR_x100.getCentroid(threshold='auto')
MR_x100.showCentroid()
MR_x100.getMask()
MR_x100.getDice()