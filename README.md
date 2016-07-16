# mrt-for-rt

-------------------------

Using MRT instead of CT for radiotherapy planning would reduce dose administered to patient before actually stating the curative treatment. However, using a distorted image for therapy planing can have fatal consequences (especially for proton & ion beam therapy). Therefore making sure the obtained images depict the patients geometry accurately is crucial.

The aim of the project is to develop a tool that calculates the distortion of an MRT image by comparing it to an CT Image.

-------------------------

The ultimate goal would be a script that automatically finds all rods (of a phantom designed for this purpose) in a 3D CT & MRT image, compares them and creates a 3D vector field of the measured distortion not only of the position, but also the rod profile.

Steps:
0) way to import DICOM images and convert into useful format for further computing using python

1) compare CT & MRT image of 1 rod in 2D
 a)	calculating the bending and displacement of the rod in MRT image relative to CT
	(important to keep in mind: CT shows plastic and filling, MRT only filling)
c)	2D → 3D
d)	automatically finding multiple rods
e)	working with the uncropped aligned CT & MRT images

2) calculating the distortion of the circular rod profile
 a)	calculating & plotting waveform along x-axis of image
	“waveform” = histogram showing signal intensity (ordinate)
			along image coordinate from left to right (abscissa)
b) 	measure along various angles (e.g. bootom left to upper right, etc → different angles)
c)	2D → 3D
d)	automatically finding multiple rods
e)	working with the uncropped aligned CT & MRT images

Helpful sites:
http://www.svds.com/jupyter-notebook-best-practices-for-data-science/
https://pyscience.wordpress.com/2014/10/19/image-segmentation-with-python-and-simpleitk/
