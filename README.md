# mrt-for-rt

-------------------------

Using MRT instead of CT for radiotherapy planning would reduce dose administered to patient before actually stating the curative treatment. However, using a distorted image for therapy planing can have fatal consequences (especially for proton & ion beam therapy). Therefore making sure the obtained images depict the patients geometry accurately is crucial.

The aim of the project is to develop a tool that calculates the distortion of an MRT image by comparing it to an CT Image.

To better show the outcome of the developed tools "jupyter notebook" is used. For more information look here:
http://www.svds.com/jupyter-notebook-best-practices-for-data-science/

-------------------------

The ultimate goal would be a script that automatically finds all rods (of a phantom designed for this purpose) in a 3D CT & MRT image, compares them and creates a 3D vector field of the measured distortion not only of the position, but also the rod profile.

Steps:
0) way to import DICOM images and convert into useful format for further computing using python
		https://pyscience.wordpress.com/2014/10/19/image-segmentation-with-python-and-simpleitk/


1) compare CT & MRT image of 1 rod in 2D
a)	calculating the bending and displacement of the rod in MRT image relative to CT
	(important to keep in mind: CT shows plastic and filling, MRT only filling)

2) calculating the distortion of the circular rod profile
a)	calculating & plotting waveform along x-axis of image
	“waveform” = histogram showing signal intensity (ordinate)
			along image coordinate from left to right (abscissa)
b) 	measure along various angles (e.g. bootom left to upper right, etc → different angles)

3) imporving scirpts:
a)	2D → 3D
b)	automatically finding multiple rods
c)	working with the uncropped aligned CT & MRT images
d)	use python to align images, get inspired here:
		https://itk.org/Wiki/SimpleITK/Tutorials/MICCAI2015
		and https://github.com/InsightSoftwareConsortium/SimpleITKTutorialMICCAI2015
		or http://insightsoftwareconsortium.github.io/SimpleITK-Notebooks/ (notebooks 60-67)
		or https://aaltodoc.aalto.fi/handle/123456789/20891
		and https://aaltodoc.aalto.fi/bitstream/handle/123456789/20891/master_Tapani_Karoliina_2016.pdf?sequence=1&isAllowed=y
		or http://kevin-keraudren.blogspot.co.at/2014/12/medical-image-analysis-ipython-tutorials.html
		and https://github.com/curiale/Medical-Image-Analysis-IPython-Tutorials
