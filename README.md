Summary
=======

This module provides one application dedicated to images mosaicking. 

Mosaic
----------------------------
This application computes mosaics of images. 

#Compositing#
##Methods##
###1.Simple compositing technique (-comp.feather none)###
Copies the last image over earlier ones in areas of overlap (Default)
###2.Feathering technique (-comp.feather large)###
Blends all images on the maximum overlapping areas (this produces
seamless mosaics, but can cause blur effect where images are not 
perfectly aligned).
###3.Feathering technique (-comp.feather slim)###
Blends the last image over earlier ones in areas of overlap, on a
given transition distance (seam can be visible from a certain zoom 
level, but it does not cause blur effect when images are not perfectly
aligned)

##Performance tuning##
###Distance map images samplig ratio (-distancemap.sr)###
For performance issues, distance map images (used by feathering methods)
size can be reduced. Indeed, distance map images are computed into a temporary
file before the compositing process, because it couldn't be done in a streamable
way. That is, the RAM must be enough to compute the entire distance map image
from each single input image. Thus, the application provides an option to
reduce the size of the distance map images. 
Setting -distancemap.sr 10 will induce a 0.1x wider distance map image, and
will save 100x less RAM. To keep the original input size, just set 
-alphamasks.spacing 1.0 (Default distance map images size reduction ratio is 10.)
###Distance map images temporary location (-tmpdir)###
Distance map images are temporary stored in the given temporary directory
(Default is system directory).

##Interpolator##
Supported interpolators are
###Nearest Neighborhood (-interpolator nn) (Default)###
###Bicubic (-interpolator bco)###
###Linear (-interpolator linear)###

##Output spacing##
It can be changed using -output.spacing

##Cutline##
Input vector data can additionaly feed the process for cutline selection (-vdcut)
One vector data is required for each input image, and must appear in the same order.
 
#Input images dynamic harmonization#
##Methods##
###None### (-harmo.method none)(Default)
No harmonization. Images pixels values are untouched.
###Band-by-band (-harmo.method none)###
Consists in minimizing a cost function based on images statistics in overlapping 
areas, for each band independently.
###True color harmonization (-harmo.method rgb)###
Consists in minimizing a cost function based on images statistics in overlapping
areas, in a decorreleted color space suitable for true color processing
(works only on true color images, i.e. RGB). Only the first 3 bands are processed
by this method. More information about the rvb method can be found in the paper
"Natural Color Satellite Image Mosaicking Using Quadratic Programming in 
Decorrelated Color Space" Cresson & Saint-Geours, July 2015, 
IEEE JSTARS Volume PP Issue 99

##Cost function##
Various cost function can be used for harmonization
###rmse (-harmo.cost rmse)###
Root mean squared error based cost function
###musig (-harmo.cost musig)###
Mean and Standard deviation based cost function
###mu (-harmo.cost mu)###
Mean based cost function

Input vector data can additionaly feed the process (need one vector data 
for one input image, in the same order as appearing in input):
  For statistics computation (-vdstats)
  For cutline selection (-vdcut)

##Statistics masks##
Input vector data can additionaly feed the process for statistics computation, 
e.g. to mask clouds or water using -vdstats.
One vector data is required for each input image, and must appear in the same order.

Licence
=======

This code is provided under the CeCILL-B free software license agreement.

Contact
=======

For any issues regarding this module please contact Rémi Cresson.

remi.cresson@teledetection.fr

Irstea ((French) National Research Institute of Science and Technology for Environment and Agriculture)
www.irstea.fr
