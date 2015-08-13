Summary
=======

This module provides one application dedicated to images mosaicking. 

Mosaic
----------------------------
This application computes mosaics of images. 

Compositing methods:
  1.Simple compositing technique (-comp.feather none) which copies the last image
over earlier ones in areas of overlap (Default)
  2.Feathering technique (-comp.feather large), which blends all images on the maximum 
overlapping areas (this produces seamless mosaics, but can cause blur
effect where images are not perfectly aligned).
  3.Feathering technique (-comp.feather slim), which blends the last image over earlier 
ones in areas of overlap, on a given transition distance (seam can be visible
from a certain zoom level, but it does not cause blur effect when images are
not perfectly aligned)

Feathering compositing methods performance tuning:
  For performance issues, distance map images (used by feathering methods)
size can be reduced. Indeed, distance map images are computed into a temporary
file before the compositing process, because it couldn't be done in a streamable
way. That is, the RAM must be enough to compute the entire distance map image
from every single input image. Thus, the application provides an option to
reduce the size of the distance map images, by multiplying the original spacing
of images. Setting -alphamasks.spacing 10 will induce a 0.1x wider distance
map image. To keep the original input size, just set -alphamasks.spacing 1.0.
Default distance map images size reduction ratio is 10.
Distance map images are temporary stored in the given temporary directory
(-tmpdir). Default is system directory.

Interpolator:
Supported interpolators:
  Nearest Neighborhood (Default)
  Bicubic
  Linear

Output spacing of the output mosaic can be changed using -output.spacing

Input images dynamic harmonization:
  None (Default)
  Band-by-band harmonization: minimize a cost function based on images 
statistics in overlapping areas, for each band independently.
  True color harmonization: minimize a cost function based on images 
statistics in overlapping areas, in a decorreleted color space suitable
for true color processing (works only on true color images, i.e. RGB).
Only the first 3 bands are processed by this method. More information 
about the rvb method can be found in the paper "Natural Color Satellite
Image Mosaicking Using Quadratic Programming in Decorrelated Color Space"
Cresson & Saint-Geours, July 2015, IEEE JSTARS Volume PP Issue 99

Input images dynamic harmonization cost functions:
  rmse: Root mean squared error based cost function
  musig: Mean and Standard deviation based cost function
  mu: Mean based cost function

Input vector data can additionaly feed the process (need one vector data 
for one input image, in the same order as appearing in input):
  For statistics computation, e.g. to mask clouds or water (-vdstats)
  For cutline selection (-vdcut)


Licence
=======

This code is provided under the CeCILL-B free software license agreement.

Contact
=======

For any issues regarding this module please contact Rémi Cresson.

remi.cresson@teledetection.fr

Irstea ((French) National Research Institute of Science and Technology for Environment and Agriculture)
www.irstea.fr
