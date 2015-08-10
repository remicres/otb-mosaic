Summary
=======

This module provides one application dedicated to images mosaicking. 

Mosaic
----------------------------
This application computes mosaics of images. 

*Input images can be mosaicked in 3 possible ways:
1.Using a simple compositing technique which copies the last image
over earlier ones in areas of overlap
2.Using a feathering technique, which blends all images on the maximum 
overlapping areas (this produces seamless mosaics, but can cause blur
effect where images are not perfectly aligned).
3.Using a feathering technique, which blends the last image over earlier 
ones in areas of overlap, on a given transition distance (seam can be visible
from a certain zoom level, but it does not cause blur effect when images are
not perfectly aligned)

For performance issues, alpha masks spacing (used by feathering techniques)
can be multiplied to speed up the process (alphamasks.spacing>>1), or have
a nicer blending (alphamasks.spacing equal to 1). Indeed, alpha masks are
computed before the compositing, using a non streamable filter.
That's why in certain cases (big images), we have to reduce the alpha masks
images size by multiplying its spacing.
Alpha masks are temporary stored in the given temporary directory (default is
system directory).

*Interpolator can be changed (Nearest Neighborhood, Linear, Bicubic)
*Output spacing can be chosen too.

*Input images dynamic can be managed using 3 possible methods:
1.Keep the original dynamic of images (fastest)
2.Band-by-band harmonization: minimize a cost function based on images 
statistics in overlapping areas, for each band independently.
3.True color harmonization: minimize a cost function based on images 
statistics in overlapping areas, in a decorreleted color space suitable
for true color processing (works only on true color images, i.e. RGB).
Only the first 3 bands are processed by this method.

More information about the rvb method can be found in the paper
"Natural Color Satellite Image Mosaicking Using Quadratic Programming
in Decorrelated Color Space" Cresson & Saint-Geours, July 2015, IEEE
JSTART Volume PP Issue 99

Remarks:
-For methods 2 and 3, images statistics has to be computed
-For methods 2 and 3, cost function can be chosen (rmse based, mean based,
mean and std based)

*Two kind of vector data can additionaly feed the process:
1.For statistics computation (e.g. to mask clouds or water)
2.For cutline selection



Licence
=======

This code is provided under the CeCILL-B free software license agreement.

Contact
=======

For any issues regarding this module please contact Rémi Cresson.

remi.cresson@teledetection.fr

Irstea ((French) National Research Institute of Science and Technology for Environment and Agriculture)
www.irstea.fr
