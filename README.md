Summary
=======

This module provides one application dedicated to images mosaicking. 

Mosaic
----------------------------
This application computes mosaics of images. 

###Compositing###
######Methods######
* Simple compositing technique (-comp.feather none): Copies the last image over earlier ones in areas of overlap (Default)
* Feathering technique (-comp.feather large): Blends all images on the maximum overlapping areas (this produces seamless mosaics, but can cause blur effect where images are not perfectly aligned).
* Feathering technique (-comp.feather slim): Blends the last image over earlier ones in areas of overlap, on a given transition distance (seam can be visible from a certain zoom level, but it does not cause blur effect when images are not perfectly aligned)

######Performance tuning: Distance map images samplig ratio (-distancemap.sr)######
For performance issues, distance map images (used by feathering methods) size can be reduced. Indeed, distance map images are computed into a temporary file before the compositing process, because it couldn't be done in a streamable way. That is, the RAM must be enough to compute the entire distance map image from each single input image. Thus, the application provides an option to reduce the size of the distance map images. Setting -distancemap.sr 10 will induce a 0.1x wider distance map image, and will save 100x less RAM. To keep the original input size, just set -alphamasks.spacing 1.0 (Default distance map images size reduction ratio is 10.)

######Interpolator######
Supported interpolators are
* Nearest Neighborhood (-interpolator nn) (Default)
* Bicubic (-interpolator bco)
* Linear (-interpolator linear)

######Output spacing######
Spacing of the output mosaic can be changed using -output.spacing

######Cutline (vector data)######
Input vector data can additionaly feed the process for cutline selection (-vdcut)
One vector data is required for each input image, and must appear in the same order.
 
###Input images dynamic harmonization###
######Methods######
* None (-harmo.method none): No harmonization. Images pixels values are untouched.(Default)
* Band-by-band (-harmo.method none): Consists in minimizing a cost function based on images statistics in overlapping areas, for each band independently.
* True color harmonization (-harmo.method rgb): Consists in minimizing a cost function based on images statistics in overlapping areas, in a decorreleted color space suitable for true color processing (works only on true color images, i.e. RGB). Only the first 3 bands are processed by this method. More information about the rvb method can be found in the paper "_Natural Color Satellite Image Mosaicking Using Quadratic Programming in Decorrelated Color Space_" Cresson & Saint-Geours, July 2015, IEEE JSTARS Volume PP Issue 99

######Cost function######
Various cost function can be used for harmonization
* rmse (-harmo.cost rmse): Root mean squared error based cost function
* musig (-harmo.cost musig): Mean and Standard deviation based cost function
* mu (-harmo.cost mu): Mean based cost function

######Masks for statistics (vector data)######
Input vector data can additionaly feed the process for statistics computation, e.g. to mask clouds or water using -vdstats. One vector data is required for each input image, and must appear in the same order.

###Temporary files location (-tmpdir)###
Distance map images, and binary masks are temporary stored in the given temporary directory (Default is system directory).

Licence
=======

This code is provided under the CeCILL-B free software license agreement.

Contact
=======

For any issues regarding this module please contact Rémi Cresson.

remi.cresson@teledetection.fr

Irstea ((French) National Research Institute of Science and Technology for Environment and Agriculture)
www.irstea.fr
