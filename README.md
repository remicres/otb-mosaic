Summary
=======

This module provides one application dedicated to images mosaicking. 

# Compositing

## Methods

* Simple compositing technique (__-comp.feather none__): Copies the last image over earlier ones in areas of overlap (Default)
* Feathering technique (__-comp.feather large__): Blends all images on the maximum overlapping areas (this produces seamless mosaics, but can cause blur effect where images are not perfectly aligned).
* Feathering technique (__-comp.feather slim__): Blends the last image over earlier ones in areas of overlap, on a given transition distance (seam can be visible from a certain zoom level, but it does not cause blur effect when images are not perfectly aligned)

## Performance tuning: Distance map images sampling ratio

For performance issues, distance map images (used by feathering methods) size can be reduced using __-distancemap.sr__. Indeed, distance map images are computed into a temporary file before the compositing process, because it couldn't be done in a streamable way. That is, the RAM must be enough to compute the entire distance map image from each single input image. Thus, the application provides an option to reduce the size of the distance map images. Setting -distancemap.sr 10 will induce a 0.1x wider distance map image, and will save 100x less RAM. To keep the original input size, just set -alphamasks.spacing 1.0 (Default distance map images size reduction ratio is 10.)

## Interpolator

Supported interpolators are
* Nearest Neighborhood (__-interpolator nn__) (Default)
* Bicubic (__-interpolator bco__)
* Linear (__-interpolator linear__)

## Output spacing

Spacing of the output mosaic can be changed using __-output.spacingx__ and __-output.spacingy__

## No-data value

The no-data value, used in both harmonization process and compositing process, can be changed using the __-nodata__ parameter

## Cutline (vector data)

Input vector data can additionaly feed the process for cutline selection using __-vdcut__. One vector data is required for each input image, and must appear in the same order.
 
# Radiometric/Colorimetric harmonization

More information about the method can be found in the paper "_Natural Color Satellite Image Mosaicking Using Quadratic Programming in Decorrelated Color Space_" Cresson & Saint-Geours, July 2015, IEEE JSTARS Volume 8 Issue 8 (https://doi.org/10.1109/JSTARS.2015.2449233)

## Methods

* None (__-harmo.method none__): No harmonization. Images pixels values are untouched.(Default)
* Band-by-band (__-harmo.method none__): Consists in minimizing a cost function based on images statistics in overlapping areas, for each band independently.
* True color harmonization (__-harmo.method rgb__): Consists in minimizing a cost function based on images statistics in overlapping areas, in a decorreleted color space suitable for true color processing (works only on true color images, i.e. RGB). Only the first 3 bands are processed by this method. 

## Cost function

Various cost function can be used for harmonization
* rmse (__-harmo.cost rmse__): Root mean squared error based cost function
* musig (__-harmo.cost musig__): Mean and Standard deviation based cost function
* mu (__-harmo.cost mu__): Mean based cost function

## Masks for statistics (vector data)

Input vector data can additionaly feed the process for statistics computation, e.g. to mask clouds or water, using __-vdstats__. One vector data is required for each input image, and must appear in the same order.

## Temporary files

Distance map images, and binary masks are temporary stored in the temporary directory specified by parameter __-tmpdir__ (Default is application directory).

How to use it?
=======

Mosaic can be used as any OTB application (gui, command line, python, c++, ...).

Licence
=======

This code is provided under the CeCILL-B free software license agreement.

Contact
=======

For any issues regarding this module please contact Rémi Cresson.

remi.cresson at irstea.fr

Irstea ((French) National Research Institute of Science and Technology for Environment and Agriculture)
www.irstea.fr
