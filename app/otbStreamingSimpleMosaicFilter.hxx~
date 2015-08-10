/*=========================================================================
 *
 *  Copyright 2015 Irstea
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/


#ifndef __StreamingSimpleMosaicFilter_hxx
#define __StreamingSimpleMosaicFilter_hxx

#include "otbStreamingSimpleMosaicFilter.h"

namespace otb {

/**
 * Processing
 */
template <class TInputImage, class TOutputImage, class TInternalValueType>
void
StreamingSimpleMosaicFilter<TInputImage, TOutputImage, TInternalValueType>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, itk::ThreadIdType threadId)
{

  // Debug info
  itkDebugMacro(<<"Actually executing thread " << threadId << " in region " << outputRegionForThread);

  // Support progress methods/callbacks
  itk::ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels() );

  // Get output pointer
  OutputImageType * mosaicImage = this->GetOutput();

  // Get number of used inputs
  const unsigned int nbOfUsedInputImages = Superclass::GetNumberOfUsedInputImages();

  // Get number of bands
  const unsigned int nBands = Superclass::GetNumberOfBands();

  // Iterate through the thread region
  IteratorType outputIt(mosaicImage, outputRegionForThread);

  // Prepare interpolated pixel
  InternalPixelType interpolatedMathPixel;
  interpolatedMathPixel.SetSize(nBands);

  // Instanciate interpolators which are DEDICATED TO THE THREAD ! (so need to
  // copy the m_interpolator)
  InterpolatorPointerType interps[nbOfUsedInputImages];
  InputImageType *        currentImage[nbOfUsedInputImages];
  PrepareImageAccessors(currentImage, interps);

  // Container for geo coordinates
  OutputImagePointType geoPoint;

  for ( outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt )
    {
    // Prepare output pixel
    OutputImagePixelType outputPixel(Superclass::GetNoDataOutputPixel() );

    // Current pixel --> Geographical point
    mosaicImage->TransformIndexToPhysicalPoint (outputIt.GetIndex(), geoPoint) ;

    // Loop on used input images
    for (unsigned int i = 0 ; i < nbOfUsedInputImages ; i++)
      {
      // Get the input image pointer
      unsigned int imgIndex = Superclass::GetUsedInputImageIndice(i);

      // Check if the point is inside the transformed thread region
      // (i.e. the region in the current input image which match the thread
      // region)
      if (interps[i]->IsInsideBuffer(geoPoint) )
        {

        // Compute the interpolated pixel value
        InputImagePixelType interpolatedPixel = interps[i]->Evaluate(geoPoint);

        // Check that interpolated pixel is not empty
        if (IsPixelNotEmpty(interpolatedPixel) )
          {
          // Update the output pixel
          for (unsigned int band = 0 ; band < nBands ; band++)
            {
            if (this->GetShiftScaleInputImages() )
              {
              InternalValueType pixelValue = static_cast<InternalValueType>(interpolatedPixel[band]);
              this->ShiftScaleValue(pixelValue, imgIndex, band);
              outputPixel[band] = static_cast<OutputImageInternalPixelType>(pixelValue);
              }
            else
              {
              outputPixel[band] = static_cast<OutputImageInternalPixelType>(interpolatedPixel[band]);
              }
            }
          } // Interpolated pixel is not empty
        }   // point inside buffer
      }     // next image

    // Update output pixel value
    outputIt.Set(outputPixel);

    // Update progress
    progress.CompletedPixel();

    } // next output pixel

}

} // end namespace gtb

#endif
