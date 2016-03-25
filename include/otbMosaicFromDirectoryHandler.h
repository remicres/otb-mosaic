/*
 * otbMosaicFromDirectoryHandler.h
 *
 *  Created on: 24 mars 2016
 *      Author: cresson
 */

#ifndef MODULES_REMOTE_OTB_MosaicFromDirectoryHandler_INCLUDE_OTBMosaicFromDirectoryHandler_H_
#define MODULES_REMOTE_OTB_MosaicFromDirectoryHandler_INCLUDE_OTBMosaicFromDirectoryHandler_H_

#include "itkImageSource.h"
#include "itkExceptionObject.h"
#include "itkImageRegion.h"

#include "otbStreamingSimpleMosaicFilter.h"

#include "otbImageFileReader.h"
#include "itkDirectory.h"
#include "otbImageIOBase.h"
#include "otbImageIOFactory.h"

#include "otbMultiToMonoChannelExtractROI.h"

namespace otb
{
/** \class MosaicFromDirectoryHandler
 * \brief  Reads mask data.
 *
 * bla bla
 *
 *
 * \ingroup IOFilters
 *
 */
template <class TOutputImage>
class ITK_EXPORT MosaicFromDirectoryHandler : public itk::ImageSource<TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef MosaicFromDirectoryHandler                    Self;
  typedef itk::ImageSource<TOutputImage>     Superclass;
  typedef itk::SmartPointer<Self>            Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MosaicFromDirectoryHandler, ImageSource);

  /** The size of the output image. */
  typedef typename TOutputImage::SizeType  SizeType;

  /** The size of the output image. */
  typedef typename TOutputImage::IndexType  IndexType;

  /** The spacing of the output image. */
  typedef typename TOutputImage::SpacingType  SpacingType;

  /** The coordinate of the output image. */
  typedef typename TOutputImage::PointType  PointType;

  /** The region of the output image. */
  typedef typename TOutputImage::RegionType  ImageRegionType;

  /** The pixel type of the output image. */
  typedef typename TOutputImage::InternalPixelType OutputImagePixelType;

  /** Typedefs for mosaic filter */
  typedef otb::VectorImage<OutputImagePixelType> InternalMaskImageType;
  typedef otb::StreamingSimpleMosaicFilter<InternalMaskImageType> MosaicFilterType;
  typedef typename MosaicFilterType::Pointer MosaicFilterPointerType;

  /** Image reader */
  typedef otb::ImageFileReader<InternalMaskImageType> ReaderType;
  typedef typename ReaderType::Pointer ReaderPointerType;

  /** Cast filter */
  typedef otb::MultiToMonoChannelExtractROI<OutputImagePixelType, OutputImagePixelType> CastFilterType;
  typedef typename CastFilterType::Pointer CastFilterPointerType;

  /** Input directory accessors */
  itkGetMacro(Directory, std::string);
  void SetDirectory(std::string directory);

  /** Output parameters setters */
  itkSetMacro(OutputSpacing, SpacingType);
  itkSetMacro(OutputSize, SizeType);
  itkSetMacro(OutputOrigin, PointType);

  /** Prepare image allocation at the first call of the pipeline processing */
  virtual void GenerateOutputInformation(void);

  /** Does the real work. */
   virtual void GenerateData();

protected:
  MosaicFromDirectoryHandler();
  virtual ~MosaicFromDirectoryHandler();

  // Masks directory
  std::string m_Directory;

  // Output parameters
  SpacingType m_OutputSpacing;
  SizeType m_OutputSize;
  PointType m_OutputOrigin;

  // Internal filters
  MosaicFilterPointerType mosaicFilter;
  CastFilterPointerType castFilter;
  std::vector<ReaderPointerType> readers;


private:

  MosaicFromDirectoryHandler(const Self &); //purposely not implemented
  void operator =(const Self&); //purposely not implemented

};

} //namespace otb

#ifndef OTB_MANUAL_INSTANTIATION
#include <otbMosaicFromDirectoryHandler.txx>
#endif

#endif /* MODULES_REMOTE_OTB_MosaicFromDirectoryHandler_INCLUDE_OTBMosaicFromDirectoryHandler_H_ */
