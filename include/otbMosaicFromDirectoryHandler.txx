/*
 * otbMosaicFromDirectoryHandler.hxx
 *
 *  Created on: 24 mars 2016
 *      Author: cresson
 */

#ifndef MODULES_REMOTE_OTB_MosaicFromDirectoryHandler_INCLUDE_OTBMosaicFromDirectoryHandler_HXX_
#define MODULES_REMOTE_OTB_MosaicFromDirectoryHandler_INCLUDE_OTBMosaicFromDirectoryHandler_HXX_

#include "otbMosaicFromDirectoryHandler.h"
#include "otbImageFileWriter.h"

namespace otb
{

template <class TOutputImage>
MosaicFromDirectoryHandler<TOutputImage>
::MosaicFromDirectoryHandler()
 {
  mosaicFilter = MosaicFilterType::New();
  castFilter = CastFilterType::New();
 }

template <class TOutputImage>
MosaicFromDirectoryHandler<TOutputImage>
::~MosaicFromDirectoryHandler()
{
}

template <class TOutputImage>
void
MosaicFromDirectoryHandler<TOutputImage>
::SetDirectory(std::string directory)
 {
  if (directory[directory.size()-1] != '/')
    {
    // If not, we add the separator
    directory.append("/");
    }

  // Get the list of files in the directory
  itk::Directory::Pointer dir = itk::Directory::New();
  if (!dir->Load(directory.c_str()))
    {
    itkExceptionMacro(<< "Unable to browse directory " << directory);
    }

  // Instanciate a new mosaic filter
  mosaicFilter = MosaicFilterType::New();
  readers.clear();

  // Browse the directory
  for (unsigned int i = 0; i < dir->GetNumberOfFiles(); i++)
    {
    const char *filename = dir->GetFile(i);
    std::string sfilename(filename);
    sfilename = directory + sfilename;

    // Try to read the file
    otb::ImageIOBase::Pointer imageIO =
        otb::ImageIOFactory::CreateImageIO(sfilename.c_str(),otb::ImageIOFactory::ReadMode);
    if( imageIO.IsNotNull() )
      {
      ReaderPointerType reader = ReaderType::New();
      reader->SetFileName(sfilename);
      reader->UpdateOutputInformation();
      readers.push_back(reader);
      mosaicFilter->PushBackInput(reader->GetOutput());
      }
    else
      {
      //      itkWarningMacro(<<"Unable to read file " << sfilename);
      }

    }

 }

template <class TOutputImage>
void
MosaicFromDirectoryHandler<TOutputImage>
::GenerateOutputInformation()
 {
  mosaicFilter->SetOutputOrigin(m_OutputOrigin);
  mosaicFilter->SetOutputSpacing(m_OutputSpacing);
  mosaicFilter->SetOutputSize(m_OutputSize);
  mosaicFilter->SetAutomaticOutputParametersComputation(false);

  castFilter->SetInput(mosaicFilter->GetOutput());

  castFilter->GraftOutput( this->GetOutput() );
  castFilter->UpdateOutputInformation();
  this->GraftOutput( castFilter->GetOutput() );
 }

template <class TOutputImage>
void
MosaicFromDirectoryHandler<TOutputImage>
::GenerateData()
 {
  castFilter->GraftOutput( this->GetOutput() );
  castFilter->Update();
  this->GraftOutput( castFilter->GetOutput() );
 }

} // end namespace otb

#endif /* MODULES_REMOTE_OTB_MosaicFromDirectoryHandler_INCLUDE_OTBMosaicFromDirectoryHandler_HXX_ */
