#include "itkFixedArray.h"

// Elevation handler
#include "otbWrapperElevationParametersHandler.h"
#include "otbWrapperApplicationFactory.h"

// Vector data rasterization
#include "otbVectorDataToLabelImageFilter.h"

// VectorData
#include "otbVectorDataIntoImageProjectionFilter.h"

// Application engine
#include "otbStandardFilterWatcher.h"
#include "itkFixedArray.h"

// Mosaic filters
#include "otbStreamingStatisticsMosaicFilter.h"
#include "otbStreamingLargeFeatherMosaicFilter.h"
#include "otbStreamingSimpleMosaicFilter.h"
#include "otbStreamingFeatherMosaicFilter.h"

// Resample filter
#include "otbStreamingResampleImageFilter.h"

// Pad image filter
#include "itkConstantPadImageFilter.h"

// otb ROI filter
#include "otbExtractROI.h"

// danielson distance image
#include "itkDanielssonDistanceMapImageFilter.h"

// Interpolators
#include "itkLinearInterpolateImageFunction.h"
#include "otbBCOInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

// Functors for lab <--> rgb color spaces
#include "otbUnaryFunctorImageFilter.h"

// Masks
#include "itkMaskImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "otbVectorImageToAmplitudeImageFilter.h"

// Solver for mosaic harmonization
#include "otbQuadraticallyConstrainedSimpleSolver.h"

// maths
#include "vcl_complex.h"

// temp filename
#include <boost/filesystem.hpp>

using namespace std;

/**
 * \class RGB2LABFunctor
 * \brief Base class for converting RGB into LAB color space (Ruderman et al.)
 */
template< class TInput, class TOutput>
class RGB2LAB
{
public:
  RGB2LAB() {
    M.set_size(3,3);
    M[0][0] = 0.3811; M[0][1] = 0.5783; M[0][2] = 0.0406;
    M[1][0] = 0.1967; M[1][1] = 0.7244; M[1][2] = 0.0790;
    M[2][0] = 0.0241; M[2][1] = 0.1288; M[2][2] = 0.8531;

    D1.set_size(3,3);
    D1.fill(0.0);
    D1[0][0] = 1.0 / vcl_sqrt(3.0);
    D1[1][1] = 1.0 / vcl_sqrt(6.0);
    D1[2][2] = 1.0 / vcl_sqrt(2.0);

    D2.set_size(3,3);
    D2.fill(1.0);
    D2[1][2] = -2.0;
    D2[2][1] = -1.0;
    D2[2][2] = 0.0;

  }

  ~RGB2LAB() {
  }

  bool operator!=( const RGB2LAB & ) const {
    return false;
  }

  bool operator==( const RGB2LAB & other ) const {
    return !(*this != other);
  }

  inline TOutput operator()( const TInput & A ) const
  {
    TOutput output;

    output.SetSize(3);
    if (A[0] == 0 && A[1] == 0 && A[2] == 0)
      {
      output.Fill(0);
      return output;
      }

    // RGB
    vnl_matrix<double> rgb(3,1);
    rgb[0][0] = A[0];
    rgb[1][0] = A[1];
    rgb[2][0] = A[2];

    // LMS
    vnl_matrix<double> lms(3,1);
    lms = M*rgb;

    // LMS (log10)
    const double log10 = vcl_log(10);
    lms[0][0] = vcl_log( lms[0][0] ) / log10;
    lms[1][0] = vcl_log( lms[1][0] ) / log10;
    lms[2][0] = vcl_log( lms[2][0] ) / log10;

    // LAB
    vnl_matrix<double> lab(3,1);
    lab = D1*(D2*lms);

    output[0] = lab[0][0];
    output[1] = lab[1][0];
    output[2] = lab[2][0];

    return output;
  }

  inline unsigned int GetOutputSize(){
    return 3;
  }

private:
  vnl_matrix<double> M;
  vnl_matrix<double> D1;
  vnl_matrix<double> D2;

};

/**
 * \class LAB2RGB Functor
 * \brief Base class for converting LAB into RGB color space (Ruderman et al.)
 */
template< class TInput, class TOutput>
class LAB2RGB
{
public:
  LAB2RGB() {
    M.set_size(3,3);
    M[0][0] = 4.4687; M[0][1] = -3.5887; M[0][2] = 0.1197;
    M[1][0] = -1.2197; M[1][1] = 2.3831; M[1][2] = -0.1626;
    M[2][0] = 0.0579; M[2][1] = -0.2584; M[2][2] = 1.1934;

    D1.set_size(3,3);
    D1.fill(0.0);
    D1[0][0] = 1.0 / vcl_sqrt(3.0);
    D1[1][1] = 1.0 / vcl_sqrt(6.0);
    D1[2][2] = 1.0 / vcl_sqrt(2.0);

    D2.set_size(3,3);
    D2.fill(1.0);
    D2[1][2] = -1.0;
    D2[2][1] = -2.0;
    D2[2][2] = 0.0;

  }

  ~LAB2RGB() {
  }

  bool operator!=( const LAB2RGB & ) const
  {
    return false;
  }

  bool operator==( const LAB2RGB & other ) const
  {
    return !(*this != other);
  }

  inline TOutput operator()( const TInput & A ) const
  {
    TOutput output;

    output.SetSize(3);

    if (A[0] == 0 && A[1] == 0 && A[2] == 0)
      {
      output.Fill(0);
      return output;
      }
    // LAB
    vnl_matrix<double> lab(3,1);
    lab[0][0] = A[0];
    lab[1][0] = A[1];
    lab[2][0] = A[2];

    // LMS
    vnl_matrix<double> lms(3,1);
    lms = D2*(D1*lab);
    lms[0][0] = vcl_pow(10.0, lms[0][0]);
    lms[1][0] = vcl_pow(10.0, lms[1][0]);
    lms[2][0] = vcl_pow(10.0, lms[2][0]);

    // RGB
    vnl_matrix<double> rgb(3,1);
    rgb = M*lms;

    output[0] = rgb[0][0];
    output[1] = rgb[1][0];
    output[2] = rgb[2][0];

    return output;
  }

  inline unsigned int GetOutputSize(){
    return 3;
  }

private:
  vnl_matrix<double> M;
  vnl_matrix<double> D1;
  vnl_matrix<double> D2;

};

namespace otb
{

enum
  {
  Interpolator_NNeighbor,
  Interpolator_BCO,
  Interpolator_Linear
  };

enum
  {
  Composition_Method_none,
  Composition_Method_large,
  Composition_Method_slim
  };

enum
  {
  Harmonisation_Method_none,
  Harmonisation_Method_bands,
  Harmonisation_Method_rgb
  };

namespace Wrapper
{

class Mosaic : public Application
{
public:
  /** Standard class typedefs. */
  typedef Mosaic                        Self;
  typedef Application                   Superclass;
  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Standard macro */
  itkNewMacro(Self);
  itkTypeMacro(Mosaic, Application);

  /** Solver (Cresson & St Geours) */
  typedef DoubleImageType::InternalPixelType                         SolverValueType;
  typedef otb::QuadraticallyConstrainedSimpleSolver<SolverValueType> SolverType;

  /* Mosaic Filters */
  typedef otb::StreamingSimpleMosaicFilter<FloatVectorImageType,
                                           FloatVectorImageType>      SimpleMosaicFilterType;
  typedef otb::StreamingLargeFeatherMosaicFilter<FloatVectorImageType, FloatVectorImageType,
                                                 DoubleImageType>     LargeFeatherMosaicFilterType;
  typedef otb::StreamingFeatherMosaicFilter<FloatVectorImageType, FloatVectorImageType,
                                                     DoubleImageType> SlimFeatherMosaicFilterType;
  typedef otb::StreamingStatisticsMosaicFilter<FloatVectorImageType,FloatVectorImageType,
                                               SolverValueType>       StatisticsMosaicFilterType;

  /* Binary masks */
  typedef otb::Image<bool>                                                                MaskImageType;
  typedef otb::Image<unsigned int>                                                        LabelImageType;
  typedef itk::MaskImageFilter<FloatVectorImageType, MaskImageType, FloatVectorImageType> MaskImageFilterType;
  typedef otb::ImageFileReader<MaskImageType>                                             MaskReaderType;

  /* UInt8 binary masks */
  typedef otb::Image<unsigned char>                                                    UInt8MaskImageType;
  typedef otb::ImageFileWriter<UInt8MaskImageType>                                     UInt8MaskWriterType;
  typedef otb::ImageFileReader<UInt8MaskImageType>                                     UInt8MaskReaderType;
  typedef itk::BinaryThresholdImageFilter<LabelImageType, UInt8MaskImageType>          LabelImageThresholdFilterType;
  typedef otb::VectorImageToAmplitudeImageFilter<FloatVectorImageType, FloatImageType> VectorImageToAmplitudeFilterType;
  typedef itk::BinaryThresholdImageFilter<FloatImageType, UInt8MaskImageType>          ImageThresholdFilterType;

  /* Alpha channel writer typedef */
  typedef otb::ImageFileReader<DoubleImageType> AlphaReaderType;

  /* Vector data filters typedefs */
  typedef otb::VectorDataIntoImageProjectionFilter<VectorDataType,
                                                   FloatVectorImageType> VectorDataReprojFilterType;
  typedef otb::VectorDataToLabelImageFilter<VectorDataType, LabelImageType> RasterizerType;

  /* Decorrelated color space <--> RGB functors typedefs */
  typedef LAB2RGB<FloatVectorImageType::PixelType, FloatVectorImageType::PixelType>               LAB2RGBFunctor;
  typedef RGB2LAB<FloatVectorImageType::PixelType, FloatVectorImageType::PixelType>               RGB2LABFunctor;
  typedef otb::UnaryFunctorImageFilter<FloatVectorImageType,FloatVectorImageType, RGB2LABFunctor> RGB2LABFilterType;
  typedef otb::UnaryFunctorImageFilter<FloatVectorImageType,FloatVectorImageType, LAB2RGBFunctor> LAB2RGBFilterType;

  /** Interpolators typedefs */
  typedef itk::LinearInterpolateImageFunction<FloatVectorImageType, double>          LinearInterpolationType;
  typedef itk::NearestNeighborInterpolateImageFunction<FloatVectorImageType, double> NearestNeighborInterpolationType;
  typedef otb::BCOInterpolateImageFunction<FloatVectorImageType>                     BCOInterpolationType;
  typedef itk::NearestNeighborInterpolateImageFunction<UInt8MaskImageType>           UInt8NNInterpolatorType;

  /** UInt8 filters typedefs */
  typedef otb::StreamingResampleImageFilter<UInt8MaskImageType, UInt8MaskImageType> UInt8ResampleImageFilterType;

private:

// Macro used to convert number to string
#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
    ( std::ostringstream() << std::dec << x ) ).str()

  /*
   * Create a reader array from a list of files
   */
  template<class TReaderType>
  vector<typename TReaderType::Pointer>
  CreateReaderArray(vector<string> inputFileNames, bool releaseDataFlag=true)
  {
    vector<typename TReaderType::Pointer> array;
    for (unsigned int i = 0 ; i < inputFileNames.size(); i++)
      {
      typename TReaderType::Pointer reader = TReaderType::New();
      reader->SetFileName(inputFileNames.at(i) );
      reader->UpdateOutputInformation();
      reader->SetReleaseDataFlag(releaseDataFlag);
      array.push_back(reader);
      }
    return array;
  }

  /*
   * Create an array of filters, which are all connected to the
   * input image list
   */
  template<class TOutputFilterType>
  vector<typename TOutputFilterType::Pointer>
  CreateConnectedFilterArrayToInputs(bool releaseDataFlag=true)
  {
    // Get the input image list
    FloatVectorImageListType::Pointer inputArray = this->GetParameterImageList("il");

    vector<typename TOutputFilterType::Pointer> outputArray;
    if (inputArray->Size() ==0)
      {
      otbAppLogFATAL("Filter array have wrong number of elements");
      }
    else
      {
      for (unsigned int i = 0 ; i < inputArray->Size() ; i++)
        {
        typename TOutputFilterType::Pointer filter = TOutputFilterType::New();
        filter->SetInput(inputArray->GetNthElement(i) );
        filter->SetReleaseDataFlag(releaseDataFlag);
        outputArray.push_back(filter);
        }
      }
    return outputArray;
  }

  /*
   * Create an array of filters, which are all connected to the given
   * input array of filters
   */
  template<class TFilterType1, class TOutputFilterType>
  vector<typename TOutputFilterType::Pointer>
  CreateConnectedFilterArray(
    vector<typename TFilterType1::Pointer>& inputArray,
    bool releaseDataFlag=true){
    vector<typename TOutputFilterType::Pointer> outputArray;
    if (inputArray.size() ==0)
      {
      otbAppLogFATAL("Filter array have wrong number of elements");
      }
    else
      {
      for (unsigned int i = 0 ; i < inputArray.size() ; i++)
        {
        typename TOutputFilterType::Pointer filter = TOutputFilterType::New();
        filter->SetInput(inputArray.at(i)->GetOutput() );
        filter->SetReleaseDataFlag(releaseDataFlag);
        outputArray.push_back(filter);
        }
      }
    return outputArray;
  }

  /*
   * Create an array of filters, which are all connected to the given
   * input array of filters (2 inputs)
   */
  template<class TFilterType1, class TFilterType2, class TOutputFilterType>
  vector<typename TOutputFilterType::Pointer>
  CreateConnectedFilterArray(
    vector<typename TFilterType1::Pointer>& inputArray1,
    vector<typename TFilterType2::Pointer>& inputArray2,
    bool releaseDataFlag=true){
    vector<typename TOutputFilterType::Pointer> outputArray;
    if (inputArray1.size() != inputArray2.size() || inputArray1.size()==0)
      {
      otbAppLogFATAL("Filter array have wrong number of elements");
      }
    else
      {
      for (unsigned int i = 0 ; i < inputArray1.size() ; i++)
        {
        typename TOutputFilterType::Pointer filter = TOutputFilterType::New();
        filter->SetInput1(inputArray1.at(i)->GetOutput() );
        filter->SetInput2(inputArray2.at(i)->GetOutput() );
        filter->SetReleaseDataFlag(releaseDataFlag);
        outputArray.push_back(filter);
        }
      }
    return outputArray;
  }

  /*
   * Create an array of filters, which are all connected to the given
   * input array of filters (2 inputs)
   */
  template<class TFilterType2, class TOutputFilterType>
  vector<typename TOutputFilterType::Pointer>
  CreateConnectedFilterArrayToInput(
    vector<typename TFilterType2::Pointer>& inputArray2,
    bool releaseDataFlag=true){

    // Get the input image list
    FloatVectorImageListType::Pointer inputArray = this->GetParameterImageList("il");

    vector<typename TOutputFilterType::Pointer> outputArray;
    if (inputArray->Size() != inputArray2.size() || inputArray2.size()==0)
      {
      otbAppLogFATAL("Filter array have wrong number of elements");
      }
    else
      {
      for (unsigned int i = 0 ; i < inputArray->Size() ; i++)
        {
        typename TOutputFilterType::Pointer filter = TOutputFilterType::New();
        filter->SetInput1(inputArray->GetNthElement(i) );
        filter->SetInput2(inputArray2.at(i)->GetOutput() );
        filter->SetReleaseDataFlag(releaseDataFlag);
        outputArray.push_back(filter);
        }
      }
    return outputArray;
  }

  /*
   * Create a mosaic filter, which is connected to the input array of filters
   */
  template<class TMosaicFilterType, class TFilterType>
  typename TMosaicFilterType::Pointer
  CreateConnectedMosaicFilter(vector<typename TFilterType::Pointer>& inputArray)
  {
    typename TMosaicFilterType::Pointer mosaicFilter = TMosaicFilterType::New();
    if (inputArray.size() ==0)
      {
      otbAppLogFATAL("Filter array have wrong number of elements");
      }
    else
      {
      for (unsigned int i = 0 ; i < inputArray.size() ; i++)
        {
        inputArray.at(i)->UpdateOutputInformation();
        mosaicFilter->PushBackInput(inputArray.at(i)->GetOutput() );
        }
      }
    return mosaicFilter;
  }

  /*
   * Create a mosaic filter, which is connected to the input arrays of filters
   */
  template<class TMosaicFilterType, class TFilterType1, class TFilterType2>
  typename TMosaicFilterType::Pointer
  CreateConnectedMosaicFilter(
    vector<typename TFilterType1::Pointer>& inputArray1,
    vector<typename TFilterType2::Pointer>& inputArray2)
  {
    typename TMosaicFilterType::Pointer mosaicFilter = TMosaicFilterType::New();
    if (inputArray1.size() != inputArray2.size() || inputArray1.size()==0)
      {
      otbAppLogFATAL("Filter array have wrong number of elements");
      }
    else
      {
      for (unsigned int i = 0 ; i < inputArray1.size() ; i++)
        {
        inputArray1.at(i)->UpdateOutputInformation();
        inputArray2.at(i)->UpdateOutputInformation();
        mosaicFilter->PushBackInputs(
          inputArray1.at(i)->GetOutput(),
          inputArray2.at(i)->GetOutput() );
        }
      }
    return mosaicFilter;
  }

  /*
   * Create a mosaic filter, which is connected to the inputs array
   */
  template<class TMosaicFilterType>
  typename TMosaicFilterType::Pointer
  CreateConnectedMosaicFilterToInputs()
  {
    // Get the input image list
    FloatVectorImageListType::Pointer inputArray = this->GetParameterImageList("il");

    typename TMosaicFilterType::Pointer mosaicFilter = TMosaicFilterType::New();
    if (inputArray->Size() ==0)
      {
      otbAppLogFATAL("ERROR (CreateConnectedMosaicFilterToInputs): Filter array have wrong number of elements");
      }
    else
      {
      for (unsigned int i = 0 ; i < inputArray->Size() ; i++)
        {
        mosaicFilter->PushBackInput(inputArray->GetNthElement(i) );
        }
      }
    return mosaicFilter;
  }

  /*
   * Create a mosaic filter, which is connected to the inputs array , and to
   * an array of filters
   */
  template<class TMosaicFilterType, class TFilterType>
  typename TMosaicFilterType::Pointer
  CreateConnectedMosaicFilterToInputs(
    vector<typename TFilterType::Pointer>& inputArray2)
  {
    // Get the input image list
    FloatVectorImageListType::Pointer inputArray = this->GetParameterImageList("il");

    typename TMosaicFilterType::Pointer mosaicFilter = TMosaicFilterType::New();
    if (inputArray->Size() != inputArray2.size() || inputArray2.size()==0)
      {
      otbAppLogFATAL("ERROR (CreateConnectedMosaicFilterToInputs): Filter array have wrong number of elements");
      }
    else
      {
      for (unsigned int i = 0 ; i < inputArray->Size() ; i++)
        {
        inputArray2.at(i)->UpdateOutputInformation();
        mosaicFilter->PushBackInputs(
          inputArray->GetNthElement(i),
          inputArray2.at(i)->GetOutput() );
        }
      }
    return mosaicFilter;
  }

  /*
   * Delete a temporary file
   */
  void deleteFile(string filename)
  {
    if( remove(filename.c_str() ) != 0 )
      {
      otbAppLogWARNING( "Error deleting file " << filename );
      }
    else
      {
      otbAppLogDEBUG( "File " << filename << " successfully deleted" );
      }
  }

  /*
   * Delete a list of temporary files
   */
  void deleteFiles(vector<string> fileList)
  {
    for (unsigned int i = 0 ; i < fileList.size() ; i++)
      {
      deleteFile(fileList.at(i) );
      }
  }

  /*
   * This function generates a filename that looks like:
   * <m_TempDirectory>/<tag>_<id>.tif
   */
  string GenerateFileName(string tag, int id)
  {
    // Create a filename
    string outputFile;

    if (m_TempDirectory.size() != 0)
      {
      outputFile += m_TempDirectory + "/";
      }
    outputFile += tag + "_" + SSTR(id) + ".tif";
    return outputFile;
  }

  void DoInit()
  {
    SetName("Mosaic");
    SetDescription("Perform mosaicking of input images");

    // Documentation
    SetDocName("Mosaic");
    SetDocLongDescription("This application performs mosaicking of images");
    SetDocLimitations("None");
    SetDocAuthors("Remi Cresson");
    SetDocSeeAlso(" ");

    AddDocTag(Tags::Manip);

    // Input image
    AddParameter(ParameterType_InputImageList,  "il",   "Input Images");
    SetParameterDescription("il", "Input images to mosaic");

    // Input vector data (cutline)
    AddParameter(ParameterType_InputVectorDataList, "vdcut", "Input VectorDatas for composition");
    SetParameterDescription("vdcut", "VectorData files to be used for cut input images");
    MandatoryOff("vdcut");

    // Input vector data (statistics masks)
    AddParameter(ParameterType_InputVectorDataList, "vdstats", "Input VectorDatas for statistics");
    SetParameterDescription("vdstats", "VectorData files to be used for statistics computation (harmonization)");
    MandatoryOff("vdstats");

    // comp (compositing)
    AddParameter(ParameterType_Group,"comp","Mosaic compositing mode");
    SetParameterDescription("comp","This group of parameters allow to set the mosaic compositing");

    // comp.feathering
    AddParameter(ParameterType_Choice,"comp.feather","Feathering method");
    SetParameterDescription("comp.feather","Set the feathering method for composition");

    // comp.feather.none
    AddChoice("comp.feather.none","The simplest composition mode");
    SetParameterDescription("comp.feather.none",
                            "No feathering method is used (Very fast). Images are stacked in overlaps");

    // comp.feather.large
    AddChoice("comp.feather.large","The large blending composition mode");
    SetParameterDescription("comp.feather.large",
                            "Blends all images on the maximum overlapping areas. May generate blur when inputs are not well aligned");

    // comp.feather.slim
    AddChoice("comp.feather.slim","The fine blending composition mode");
    SetParameterDescription("comp.feather.slim",
                            "Blends the last image over earlier ones in areas of overlap, on a given transition distance");
    // comp.feather.slim.exponent (i.e. blending smoothness)
    AddParameter(ParameterType_Float, "comp.feather.slim.exponent", "Transition smoothness (Unitary exponent = linear transition)");
    SetDefaultParameterFloat("comp.feather.slim.exponent", 1.0);
    SetMinimumParameterFloatValue("comp.feather.slim.exponent", 0);
    MandatoryOff("comp.feather.slim.exponent");
    // comp.feather.slim.lenght (i.e. blending lenght)
    AddParameter(ParameterType_Float, "comp.feather.slim.lenght", "Transition lenght (In cartographic units)");
    MandatoryOn("comp.feather.slim.lenght");
    SetMinimumParameterFloatValue("comp.feather.slim.lenght", 0);
    MandatoryOff("comp.feather.slim.lenght");

    // harmo (harmonization)
    AddParameter(ParameterType_Group,"harmo","Spectral bands harmonization mode");
    SetParameterDescription("harmo","This group of parameters allow to tune the mosaic harmonization behavior");

    // harmo.method
    AddParameter(ParameterType_Choice,"harmo.method","harmonization method");
    SetParameterDescription("harmo.method","Set the harmonization method");
    AddChoice("harmo.method.none","None");
    SetParameterDescription("harmo.method.none","No automatic harmonization is done");
    AddChoice("harmo.method.band","Spectral bands");
    SetParameterDescription("harmo.method.band",
                            "Automatic harmonization is performed on each individual spectral band");
    AddChoice("harmo.method.rgb","True colors");
    SetParameterDescription("harmo.method.rgb",
                            "Works on RGB colors only. Use an harmonization method based on quadratic programming in decorrelated color space");

    // harmo.cost (harmonization cost function)
    AddParameter(ParameterType_Choice,"harmo.cost","harmonization cost function");
    SetParameterDescription("harmo.cost","Set the harmonization cost function");
    AddChoice("harmo.cost.rmse","RMSE based");
    AddChoice("harmo.cost.musig","Mean+Std based");
    AddChoice("harmo.cost.mu","Mean based");

    // Output image
    AddParameter(ParameterType_OutputImage,  "out",   "Output image");
    SetParameterDescription("out"," Output image.");

    // Interpolators
    AddParameter(ParameterType_Choice,   "interpolator", "Interpolation");
    SetParameterDescription("interpolator",
                            "This group of parameters allows to define how the input image will be interpolated during resampling.");
    MandatoryOff("interpolator");
    AddChoice("interpolator.nn",     "Nearest Neighbor interpolation");
    SetParameterDescription("interpolator.nn",
                            "Nearest neighbor interpolation leads to poor image quality, but it is very fast.");
    AddChoice("interpolator.bco",    "Bicubic interpolation");
    AddParameter(ParameterType_Radius, "interpolator.bco.radius", "Radius for bicubic interpolation");
    SetParameterDescription("interpolator.bco.radius",
                            "This parameter allows to control the size of the bicubic interpolation filter. If the target pixel size is higher than the input pixel size, increasing this parameter will reduce aliasing artefacts.");
    AddChoice("interpolator.linear", "Linear interpolation");
    SetParameterDescription("interpolator.linear",
                            "Linear interpolation leads to average image quality but is quite fast");
    SetDefaultParameterInt("interpolator.bco.radius", 2);

    // Spacing of the output image
    AddParameter(ParameterType_Group, "output", "Output Image Grid");
    SetParameterDescription("output",
                            "This group of parameters allows to define the grid on which the input image will be resampled.");
    AddParameter(ParameterType_Float, "output.spacing", "Pixel Size");
    SetParameterDescription("output.spacing",
                            "Size of each pixel along X axis (meters for cartographic projections, degrees for geographic ones)");
    MandatoryOff("output.spacing");

    // temp. dir
    AddParameter(ParameterType_String, "tmpdir", "Temporary directory storing masks");
    SetParameterDescription("tmpdir","Temporary directory storing masks");
    MandatoryOff("tmpdir");

    AddParameter(ParameterType_Group, "alphamasks", "Alpha masks computation");
    AddParameter(ParameterType_Float, "alphamasks.spacing", "Alpha masks spacing, in multiple of the corresponding input image spacing");
    SetParameterDescription("alphamasks.spacing",
                            "Can be increased if input images are too big, or in order to speed up the process");
    SetDefaultParameterFloat("alphamasks.spacing", 10);

    AddRAMParameter();

  }

  void DoUpdateParameters()
  {
    // Nothing to do here : all parameters are independent
  }

  /*
   * Write a binary mask to disk from a vector data
   */
  void RasterizeBinaryMask(VectorDataType * vd, FloatVectorImageType * reference,
                           string outputFileName, double spacingRatio, bool invert=false){

    // Reproject VectorData
    VectorDataReprojFilterType::Pointer vdReproj = VectorDataReprojFilterType::New();

    vdReproj->SetInputVectorData(vd);
    vdReproj->SetInputImage(reference);
    vdReproj->Update();

    // Rasterize vector data
    GDALAllRegister();
    RasterizerType::Pointer rasterizer = RasterizerType::New();
    rasterizer->AddVectorData(vdReproj->GetOutput() );
    rasterizer->SetOutputOrigin(reference->GetOrigin() );
    LabelImageType::SizeType outputSize = reference->GetLargestPossibleRegion().GetSize();
    outputSize[0] /= spacingRatio;
    outputSize[1] /= spacingRatio;
    rasterizer->SetOutputSize(outputSize);
    LabelImageType::SpacingType outputSpacing = reference->GetSpacing();
    outputSpacing[0] *= spacingRatio;
    outputSpacing[1] *= spacingRatio;
    rasterizer->SetOutputSpacing(outputSpacing);
    rasterizer->SetOutputProjectionRef(reference->GetProjectionRef() );
    rasterizer->SetBurnAttribute("________");
    rasterizer->SetGlobalWarningDisplay(false);

    // Threshold background
    LabelImageThresholdFilterType::Pointer labelThreshold = LabelImageThresholdFilterType::New();
    labelThreshold->SetInput(rasterizer->GetOutput() );
    if (invert)
      {
      labelThreshold->SetInsideValue(itk::NumericTraits<UInt8MaskImageType::InternalPixelType>::max() );
      labelThreshold->SetOutsideValue(itk::NumericTraits<UInt8MaskImageType::InternalPixelType>::Zero);
      }
    else
      {
      labelThreshold->SetInsideValue(itk::NumericTraits<UInt8MaskImageType::InternalPixelType>::Zero);
      labelThreshold->SetOutsideValue(itk::NumericTraits<UInt8MaskImageType::InternalPixelType>::max() );
      }
    labelThreshold->SetLowerThreshold(1);
    labelThreshold->SetUpperThreshold(itk::NumericTraits<LabelImageType::InternalPixelType>::max() );

    // Write file
    UInt8MaskWriterType::Pointer writer = UInt8MaskWriterType::New();
    writer->SetInput(labelThreshold->GetOutput() );
    writer->SetFileName(outputFileName);
    AddProcess(writer,"Writing binary mask (from vector data) "+outputFileName);
    writer->Update();

  }

  /*
   * Write on disk the distance image
   * Inputs:
   * -input filename (binary mask)
   * -output filename (distance image)
   */
  void WriteDistanceImage(string inputBinaryMaskFileName, string outputDistanceImageFileName)
  {
    /** Typedefs */
    typedef itk::DanielssonDistanceMapImageFilter<UInt8MaskImageType,
                                                  DoubleImageType,
                                                  DoubleImageType> ApproximateSignedDistanceMapImageFilterType;
    typedef otb::ImageFileWriter<DoubleImageType>                                         WriterType;
    typedef itk::ConstantPadImageFilter<UInt8MaskImageType, UInt8MaskImageType>           PadFilterType;
    typedef otb::ExtractROI<UInt8MaskImageType::PixelType, UInt8MaskImageType::PixelType> ExtractFilterType;

    // Read the binary mask image
    UInt8MaskReaderType::Pointer reader = UInt8MaskReaderType::New();
    reader->SetFileName(inputBinaryMaskFileName);

    // Pad the image
    const unsigned int paddingRadius = 2;
    PadFilterType::Pointer padFilter = PadFilterType::New();
    padFilter->SetConstant(itk::NumericTraits<UInt8MaskImageType::InternalPixelType>::max() );
    padFilter->SetInput(reader->GetOutput() );
    UInt8MaskImageType::SizeType padLowerBound; padLowerBound.Fill(paddingRadius);
    UInt8MaskImageType::SizeType padUpperBound; padUpperBound.Fill(paddingRadius);
    padFilter->SetPadLowerBound(padLowerBound);
    padFilter->SetPadUpperBound(padUpperBound);
    padFilter->UpdateOutputInformation();

    // Extract region
    ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
    extractFilter->SetInput(padFilter->GetOutput() );
    extractFilter->SetExtractionRegion(padFilter->GetOutput()->GetLargestPossibleRegion() );

    // Compute the approximate signed distance image
    ApproximateSignedDistanceMapImageFilterType::Pointer approximateSignedDistanceMapImageFilter =
      ApproximateSignedDistanceMapImageFilterType::New();
    approximateSignedDistanceMapImageFilter->SetInput(extractFilter->GetOutput() );
    approximateSignedDistanceMapImageFilter->SetInputIsBinary(true);
    approximateSignedDistanceMapImageFilter->SetUseImageSpacing(true);
    approximateSignedDistanceMapImageFilter->Update();

    // Write the distance image
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outputDistanceImageFileName);
    writer->SetInput(approximateSignedDistanceMapImageFilter->GetOutput() );
    AddProcess(writer,"Writing distance image "+outputDistanceImageFileName);
    writer->Update();

  }

  /*
   * Write a binary mask from an input image
   */
  void WriteBinaryMask(FloatVectorImageType* referenceImage,
                       string outputFileName, double spacingRatio = 1.0)
  {
    // Vector image to amplitude image
    VectorImageToAmplitudeFilterType::Pointer ampFilter = VectorImageToAmplitudeFilterType::New();

    ampFilter->SetInput(referenceImage);

    // Threshold image
    ImageThresholdFilterType::Pointer thresholdFilter = ImageThresholdFilterType::New();
    thresholdFilter->SetInput(ampFilter->GetOutput() );
    thresholdFilter->SetOutsideValue(itk::NumericTraits<UInt8MaskImageType::InternalPixelType>::Zero);
    thresholdFilter->SetInsideValue(itk::NumericTraits<UInt8MaskImageType::InternalPixelType>::max() );
    thresholdFilter->SetLowerThreshold(itk::NumericTraits<FloatVectorImageType::InternalPixelType>::Zero);
    thresholdFilter->SetUpperThreshold(itk::NumericTraits<FloatVectorImageType::InternalPixelType>::Zero);
    thresholdFilter->UpdateOutputInformation();

    // Resample image
    UInt8ResampleImageFilterType::Pointer resampler = UInt8ResampleImageFilterType::New();
    resampler->SetInput(thresholdFilter->GetOutput() );
    LabelImageType::SizeType outputSize = thresholdFilter->GetOutput()->GetLargestPossibleRegion().GetSize();
    outputSize[0] = outputSize[0] / spacingRatio + 1;
    outputSize[1] = outputSize[1] / spacingRatio + 1;
    resampler->SetOutputSize(outputSize);
    LabelImageType::SpacingType outputSpacing = thresholdFilter->GetOutput()->GetSpacing();
    outputSpacing[0] *= spacingRatio;
    outputSpacing[1] *= spacingRatio;
    resampler->SetOutputSpacing(outputSpacing);
    resampler->SetOutputOrigin(thresholdFilter->GetOutput()->GetOrigin() );

    // Write image
    UInt8MaskWriterType::Pointer writer = UInt8MaskWriterType::New();
    writer->SetInput(resampler->GetOutput() );
    writer->SetFileName(outputFileName);
    AddProcess(writer,"Writing binary mask (from image boundaries) "+outputFileName);
    writer->Update();
  }

  /*
   * Write the distance image of the input image #id
   */
  void WriteDistanceImageFromCutline(FloatVectorImageType* image, VectorDataType * vd, string outputFileName)
  {

    // Generate a temporary filenames for the resampled mask
    string temporaryFileName = GenerateFileName("tmp_binary_rasterized_mask", 0);

    // Write a binary mask
    RasterizeBinaryMask(vd, image, temporaryFileName, m_AlphaMasksSpacingMultiplicator);

    // Create distance image
    WriteDistanceImage(temporaryFileName, outputFileName);

    // Delete the temporary file
    deleteFile(temporaryFileName);

  }

  /*
   * Write the distance image of the input image #id
   */
  string WriteDistanceImageFromBoundaries(FloatVectorImageType* image, string outputFileName)
  {

    // Generate a temporary filenames for the resampled mask
    string temporaryFileName = GenerateFileName("tmp_binary_mask", 0);

    // Write a temporary binary mask
    WriteBinaryMask(image, temporaryFileName, m_AlphaMasksSpacingMultiplicator);

    // Create distance image
    WriteDistanceImage(temporaryFileName, outputFileName);

    // Delete the temporary binary mask file
    deleteFile(temporaryFileName);

    // Prepare image for new data
    image->PrepareForNewData();

    return outputFileName;

  }

  /*
   * Set the correction model to the mosaic filter
   */
  template <class TMosaicFilterType>
  void SetCorrectionModel(typename TMosaicFilterType::Pointer& filter)
  {
    // If no harmonization is required, we return
    if (this->GetParameterInt("harmo.method")==Harmonisation_Method_none)
      {
      return;
      }

    // Set the solver objective function (aka "correction cost function")
    SolverType::Pointer solver = SolverType::New();
    switch ( GetParameterInt("harmo.cost") )
      {
      case SolverType::Cost_Function_mu:
        solver->SetMeanBased();
        otbAppLogINFO("Correction cost function: mean based cost function");
        break;
      case SolverType::Cost_Function_musig:
        solver->SetMeanAndStandardDeviationBased();
        otbAppLogINFO("Correction cost function: mean-stdev based cost function");
        break;
      case SolverType::Cost_Function_rmse:
        solver->SetRMSEBased();
        otbAppLogINFO("Correction cost function: RMSE based cost function");
        break;
      default:
        otbWarningMacro("Unknow correction cost function. Setting to RMSE");
        solver->SetRMSEBased();
        break;
      }

    vnl_matrix<long> area = m_StatsFilter->GetAreaInPixels();

    // Colorimetric correction (Cresson & St Geours)
    filter->UpdateOutputInformation();
    const unsigned int          nImages = this->GetParameterImageList("il")->Size();
    const unsigned int          nBands = filter->GetOutput()->GetNumberOfComponentsPerPixel();
    vnl_matrix<SolverValueType> scales(nImages, nBands, 1);
    solver->SetAreaInOverlaps(m_StatsFilter->GetAreaInPixels() );

    for (unsigned int band = 0 ; band < nBands ; band++)
      {
      otbAppLogINFO("computing correction model for band " << band );
      solver->SetMeanInOverlaps(m_StatsFilter->GetMean(band) );
      solver->SetStandardDeviationInOverlaps(m_StatsFilter->GetStDev(band) );
      solver->SetMeanOfProductsInOverlaps(m_StatsFilter->GetProdMean(band) );
      solver->Solve();

      // Keep shifts & scales
      scales.set_column(band, solver->GetOutputCorrectionModel() );

      otbAppLogINFO("\n\t[ Band " << band << " ]"
                                  << "\n\tGains  : " << solver->GetOutputCorrectionModel() );
      }

    // Set filter correction model
    filter->SetScaleMatrix(scales);
    filter->SetShiftScaleInputImages(true);
  }

  /*
   * Set interpolator to the mosaic filter
   */
  template <class TMosaicFilterType>
  void SetInterpolator(typename TMosaicFilterType::Pointer& filter)
  {
    // Get Interpolator
    switch ( GetParameterInt("interpolator") )
      {
      case Interpolator_Linear:
        {
        typedef itk::LinearInterpolateImageFunction<FloatVectorImageType,
                                                    double>          LinearInterpolationType;
        LinearInterpolationType::Pointer interpolator = LinearInterpolationType::New();
        filter->SetInterpolator(interpolator);
        }
        break;
      case Interpolator_NNeighbor:
        {
        typedef itk::NearestNeighborInterpolateImageFunction<FloatVectorImageType,
                                                             double> NearestNeighborInterpolationType;
        NearestNeighborInterpolationType::Pointer interpolator = NearestNeighborInterpolationType::New();
        filter->SetInterpolator(interpolator);
        }
        break;
      case Interpolator_BCO:
        {
        typedef otb::BCOInterpolateImageFunction<FloatVectorImageType> BCOInterpolationType;
        BCOInterpolationType::Pointer interpolator = BCOInterpolationType::New();
        interpolator->SetRadius(GetParameterInt("interpolator.bco.radius") );
        filter->SetInterpolator(interpolator);
        }
        break;
      }
  }

  /*
   * Change a mosaic filter spacing
   */
  template <class TMosaicFilterType>
  void SetSpacing(typename TMosaicFilterType::Pointer& filter)
  {
    FloatVectorImageType::SpacingType outputSpacing;

    filter->SetAutomaticOutputParametersComputation(true);
    filter->UpdateOutputInformation();
    if (this->HasValue("output.spacing") )
      {
      outputSpacing[0] = GetParameterFloat("output.spacing");
      outputSpacing[1] = -1.0*GetParameterFloat("output.spacing");
      typename TMosaicFilterType::OutputImageSpacingType spacing = filter->GetOutputSpacing();
      typename TMosaicFilterType::OutputImageSizeType size = filter->GetOutputSize();
      typename TMosaicFilterType::OutputImagePointType origin = filter->GetOutputOrigin();
      size[0]*=(spacing[0]/outputSpacing[0]);
      size[1]*=(spacing[1]/outputSpacing[1]);
      filter->SetOutputSpacing(outputSpacing);
      filter->SetOutputOrigin(origin);
      filter->SetOutputSize(size);
      filter->SetAutomaticOutputParametersComputation(false);
      filter->UpdateOutputInformation();
      }
  }

  /*
   * Configure the mosaic filter:
   * -set interpolator
   * -set output spacing
   * -set correction model
   * -...
   */
  template <class TMosaicFilterType>
  void ConfigureMosaicFilter(typename TMosaicFilterType::Pointer& filter)
  {
    SetInterpolator<TMosaicFilterType>(filter);
    SetCorrectionModel<TMosaicFilterType>(filter);
    SetSpacing<TMosaicFilterType>(filter);
  }

  template <class TMosaicFilterType>
  void ComputeDistanceOffset(typename TMosaicFilterType::Pointer& filter)
  {
	  filter->UpdateOutputInformation();
	  typename TMosaicFilterType::OutputImageSpacingType spacing = filter->GetOutputSpacing();
	  float multiplicator = GetParameterFloat("alphamasks.spacing");
	  float maxSpacing = vnl_math_max(multiplicator*vnl_math_abs(spacing[0]),
			  multiplicator*vnl_math_abs(spacing[1]));
	  filter->SetDistanceOffset(maxSpacing);
  }

  void DoExecute()
  {
    GDALAllRegister();

    // Get the input image list
    FloatVectorImageListType::Pointer imagesList = this->GetParameterImageList("il");

    // Get the input VectorData list
    VectorDataListType* statsVectorDataList = GetParameterVectorDataList("vdstats");
    VectorDataListType* cutVectorDataList = GetParameterVectorDataList("vdcut");

    // Check the number of images & masks
    const unsigned int nImages = imagesList->Size();
    const unsigned int nMasks = statsVectorDataList->Size();
    const unsigned int nCutline = cutVectorDataList->Size();
    if ( nCutline != 0 && nCutline != nImages)
      {
      otbAppLogFATAL("Number of input cutlines (" << nCutline
                                                  << ") should be equal to number of images (" << nImages << ")");
      }
    if ( nMasks != 0 && nMasks != nImages)
      {
      otbAppLogFATAL("Number of input masks (" << nMasks
                                               << ") should be equal to number of images (" << nImages << ")");
      }

    // Get parameters
    m_TempDirectory = GetParameterAsString("tmpdir");
    if (m_TempDirectory.empty())
      {
      boost::filesystem::path temp = boost::filesystem::temp_directory_path();
      m_TempDirectory = temp.native();
      }
    otbAppLogINFO(<<"Temporary directory is:"<<m_TempDirectory);
    m_AlphaMasksSpacingMultiplicator = GetParameterFloat("alphamasks.spacing");

    /////////////////////////////////////////////////////////////
    //				 Compute stats (if needed)
    /////////////////////////////////////////////////////////////

    if (this->GetParameterInt("harmo.method")==Harmonisation_Method_none)
      {
      otbAppLogINFO("No harmonisation method is selected");
      }
    else
      {
      if (nMasks != 0)
        {

        // Write binary masks used for statistics computation
        otbAppLogINFO("Computing masks for statistics... ");
        masksForStatsFileNameList.clear();
        for (unsigned int i = 0 ; i < nMasks ; i++)
          {
          string outputFileName = GenerateFileName("tmp_binary_mask_for_stats", i);
          RasterizeBinaryMask(statsVectorDataList->GetNthElement(i),
                              imagesList->GetNthElement(i), outputFileName, 1.0, true);
          masksForStatsFileNameList.push_back( outputFileName );
          }

        // Create stats masks reader array
        m_MaskReaderForStats = CreateReaderArray<MaskReaderType>(masksForStatsFileNameList);

        if (this->GetParameterInt("harmo.method")==Harmonisation_Method_rgb)
          {
          otbAppLogINFO("Computing statistics in a decorrelated colors space suitable for true-colors ");

          // rgb2lab array
          m_rgb2labFilter = CreateConnectedFilterArrayToInputs<RGB2LABFilterType>();

          // maskImageFilter array
          m_MaskImageFilter =  CreateConnectedFilterArray<RGB2LABFilterType,MaskReaderType, MaskImageFilterType>(
              m_rgb2labFilter, m_MaskReaderForStats);

          }
        else if (this->GetParameterInt("harmo.method")==Harmonisation_Method_bands)
          {
          otbAppLogINFO("Computing statistics in the radiometric color space");

          // maskImageFilter array
          m_MaskImageFilter =  CreateConnectedFilterArrayToInput<MaskReaderType, MaskImageFilterType>(
              m_MaskReaderForStats);
          }
        else
          {
          itkExceptionMacro("Unknown harmo.method parameter");
          }

        // maskImageFilter-->statsFilter
        m_StatsFilter = CreateConnectedMosaicFilter<StatisticsMosaicFilterType,MaskImageFilterType>(m_MaskImageFilter);
        }
      else // no input mask
        {
        otbAppLogINFO("No input shp for stats. Skipping masks generation.");

        if (this->GetParameterInt("harmo.method")==Harmonisation_Method_rgb)
          {
          otbAppLogINFO("Computing statistics in a decorrelated colors space suitable for true-colors ");

          // rgb2lab array
          m_rgb2labFilter = CreateConnectedFilterArrayToInputs<RGB2LABFilterType>();

          // rgb2lab-->statsFilter
          m_StatsFilter = CreateConnectedMosaicFilter<StatisticsMosaicFilterType,RGB2LABFilterType>(m_rgb2labFilter);
          }
        else
          {
          otbAppLogINFO("Computing statistics in the radiometric color space");

          // inputs-->statsFilter
          m_StatsFilter = CreateConnectedMosaicFilterToInputs<StatisticsMosaicFilterType>();
          }
        }

      // Compute stats (explicit streaming)
      m_StatsFilter->UpdateOutputInformation();
      typedef otb::RAMDrivenAdaptativeStreamingManager<FloatVectorImageType> StreamingManagerType;
      StreamingManagerType::Pointer m_StreamingManager = StreamingManagerType::New();
      m_StreamingManager->PrepareStreaming(m_StatsFilter->GetOutput(),
                                           m_StatsFilter->GetOutput()->GetLargestPossibleRegion() );
      int m_NumberOfDivisions = m_StreamingManager->GetNumberOfSplits();
      for (int m_CurrentDivision = 0;
           m_CurrentDivision < m_NumberOfDivisions;
           m_CurrentDivision++)
        {
        FloatVectorImageType::RegionType streamRegion = m_StreamingManager->GetSplit(m_CurrentDivision);
        m_StatsFilter->GetOutput()->SetRequestedRegion(streamRegion);
        m_StatsFilter->GetOutput()->PropagateRequestedRegion();
        m_StatsFilter->GetOutput()->UpdateOutputData();
        AddProcess(m_StatsFilter, "Computing stats (Tile "
                   + SSTR(m_CurrentDivision) + " on "
                   + SSTR(m_NumberOfDivisions) + ")");
        }
      }

    /////////////////////////////////////////////////////////////
    //				Instanciate the mosaic filters
    /////////////////////////////////////////////////////////////

    if (GetParameterInt("comp.feather")==Composition_Method_none)
      {
      // No need for alpha masks
      otbAppLogINFO("Composition method is set to none. Skipping alpha masks. (No cutline will be used)");

      // Use a simple mosaic filter
      if (this->GetParameterInt("harmo.method")==Harmonisation_Method_rgb)
        {
        otbAppLogINFO("Performing simple composition method in rgb color space");

        // Filter input: rgb2lab filter array
        m_simpleMosaicFilter = CreateConnectedMosaicFilter<SimpleMosaicFilterType,
                                                           RGB2LABFilterType>(m_rgb2labFilter);

        // Filter output: lab2rgb filter array
        m_lab2rgbFilter = LAB2RGBFilterType::New();
        m_lab2rgbFilter->SetInput(m_simpleMosaicFilter->GetOutput() );

        // lab2rgb filter array output: Outputs
        SetParameterOutputImage("out", m_lab2rgbFilter->GetOutput() );
        }
      else // radiometric color space
        {
        otbAppLogINFO("Performing simple composition method in radiometric color space");

        // Filter input: Inputs
        m_simpleMosaicFilter = CreateConnectedMosaicFilterToInputs<SimpleMosaicFilterType>();

        // Filter output: Outputs
        SetParameterOutputImage("out", m_simpleMosaicFilter->GetOutput() );
        }
      ConfigureMosaicFilter<SimpleMosaicFilterType>(m_simpleMosaicFilter);
      }
    else
      {
      // Compute distance images
      otbAppLogINFO("Computing alpha masks... ");
      distanceImageFileNameList.clear();
      for (unsigned int i = 0 ; i < nImages ; i++)
        {
        string outputFileName = GenerateFileName("tmp_distance_image", i);
        if (nCutline != 0)
          {
          WriteDistanceImageFromCutline(imagesList->GetNthElement(i),
                                        cutVectorDataList->GetNthElement(i), outputFileName);
          }
        else // do not use vector cutline (use boundaries)
          {
          WriteDistanceImageFromBoundaries(imagesList->GetNthElement(i), outputFileName);
          }
        distanceImageFileNameList.push_back(outputFileName);
        }

      // Create alpha masks reader array
      m_AlphaReader = CreateReaderArray<AlphaReaderType>(distanceImageFileNameList);

      // Instanciate the mosaic filter depending the harmonization mode chosen
      if (this->GetParameterInt("harmo.method")==Harmonisation_Method_rgb)
        {
        otbAppLogINFO("Mosaic compositing in decorrelated color space");

        // Last filter is the lab-->rgb functor
        m_lab2rgbFilter = LAB2RGBFilterType::New();
        SetParameterOutputImage("out", m_lab2rgbFilter->GetOutput() );

        if (GetParameterInt("comp.feather")==Composition_Method_large)
          {
          m_largeFeatherMosaicFilter = CreateConnectedMosaicFilter<LargeFeatherMosaicFilterType,
                                                                   RGB2LABFilterType, AlphaReaderType>(m_rgb2labFilter,
                                                                                                       m_AlphaReader);
          m_lab2rgbFilter->SetInput(m_largeFeatherMosaicFilter->GetOutput() );
          }
        else if (GetParameterInt("comp.feather")==Composition_Method_slim)
          {
          m_smallFeatherMosaicFilter = CreateConnectedMosaicFilter<SlimFeatherMosaicFilterType,
                                                                RGB2LABFilterType, AlphaReaderType>(m_rgb2labFilter,
                                                                                                    m_AlphaReader);
          m_lab2rgbFilter->SetInput(m_smallFeatherMosaicFilter->GetOutput() );
          }
        }
      else // radiometric color space
        {
        otbAppLogINFO("Mosaic compositing in radiometric color space");

        if (GetParameterInt("comp.feather")==Composition_Method_large)
          {
          m_largeFeatherMosaicFilter = CreateConnectedMosaicFilterToInputs<LargeFeatherMosaicFilterType,
                                                                           AlphaReaderType>(m_AlphaReader);
          SetParameterOutputImage("out", m_largeFeatherMosaicFilter->GetOutput() );
          }
        else if (GetParameterInt("comp.feather")==Composition_Method_slim)
          {
          m_smallFeatherMosaicFilter = CreateConnectedMosaicFilterToInputs<SlimFeatherMosaicFilterType,
                                                                        AlphaReaderType>(m_AlphaReader);
          SetParameterOutputImage("out", m_smallFeatherMosaicFilter->GetOutput() );
          }
        }

      if (GetParameterInt("comp.feather")==Composition_Method_large)
        {
        ConfigureMosaicFilter<LargeFeatherMosaicFilterType>(m_largeFeatherMosaicFilter);
        ComputeDistanceOffset<LargeFeatherMosaicFilterType>(m_largeFeatherMosaicFilter);
        }
      else if (GetParameterInt("comp.feather")==Composition_Method_slim)
        {
        ConfigureMosaicFilter<SlimFeatherMosaicFilterType>(m_smallFeatherMosaicFilter);
        ComputeDistanceOffset<SlimFeatherMosaicFilterType>(m_smallFeatherMosaicFilter);

        // Set transition lenght and smoothness
        m_smallFeatherMosaicFilter->SetFeatheringTransitionDistance(GetParameterFloat("comp.feather.slim.lenght"));
        m_smallFeatherMosaicFilter->SetFeatheringSmoothness(GetParameterFloat("comp.feather.slim.exponent"));
        }
      }

  }   // DOExecute()

  void AfterExecuteAndWriteOutputs()
  {
    otbAppLogINFO("Clean temporary files");
    deleteFiles(distanceImageFileNameList);
    deleteFiles(masksForStatsFileNameList);
    otbAppLogINFO("Done");
  }

  // Mosaicking filters
  SimpleMosaicFilterType::Pointer m_simpleMosaicFilter;
  LargeFeatherMosaicFilterType::Pointer m_largeFeatherMosaicFilter;
  SlimFeatherMosaicFilterType::Pointer m_smallFeatherMosaicFilter;
  StatisticsMosaicFilterType::Pointer m_StatsFilter;

  // rgb<-->lab functors
  vector<RGB2LABFilterType::Pointer> m_rgb2labFilter;
  LAB2RGBFilterType::Pointer m_lab2rgbFilter;

  // mask image filters
  vector<MaskImageFilterType::Pointer> m_MaskImageFilter;
  vector<MaskReaderType::Pointer> m_MaskReaderForStats;
  vector<AlphaReaderType::Pointer> m_AlphaReader;

  // Parameters
  string m_TempDirectory; // Temp. directory
  vector<string> distanceImageFileNameList;
  vector<string> masksForStatsFileNameList;

  double m_AlphaMasksSpacingMultiplicator; // Alpha masks spacing

};
}
}

OTB_APPLICATION_EXPORT( otb::Wrapper::Mosaic )
