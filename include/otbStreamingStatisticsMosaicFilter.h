#ifndef __StreamingStatisticsMosaicFilter_H
#define __StreamingStatisticsMosaicFilter_H

#include "otbStreamingMosaicFilterBase.h"
#include "otbPersistentMosaicFilter.h"
#include "otbPersistentFilterStreamingDecorator.h"
#include "itkArray.h"
#include "itkSimpleDataObjectDecorator.h"
#include "itkImageRegionConstIteratorWithOnlyIndex.h"

namespace otb
{
/** \class PersistentStatisticsMosaicFilter
 * \brief Computes the statistics of a mosaic from an input images set.
 * The output pixel value is equal to the number of overlaps
 *
 * Support streaming
 *
 * The pixels must support the operator ==, +, /, etc.
 * The "no data value", output spacing, interpolator can be chosen.
 * The behavior of the filter is to compute input images statistics
 * as they were in a layered fashion: interpolators are used to
 * compute pixels values of all input images for a given output pixel.
 * This is used to compute the following matrices: mean, standard
 * deviation, min, max, and means of pixels products. Let's denote
 * X one of these matrices, then X\{ij} is the statistic of the image
 * i in the overlapping area with the image j.
 *
 * \ingroup OTBMosaic
 */
template <class TInputImage, class TOutputImage=TInputImage, class TInternalValueType=double>
class ITK_EXPORT PersistentStatisticsMosaicFilter : public otb::PersistentMosaicFilter<TInputImage, TOutputImage,
                                                                                         TInternalValueType>
{
public:

  /** Standard Self typedef */
  typedef PersistentStatisticsMosaicFilter                                               Self;
  typedef PersistentMosaicFilter<TInputImage, TOutputImage, TInternalValueType> Superclass;
  typedef itk::SmartPointer<Self>                                                       Pointer;
  typedef itk::SmartPointer<const Self>                                                 ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(PersistentStatisticsMosaicFilter, PersistentMosaicFilter);

  /** Input image typedefs.  */
  typedef typename Superclass::InputImageType              InputImageType;
  typedef typename Superclass::InputImagePointer           InputImagePointer;
  typedef typename Superclass::InputImagePointType         InputImagePointType;
  typedef typename Superclass::InputImagePixelType         InputImagePixelType;
  typedef typename Superclass::InputImageIndexType         InputImageIndexType;
  typedef typename Superclass::InputImageSizeType          InputImageSizeType;
  typedef typename Superclass::InputImageSpacingType       InputImageSpacingType;
  typedef typename Superclass::InputImageInternalPixelType InputImageInternalPixelType;
  typedef typename Superclass::InputImageRegionType        InputImageRegionType;

  /** Output image typedefs.  */
  typedef typename Superclass::OutputImageType              OutputImageType;
  typedef typename Superclass::OutputImagePointer           OutputImagePointer;
  typedef typename Superclass::OutputImagePointType         OutputImagePointType;
  typedef typename Superclass::OutputImagePixelType         OutputImagePixelType;
  typedef typename Superclass::OutputImageIndexType         OutputImageIndexType;
  typedef typename Superclass::OutputImageSizeType          OutputImageSizeType;
  typedef typename Superclass::OutputImageSpacingType       OutputImageSpacingType;
  typedef typename Superclass::OutputImageInternalPixelType OutputImageInternalPixelType;
  typedef typename Superclass::OutputImageRegionType        OutputImageRegionType;

  /** Internal computing typedef support. */
  typedef typename Superclass::InternalValueType       InternalValueType;
  typedef typename Superclass::ContinuousIndexType     ContinuousIndexType;
  typedef typename Superclass::InterpolatorType        InterpolatorType;
  typedef typename Superclass::InterpolatorPointerType InterpolatorPointerType;
  typedef typename Superclass::DefaultInterpolatorType DefaultInterpolatorType;
  typedef typename Superclass::InternalImageType       InternalImageType;
  typedef typename Superclass::InternalPixelType       InternalPixelType;
  typedef typename Superclass::ConstIteratorType       ConstIteratorType;
  typedef typename Superclass::StreamingTraitsType     StreamingTraitsType;

  typedef itk::ImageRegionConstIteratorWithOnlyIndex<OutputImageType> IteratorType;

  /** Typedefs for statistics */
  typedef vnl_matrix<InternalValueType>       RealMatrixType;
  typedef vnl_vector<InternalValueType>       RealVectorType;
  typedef vnl_matrix<long>                    LongMatrixType;
  typedef vnl_vector<long>                    LongVectorType;
  typedef std::vector<RealMatrixType>         RealMatrixListType;
  typedef itk::SimpleDataObjectDecorator<RealMatrixListType>  RealMatrixListObjectType;
  typedef itk::SimpleDataObjectDecorator<LongMatrixType>      LongMatrixObjectType;

  /** Smart Pointer type to a DataObject. */
  typedef typename itk::DataObject::Pointer DataObjectPointer;
  typedef itk::ProcessObject::DataObjectPointerArraySizeType DataObjectPointerArraySizeType;

  /** Overrided methods */
  virtual void AllocateOutputs();
  virtual void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, itk::ThreadIdType threadId );
  virtual void Reset();
  virtual void Synthetize();

  /** Make a DataObject of the correct type to be used as the specified
   * output. */
  DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx) ITK_OVERRIDE;
  using Superclass::MakeOutput;

  /** Return the computed Mean. */
  RealMatrixListType GetMeans() const
  {
    return this->GetMeansOutput()->Get();
  }
  RealMatrixListObjectType* GetMeansOutput();
  const RealMatrixListObjectType* GetMeansOutput() const;

  /** Return the computed Std. */
  RealMatrixListType GetStds() const
  {
    return this->GetStdsOutput()->Get();
  }
  RealMatrixListObjectType* GetStdsOutput();
  const RealMatrixListObjectType* GetStdsOutput() const;

  /** Return the computed Min. */
  RealMatrixListType GetMins() const
  {
    return this->GetMinsOutput()->Get();
  }
  RealMatrixListObjectType* GetMinsOutput();
  const RealMatrixListObjectType* GetMinsOutput() const;

  /** Return the computed Max. */
  RealMatrixListType GetMaxs() const
  {
    return this->GetMaxsOutput()->Get();
  }
  RealMatrixListObjectType* GetMaxsOutput();
  const RealMatrixListObjectType* GetMaxsOutput() const;

  /** Return the computed MeansOfProduct. */
  RealMatrixListType GetMeansOfProducts() const
  {
    return this->GetMeansOfProductsOutput()->Get();
  }
  RealMatrixListObjectType* GetMeansOfProductsOutput();
  const RealMatrixListObjectType* GetMeansOfProductsOutput() const;

  /** Return the computed areas. */
  LongMatrixType GetAreas() const
  {
    return this->GetAreasOutput()->Get();
  }
  LongMatrixObjectType* GetAreasOutput();
  const LongMatrixObjectType* GetAreasOutput() const;

protected:
  PersistentStatisticsMosaicFilter();
  virtual ~PersistentStatisticsMosaicFilter() {}

  /** PersistentStatisticsMosaicFilter can be implemented as a multithreaded filter.
   * Therefore, this implementation provides a ThreadedGenerateData() routine
   * which is called for each processing thread. The output image data is
   * allocated automatically by the superclass prior to calling
   * ThreadedGenerateData().  ThreadedGenerateData can only write to the
   * portion of the output image specified by the parameter
   * "outputRegionForThread"
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData()  */



  /** Class for storing thread results:
   * -sum of values
   * -sum of squared values
   * -min value
   * -max value
   * -count
   */
  class ThreadResultsContainer {
public:
    /** Default constructor */
    ThreadResultsContainer(){
    }

    /* Constructor with size */
    ThreadResultsContainer(unsigned int nbOfBands, unsigned int nbOfSamples)
    {
      Clear(nbOfBands, nbOfSamples);
    }

    /* Copy constructor */
    ThreadResultsContainer(const ThreadResultsContainer& other)
    {
      m_count   = LongVectorType(other.m_count);
      m_sum     = RealMatrixType(other.m_sum);
      m_cosum   = RealMatrixType(other.m_cosum);
      m_sqSum   = RealMatrixType(other.m_sqSum);
      m_min     = RealMatrixType(other.m_min);
      m_max     = RealMatrixType(other.m_max);
    }

    /* Clear routine: Resize at the specified dimension and clear values */
    void Clear(unsigned int nbOfBands, unsigned int nbOfSamples)
    {
      InternalValueType zeroValue  = itk::NumericTraits<InternalValueType>::Zero;
      InternalValueType supValue   = itk::NumericTraits<InternalValueType>::max();
      InternalValueType infValue   = itk::NumericTraits<InternalValueType>::NonpositiveMin();

      m_count   = LongVectorType(nbOfSamples,0);
      m_sum     = RealMatrixType(nbOfBands,nbOfSamples,zeroValue);
      m_cosum   = RealMatrixType(nbOfBands,nbOfSamples,zeroValue);
      m_sqSum   = RealMatrixType(nbOfBands,nbOfSamples,zeroValue);
      m_min     = RealMatrixType(nbOfBands,nbOfSamples,supValue);
      m_max     = RealMatrixType(nbOfBands,nbOfSamples,infValue);
    }

    /* Pixel update */
    void Update( const InputImagePixelType& pixel, unsigned int sampleId)
    {
      unsigned int nbOfBands = pixel.Size();

      m_count[sampleId]++;
      for (unsigned int band = 0 ; band < nbOfBands ; band++)
        {
        // Cast
        InternalValueType pixelValue = static_cast<InternalValueType>(pixel[band]);

        // Update Min & max
        if (pixelValue < m_min[band][sampleId])
          m_min[band][sampleId] = pixelValue;
        if (pixelValue > m_max[band][sampleId])
          m_max[band][sampleId] = pixelValue;

        // Update Sums
        m_sum[band][sampleId] += pixelValue;
        m_sqSum[band][sampleId] += pixelValue*pixelValue;
        }
    }

    /* 2 pixels update */
    void Update( const InputImagePixelType& pixel_i,const InputImagePixelType& pixel_j, unsigned int sampleId)
    {
      Update(pixel_i, sampleId);
      unsigned int nbOfBands = pixel_i.Size();
      for (unsigned int band = 0 ; band < nbOfBands ; band++)
        {
        // Cast
        InternalValueType pixelValue_i = static_cast<InternalValueType>(pixel_i[band]);
        InternalValueType pixelValue_j = static_cast<InternalValueType>(pixel_j[band]);

        m_cosum[band][sampleId] += pixelValue_i*pixelValue_j;
        }
    }

    /* Self update */
    void Update(const ThreadResultsContainer& other)
    {
      unsigned int nbOfBands = other.m_sum.rows();
      unsigned int nbOfSamples = other.m_sum.cols();

      for (unsigned int sampleId = 0 ; sampleId < nbOfSamples ; sampleId++)
        {
        m_count[sampleId] += other.m_count[sampleId];
        for (unsigned int band = 0 ; band < nbOfBands ; band++)
          {
          m_sum[band][sampleId] += other.m_sum[band][sampleId];
          m_cosum[band][sampleId] += other.m_cosum[band][sampleId];
          m_sqSum[band][sampleId] += other.m_sqSum[band][sampleId];
          if (other.m_min[band][sampleId] < m_min[band][sampleId])
            m_min[band][sampleId] = other.m_min[band][sampleId];
          if (other.m_max[band][sampleId] > m_max[band][sampleId])
            m_max[band][sampleId] = other.m_max[band][sampleId];
          }
        }
    }

    RealMatrixType   m_sum;
    RealMatrixType   m_sqSum;
    RealMatrixType   m_cosum;
    RealMatrixType   m_min;
    RealMatrixType   m_max;
    LongVectorType   m_count;
  };

  // Internal threads count
  std::vector<ThreadResultsContainer> internalThreadResults;

  // Results
  RealMatrixListType m_Means;
  RealMatrixListType m_Stds;
  RealMatrixListType m_ProdMeans;
  RealMatrixListType m_Mins;
  RealMatrixListType m_Maxs;
  LongMatrixType     m_Area;

private:
  PersistentStatisticsMosaicFilter(const Self&); //purposely not implemented
  void operator=(const Self&);                  //purposely not implemented

}; // end of class



/*
 * Decorator
 */
template<class TInputImage,class TOutputImage,class TInternalValueType>
class ITK_EXPORT StreamingStatisticsMosaicFilter :
  public PersistentFilterStreamingDecorator<
  PersistentStatisticsMosaicFilter<TInputImage, TOutputImage, TInternalValueType> >
{
public:
  /** Standard Self typedef */
  typedef StreamingStatisticsMosaicFilter Self;
  typedef PersistentFilterStreamingDecorator
  <PersistentStatisticsMosaicFilter<TInputImage, TOutputImage, TInternalValueType> > Superclass;
  typedef itk::SmartPointer<Self>       Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Type macro */
  itkNewMacro(Self);

  /** Creation through object factory macro */
  itkTypeMacro(StreamingStatisticsMosaicFilter, PersistentFilterStreamingDecorator);

  /** Typedefs for statistics */
  typedef typename Superclass::FilterType::RealMatrixListType          RealMatrixListType;
  typedef typename Superclass::FilterType::LongMatrixType              LongMatrixType;
  typedef typename Superclass::FilterType::RealMatrixListObjectType    RealMatrixListObjectType;
  typedef typename Superclass::FilterType::LongMatrixObjectType        LongMatrixObjectType;

  /** Input image typedefs.  */
  typedef typename Superclass::FilterType::InputImageType              InputImageType;
  typedef typename Superclass::FilterType::InputImagePointer           InputImagePointer;
  typedef typename Superclass::FilterType::InputImagePointType         InputImagePointType;
  typedef typename Superclass::FilterType::InputImagePixelType         InputImagePixelType;
  typedef typename Superclass::FilterType::InputImageIndexType         InputImageIndexType;
  typedef typename Superclass::FilterType::InputImageSizeType          InputImageSizeType;
  typedef typename Superclass::FilterType::InputImageSpacingType       InputImageSpacingType;
  typedef typename Superclass::FilterType::InputImageInternalPixelType InputImageInternalPixelType;
  typedef typename Superclass::FilterType::InputImageRegionType        InputImageRegionType;

  /** Output image typedefs.  */
  typedef typename Superclass::FilterType::OutputImageType              OutputImageType;
  typedef typename Superclass::FilterType::OutputImagePointer           OutputImagePointer;
  typedef typename Superclass::FilterType::OutputImagePointType         OutputImagePointType;
  typedef typename Superclass::FilterType::OutputImagePixelType         OutputImagePixelType;
  typedef typename Superclass::FilterType::OutputImageIndexType         OutputImageIndexType;
  typedef typename Superclass::FilterType::OutputImageSizeType          OutputImageSizeType;
  typedef typename Superclass::FilterType::OutputImageSpacingType       OutputImageSpacingType;
  typedef typename Superclass::FilterType::OutputImageInternalPixelType OutputImageInternalPixelType;
  typedef typename Superclass::FilterType::OutputImageRegionType        OutputImageRegionType;

  /** Internal computing typedef support. */
  typedef typename Superclass::FilterType::InternalValueType       InternalValueType;
  typedef typename Superclass::FilterType::ContinuousIndexType     ContinuousIndexType;
  typedef typename Superclass::FilterType::InterpolatorType        InterpolatorType;
  typedef typename Superclass::FilterType::InterpolatorPointerType InterpolatorPointerType;
  typedef typename Superclass::FilterType::DefaultInterpolatorType DefaultInterpolatorType;
  typedef typename Superclass::FilterType::InternalImageType       InternalImageType;
  typedef typename Superclass::FilterType::InternalPixelType       InternalPixelType;
  typedef typename Superclass::FilterType::IteratorType            IteratorType;
  typedef typename Superclass::FilterType::ConstIteratorType       ConstIteratorType;
  typedef typename Superclass::FilterType::StreamingTraitsType     StreamingTraitsType;

  using Superclass::PushBackInput;
  void PushBackInput(InputImageType * input)
  {
    this->GetFilter()->PushBackInput(input);
  }

//  otbSetObjectMemberMacro(Filter, NoDataValue, RealType);
//  otbGetObjectMemberMacro(Filter, NoDataValue, RealType);

  /** Return the computed Means. */
  RealMatrixListType GetMeans() const
  {
    return this->GetFilter()->GetMeansOutput()->Get();
  }
  RealMatrixListObjectType* GetMeansOutput()
  {
    return this->GetFilter()->GetMeansOutput();
  }
  const RealMatrixListObjectType* GetMeansOutput() const
  {
    return this->GetFilter()->GetMeansOutput();
  }

  /** Return the computed Stds. */
  RealMatrixListType GetStds() const
  {
    return this->GetFilter()->GetStdsOutput()->Get();
  }
  RealMatrixListObjectType* GetStdsOutput()
  {
    return this->GetFilter()->GetStdsOutput();
  }
  const RealMatrixListObjectType* GetStdsOutput() const
  {
    return this->GetFilter()->GetStdsOutput();
  }

  /** Return the computed MeansOfProducts. */
  RealMatrixListType GetMeansOfProducts() const
  {
    return this->GetFilter()->GetMeansOfProductsOutput()->Get();
  }
  RealMatrixListObjectType* GetMeansOfProductsOutput()
  {
    return this->GetFilter()->GetMeansOfProductsOutput();
  }
  const RealMatrixListObjectType* GetMeansOfProductsOutput() const
  {
    return this->GetFilter()->GetMeansOfProductsOutput();
  }

  /** Return the computed Mins. */
  RealMatrixListType GetMins() const
  {
    return this->GetFilter()->GetMinsOutput()->Get();
  }
  RealMatrixListObjectType* GetMinsOutput()
  {
    return this->GetFilter()->GetMinsOutput();
  }
  const RealMatrixListObjectType* GetMinsOutput() const
  {
    return this->GetFilter()->GetMinsOutput();
  }

  /** Return the computed Maxs. */
  RealMatrixListType GetMaxs() const
  {
    return this->GetFilter()->GetMaxsOutput()->Get();
  }
  RealMatrixListObjectType* GetMaxsOutput()
  {
    return this->GetFilter()->GetMaxsOutput();
  }
  const RealMatrixListObjectType* GetMaxsOutput() const
  {
    return this->GetFilter()->GetMaxsOutput();
  }

  /** Return the computed Areas. */
  LongMatrixType GetAreas() const
  {
    return this->GetFilter()->GetAreasOutput()->Get();
  }
  LongMatrixObjectType* GetAreasOutput()
  {
    return this->GetFilter()->GetAreasOutput();
  }
  const LongMatrixObjectType* GetAreasOutput() const
  {
    return this->GetFilter()->GetAreasOutput();
  }

protected:
  /** Constructor */
  StreamingStatisticsMosaicFilter() {};
  /** Destructor */
  ~StreamingStatisticsMosaicFilter() ITK_OVERRIDE {}

private:
  StreamingStatisticsMosaicFilter(const Self &); //purposely not implemented
  void operator =(const Self&); //purposely not implemented
};

} // end namespace otb

#ifndef OTB_MANUAL_INSTANTIATION
#include "otbStreamingStatisticsMosaicFilter.txx"
#endif

#endif
