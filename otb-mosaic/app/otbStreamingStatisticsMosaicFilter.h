#ifndef __StreamingStatisticsMosaicFilter_H
#define __StreamingStatisticsMosaicFilter_H

#include "otbStreamingMosaicFilterBase.h"

namespace otb
{
/** \class StreamingStatisticsMosaicFilter
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
 * \ingroup Multithreaded
 */
template <class TInputImage, class TOutputImage=TInputImage, class TInternalValueType=double>
class ITK_EXPORT StreamingStatisticsMosaicFilter : public otb::StreamingMosaicFilterBase<TInputImage, TOutputImage,
                                                                                         TInternalValueType>
{
public:

  /** Standard Self typedef */
  typedef StreamingStatisticsMosaicFilter                                               Self;
  typedef otb::StreamingMosaicFilterBase<TInputImage, TOutputImage, TInternalValueType> Superclass;
  typedef itk::SmartPointer<Self>                                                       Pointer;
  typedef itk::SmartPointer<const Self>                                                 ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(StreamingStatisticsMosaicFilter, StreamingMosaicFilterBase);

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
  typedef typename Superclass::IteratorType            IteratorType;
  typedef typename Superclass::ConstIteratorType       ConstIteratorType;
  typedef typename Superclass::StreamingTraitsType     StreamingTraitsType;

  /** Get Stats methods */
  vnl_matrix<InternalValueType> GetMean(unsigned int band);

  vnl_matrix<InternalValueType> GetProdMean(unsigned int band);

  vnl_matrix<InternalValueType> GetStDev(unsigned int band);

  vnl_vector<InputImageInternalPixelType> GetMin(unsigned int band);

  vnl_vector<InputImageInternalPixelType> GetMax(unsigned int band);

  vnl_matrix<long> GetAreaInPixels();

protected:
  StreamingStatisticsMosaicFilter();
  virtual ~StreamingStatisticsMosaicFilter() {
  }

  /** StreamingStatisticsMosaicFilter can be implemented as a multithreaded filter.
   * Therefore, this implementation provides a ThreadedGenerateData() routine
   * which is called for each processing thread. The output image data is
   * allocated automatically by the superclass prior to calling
   * ThreadedGenerateData().  ThreadedGenerateData can only write to the
   * portion of the output image specified by the parameter
   * "outputRegionForThread"
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData()  */

  /** Overrided methods */
  virtual void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, itk::ThreadIdType threadId );

  virtual void GenerateOutputInformation();

  virtual void AfterThreadedGenerateData();

  virtual void BeforeThreadedGenerateData();

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
      m_count = vnl_vector<long>(other.m_count);
      m_sum = vnl_matrix<InternalValueType>(other.m_sum);
      m_cosum = vnl_matrix<InternalValueType>(other.m_cosum);
      m_sqSum = vnl_matrix<InternalValueType>(other.m_sqSum);
      m_min = vnl_matrix<InternalValueType>(other.m_min);
      m_max = vnl_matrix<InternalValueType>(other.m_max);
    }

    /* Clear routine: Resize at the specified dimension and clear values */
    void Clear(unsigned int nbOfBands, unsigned int nbOfSamples)
    {
      InternalValueType zeroValue = itk::NumericTraits<InternalValueType>::Zero;
      InternalValueType supValue = itk::NumericTraits<InternalValueType>::max();
      InternalValueType infValue = itk::NumericTraits<InternalValueType>::NonpositiveMin();

      m_count = vnl_vector<long>(nbOfSamples,0);
      m_sum = vnl_matrix<InternalValueType>(nbOfBands,nbOfSamples,zeroValue);
      m_cosum = vnl_matrix<InternalValueType>(nbOfBands,nbOfSamples,zeroValue);
      m_sqSum = vnl_matrix<InternalValueType>(nbOfBands,nbOfSamples,zeroValue);
      m_min = vnl_matrix<InternalValueType>(nbOfBands,nbOfSamples,supValue);
      m_max = vnl_matrix<InternalValueType>(nbOfBands,nbOfSamples,infValue);
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

    vnl_matrix<InternalValueType> m_sum;
    vnl_matrix<InternalValueType> m_sqSum;
    vnl_matrix<InternalValueType> m_cosum;
    vnl_matrix<InternalValueType> m_min;
    vnl_matrix<InternalValueType> m_max;
    vnl_vector<long>              m_count;
  };

  // Results
  std::vector<ThreadResultsContainer> m_ThreadResults;
  ThreadResultsContainer              m_FinalResults;

private:
  StreamingStatisticsMosaicFilter(const Self&); //purposely not implemented
  void operator=(const Self&);                  //purposely not implemented

}; // end of class

} // end namespace otb

#ifndef OTB_MANUAL_INSTANTIATION
#include "otbStreamingStatisticsMosaicFilter.hxx"
#endif

#endif
