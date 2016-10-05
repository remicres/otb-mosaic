#ifndef QuadraticallyConstrainedSimpleSolver_H_
#define QuadraticallyConstrainedSimpleSolver_H_

#include "itkObjectFactory.h"
#include "itkLightObject.h"
#include "itkNumericTraits.h"
#include <vnl/vnl_matrix.h>
#include "vnl/algo/vnl_solve_qp.h"

namespace otb {

/**
 * \class QuadraticallyConstrainedSimpleSolver
 * \brief Solves the optimisation problem for radiometric harmonisation of multiple
 * overlapping images.
 *
 * This solver inputs statistics of the overlapping images, and produces a
 * zero-y intercept linear correction model. Various cost functions can be
 * employed: RMSE based, Mean based, Mean+Standard deviation based, and Mean
 * + weighted Standard deviation bases.
 *
 * Inputs:
 * -N x N Matrix of mean of overlaps ij
 * -N x N Matrix of standard deviation of overlaps ij
 * -N x N Matrix of area of overlaps ij
 * -N x N Matrix of mean of pixels products of overlaps ij
 *
 * For all i and j, m_{ij} = stats of image i in overlap ij
 *
 * Output:
 * N x 1 Vector of scales to apply to images
 *
 *  For more details, see Cresson & St-Geours "Natural color satellite image
 *  mosaicking using quadratic programming in decorrelated color space"
 *  in IEEE JSTARS Issue 99 Volume PP, July 2015
 *
 *  http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=7154397&filter%3DAND%28p_IS_Number%3A4609444%29
 *
 *
 * \ingroup OTBMosaic
 */

template<class ValueType>
class ITK_EXPORT QuadraticallyConstrainedSimpleSolver : public itk::LightObject
{
public:

  /** Standard class typedef */
  typedef QuadraticallyConstrainedSimpleSolver Self;
  typedef itk::LightObject                     Superclass;
  typedef itk::SmartPointer<Self>              Pointer;
  typedef itk::SmartPointer<const Self>        ConstPointer;

  /** Runtime information support. */
  itkTypeMacro(QuadraticallyConstrainedSimpleSolver, LightObject);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Typedefs */
  typedef vnl_matrix<ValueType> RealMatrixType;
  typedef vnl_matrix<long>      LongMatrixType;
  typedef vnl_vector<ValueType> RealVectorType;

  /** Enum for objective function type */
  enum ObjectiveFunctionType
    {
    Cost_Function_rmse,          // Root mean square error based
    Cost_Function_musig,         // Mean and standard deviation based
    Cost_Function_mu,            // Mean based
    Cost_Function_weighted_musig // Mean and weighted standard deviation based
    };

  /** mean in overlaps matrix getters/setters */
  void SetMeanInOverlaps(RealMatrixType matrix) {
    m_MeanInOverlaps = RealMatrixType(matrix);
  }

  RealMatrixType GetMeanInOverlaps() {
    return m_MeanInOverlaps;
  }

  /** standard deviation in overlaps matrix getters/setters */
  void SetStandardDeviationInOverlaps(RealMatrixType matrix) {
    m_StandardDeviationInOverlaps = RealMatrixType(matrix);
  }

  RealMatrixType GetStandardDeviationInOverlaps() {
    return m_StandardDeviationInOverlaps;
  }

  /** area in overlaps matrix getters/setters */
  void SetAreaInOverlaps(RealMatrixType matrix) {
    m_AreaInOverlaps = RealMatrixType(matrix);
  }

  void SetAreaInOverlaps(LongMatrixType matrix) {
    const unsigned int n = matrix.cols();

    m_AreaInOverlaps = RealMatrixType (n,n,0);
    for (unsigned int i = 0 ; i < n ; i++)
      {
      for (unsigned int j = 0 ; j < n ; j++)
        m_AreaInOverlaps[i][j] = 1.0 * matrix[i][j];
      }
  }

  RealMatrixType GetAreaInOverlaps() {
    return m_AreaInOverlaps;
  }

  /** mean of pixels products in overlaps matrix getters/setters */
  void SetMeanOfProductsInOverlaps(RealMatrixType matrix) {
    m_MeanOfProductsInOverlaps = RealMatrixType(matrix);
  }

  RealMatrixType GetMeanOfProductsInOverlaps() {
    return m_MeanOfProductsInOverlaps;
  }

  /** get the output correction model */
  RealVectorType GetOutputCorrectionModel() {
    return m_OutputCorrectionModel;
  }

  /**
   * Set the STD weight in harmonization
   * if value is near 0, importance is accorded to MEAN
   * if value is 1, importance is the same than MEAN
   * if value is higher than 1, importance is accorder to STD
   */
  void SetWeightOfStandardDeviationTerm(ValueType weight) {
    m_WeightOfStandardDeviationTerm = weight;
  }

  ValueType GetWeightOfStandardDeviationTerm() {
    return m_WeightOfStandardDeviationTerm;
  }

  /** Solving routine */
  void Solve();

  /** Set the cost function type */
  void SetMeanBased(){
    oft = Cost_Function_mu;
  }

  void SetMeanAndStandardDeviationBased(){
    oft = Cost_Function_musig;
  }

  void SetRMSEBased(){
    oft = Cost_Function_rmse;
  }

  void SetWeightedMeanAndStandardDeviationBased(){
    oft = Cost_Function_weighted_musig;
  }

protected:

  QuadraticallyConstrainedSimpleSolver();
  virtual ~QuadraticallyConstrainedSimpleSolver();

private:

  // Check inputs
  void CheckInputs(void);

  // Deep First Search
  void DFS(bool * marked, unsigned int s);

  // Compute the objective matrix
  vnl_matrix<ValueType> GetQuadraticObjectiveMatrix();

  // Input
  RealMatrixType m_MeanInOverlaps;
  RealMatrixType m_StandardDeviationInOverlaps;
  RealMatrixType m_AreaInOverlaps;
  RealMatrixType m_MeanOfProductsInOverlaps;

  // Params
  ValueType m_WeightOfStandardDeviationTerm; // could be manually tuned for
                                             // different results

  // Output correction models
  RealVectorType m_OutputCorrectionModel;

  // objective funciton type (enum)
  ObjectiveFunctionType oft;

};

} /* namespace otb */
#ifndef OTB_MANUAL_INSTANTIATION
#include "otbQuadraticallyConstrainedSimpleSolver.txx"
#endif
#endif /* QuadraticallyConstrainedSimpleSolver_H_ */
