#ifndef QuadraticallyConstrainedSimpleSolver_txx_
#define QuadraticallyConstrainedSimpleSolver_txx_

#include "otbQuadraticallyConstrainedSimpleSolver.h"

namespace otb {

template<class ValueType> QuadraticallyConstrainedSimpleSolver<ValueType>
::QuadraticallyConstrainedSimpleSolver() {
  m_WeightOfStandardDeviationTerm = itk::NumericTraits<ValueType>::One;
  oft = Cost_Function_rmse;
}

template<class ValueType> QuadraticallyConstrainedSimpleSolver<ValueType>
::~QuadraticallyConstrainedSimpleSolver() {
}

/*
 * Used to check layout topology consistency (Deep First Search)
 *
 * "m_AreaInOverlaps[i][j]>0" is equivalent to "images i and j are
 * overlapping with a non empty intersection (i.e. non null data)"
 */
template<class ValueType>
void QuadraticallyConstrainedSimpleSolver<ValueType>
::DFS(std::vector<bool> & marked, unsigned int s)
{

  // mark the s vertex
  marked[s] = true;

  // Get neighborhood
  for (unsigned int i = 0 ; i < m_AreaInOverlaps.rows() ; i++)
    {
    if (s!=i && m_AreaInOverlaps[s][i]>0 && !marked[i])
      {
      DFS(marked, i);
      }
    }

}

/*
 * Check input matrices dimensions, and layout consistency
 */
template<class ValueType>
void QuadraticallyConstrainedSimpleSolver<ValueType>
::CheckInputs() {

  // Check area matrix is non empty
  const unsigned int n = m_AreaInOverlaps.cols();
  if (n == 0)
    {
    itkExceptionMacro( << "Input area matrix has 0 elements");
    }

  bool inputMatricesAreValid = true;

  // Check "areas" and "means" matrices size
  if ((n != m_AreaInOverlaps.rows()) || (n != m_MeanInOverlaps.cols()) || (n != m_MeanInOverlaps.rows()))
    {
    inputMatricesAreValid = false;
    }

  // Check "std" matrix size
  if ((oft == Cost_Function_musig) || (oft==Cost_Function_weighted_musig) || (oft==Cost_Function_rmse))
    {
    if ((n != m_StandardDeviationInOverlaps.cols()) || (n != m_StandardDeviationInOverlaps.rows()))
      {
      inputMatricesAreValid = false;
      }
    }

  // Check "means of products" matrix size
  if (oft == Cost_Function_musig)
    {
    if ((n != m_MeanOfProductsInOverlaps.cols()) || (n != m_MeanOfProductsInOverlaps.rows()))
      {
      inputMatricesAreValid = false;
      }
    }

  if (!inputMatricesAreValid)
    {
    itkExceptionMacro( << "Input matrices must be square and have the same number of elements.");
    }

  // Images layout topology Check (using Depth First Search)
  std::vector<bool> marked;
  for (unsigned int i = 0 ; i < m_AreaInOverlaps.rows() ; i++)
    marked.push_back(false);
  bool valid = true;
  DFS(marked, 0);
  for (unsigned int i = 0 ; i < m_AreaInOverlaps.rows() ; i++)
    valid&=marked[i];

  if (!valid)
    {
    itkExceptionMacro( << "Inconsistent images layout: All images must be connected, at least with one overlap.")
    }
}

/*
 * Compute the objective function
 *
 * VNL is not sufficient: it has weak solving routines, and can not deal with QP subject to
 * a linear equality constraint plus lower limits (that is < and = linear constraints)
 * With vnl, we keep it simple and solve only the zero-y intercept linear case
 * but sometimes it fails because numerical instabilities. (vnl quadratic routines are not very reliable)
 *
 * With a good quadratic solver, we could e.g. use a general linear model (Xout = Xin*a+b)
 * Unfortunately it can be done with VNL so far. Tested & implemented successfully with
 * OOQP (Fastest, o(n)) and Quadprog++ (Fast, o(n)), and CGAL exact type solver (very slow, o(n^a) with a>1)
 * but has to rely on external dependencies...
 *
 * (see Cresson & St-Geours in IEEE JSTARS, "Natural color satellite image
 *  mosaicking using quadratic programming in decorrelated color space")
 *
 *
 */
template<class ValueType>
typename QuadraticallyConstrainedSimpleSolver<ValueType>::RealMatrixType
QuadraticallyConstrainedSimpleSolver<ValueType>
::GetQuadraticObjectiveMatrix()
{
  // Set STD matrix weight
  ValueType w;

  if (oft == Cost_Function_mu)
    {
    w = 0.0;
    }
  if (oft == Cost_Function_musig)
    {
    w = 1.0;
    }
  if (oft ==  Cost_Function_weighted_musig)
    {
    w = (ValueType) m_WeightOfStandardDeviationTerm;
    }

  const unsigned int n = m_MeanInOverlaps.cols();

  // Temporary matrices H, K, L
  RealMatrixType H(n,n,0), K(n,n,0), L(n,n,0);
  RealMatrixType H_RMSE(n,n,0);
  for (unsigned int i = 0 ; i < n ; i++)
    {
    for (unsigned int j = 0 ; j < n ; j++)
      {
      if (i==j)
        {
        // Diag (i=j)
        for (unsigned int k = 0 ; k < n ; k++)
          {
          if (i!=k)
            {
            H[i][j] += m_AreaInOverlaps[i][k] *
              (m_MeanInOverlaps[i][k]*m_MeanInOverlaps[i][k] + w*m_StandardDeviationInOverlaps[i][k]*
               m_StandardDeviationInOverlaps[i][k]);
            K[i][j] += m_AreaInOverlaps[i][k] * m_MeanInOverlaps[i][k];
            L[i][j] += m_AreaInOverlaps[i][k];
            H_RMSE[i][j] += m_AreaInOverlaps[i][k]*
              (m_MeanInOverlaps[i][k]*m_MeanInOverlaps[i][k]+m_StandardDeviationInOverlaps[i][k]*
               m_StandardDeviationInOverlaps[i][k]);
            }
          }
        }
      else
        {
        // Other than diag (i!=j)
        H[i][j] = -m_AreaInOverlaps[i][j] *
          (m_MeanInOverlaps[i][j]*m_MeanInOverlaps[j][i] + w*m_StandardDeviationInOverlaps[i][j]*
           m_StandardDeviationInOverlaps[j][i]);
        K[i][j] = -m_AreaInOverlaps[i][j] * m_MeanInOverlaps[i][j];
        L[i][j] = -m_AreaInOverlaps[i][j];
        H_RMSE[i][j] = -m_AreaInOverlaps[i][j] * m_MeanOfProductsInOverlaps[i][j];
        }
      }
    }

  if (oft == Cost_Function_rmse)
    {
    H = H_RMSE;
    }

  return H;

}

/*
 * QP Solving using vnl
 */
template<class ValueType>
void
QuadraticallyConstrainedSimpleSolver<ValueType>
::Solve()
{
  CheckInputs();

  // Number of images
  const unsigned int n = m_MeanInOverlaps.cols();

  // Objective function
  RealMatrixType Q = GetQuadraticObjectiveMatrix();
  RealVectorType g(n,0);

  // Constraint (Energy conservation)
  RealMatrixType A(1,n);
  RealVectorType b(1,0);
  for (unsigned int i = 0 ; i < n ; i++)
    {
    b[0] += m_AreaInOverlaps[i][i]*m_MeanInOverlaps[i][i];
    A[0][i] = m_AreaInOverlaps[i][i]*m_MeanInOverlaps[i][i];
    }

  RealVectorType x(n,1);
  // Change tol. to 0.01 is a quick hack to avoid numerical instability...
  bool solv = vnl_solve_qp_with_non_neg_constraints(Q,g,A,b,x,0.01);
  if (solv)
    {
    m_OutputCorrectionModel = RealVectorType(x);
    }
  else
    {
    itkWarningMacro( "vnl_solve_qp_with_non_neg_constraints failed." );
    m_OutputCorrectionModel.set_size(n);
    m_OutputCorrectionModel.fill(itk::NumericTraits<ValueType>::One);
    }
}

}

#endif
