/*----------------------------------------------------------------------*/
/*!
\file linalg_gauss.H

\brief Gauss elimination for small nxn systems

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#ifndef LINALG_GAUSS_H
#define LINALG_GAUSS_H

#include "../linalg/linalg_fixedsizematrix.H"

namespace LINALG
{

/*!
  \brief computes a Gaussian elimination for a linear system of equations

  \tparam do_piv   (in)    : do_piv = true does pivoting, do_piv = false does not do pivoting
  \tparam dim      (in)    : dimension of the matrix
  \return determinant of system matrix
*/
template<bool do_piv, unsigned dim>
double gaussElimination(
  typename LINALG::Matrix<dim, dim> & A,  ///< (in)    : system matrix
  typename LINALG::Matrix<dim, 1>   & b,  ///< (in)    : right-hand-side
  typename LINALG::Matrix<dim, 1>   & x   ///< (out)   : solution vector
  );


/*!
  \brief computes a Gaussian elimination for a linear system of equations after infnorm scaling

  \tparam dim      (in)    : dimension of the matrix
  \return determinant of system matrix
*/
template<unsigned dim>
double scaledGaussElimination(
  typename LINALG::Matrix<dim, dim> & A,  ///< (in)    : system matrix
  typename LINALG::Matrix<dim, 1>   & b,  ///< (in)    : right-hand-side
  typename LINALG::Matrix<dim, 1>   & x   ///< (out)   : solution vector
  );


}

#endif // #ifndef LINALG_GAUSS_H
