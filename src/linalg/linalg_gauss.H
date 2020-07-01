/*----------------------------------------------------------------------*/
/*! \file

\brief Gauss elimination for small nxn systems

\level 1

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
    \tparam valtype  (in)    : type of values in the matrix - standart floating point or long
    precision
    \return determinant of system matrix
  */
  template <bool do_piv, unsigned dim, typename valtype>
  valtype gaussElimination(LINALG::Matrix<dim, dim, valtype>& A,  ///< (in)    : system matrix
      LINALG::Matrix<dim, 1, valtype>& b,                         ///< (in)    : right-hand-side
      LINALG::Matrix<dim, 1, valtype>& x                          ///< (out)   : solution vector
  );


  /*!
    \brief computes a Gaussian elimination for a linear system of equations after infnorm scaling

    \tparam dim      (in)    : dimension of the matrix
    \return determinant of system matrix
  */
  template <unsigned dim>
  double scaledGaussElimination(LINALG::Matrix<dim, dim>& A,  ///< (in)    : system matrix
      LINALG::Matrix<dim, 1>& b,                              ///< (in)    : right-hand-side
      LINALG::Matrix<dim, 1>& x                               ///< (out)   : solution vector
  );


}  // namespace LINALG

#endif  // #ifndef LINALG_GAUSS_H
