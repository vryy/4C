/*----------------------------------------------------------------------*/
/*! \file

\brief Gauss elimination for small nxn systems

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_LINALG_GAUSS_HPP
#define BACI_LINALG_GAUSS_HPP

#include "baci_config.hpp"

#include "baci_linalg_fixedsizematrix.hpp"

BACI_NAMESPACE_OPEN

namespace CORE::LINALG
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
  valtype gaussElimination(CORE::LINALG::Matrix<dim, dim, valtype>& A,  ///< (in)    : system matrix
      CORE::LINALG::Matrix<dim, 1, valtype>& b,  ///< (in)    : right-hand-side
      CORE::LINALG::Matrix<dim, 1, valtype>& x   ///< (out)   : solution vector
  );


  /*!
    \brief computes a Gaussian elimination for a linear system of equations after infnorm scaling

    \tparam dim      (in)    : dimension of the matrix
    \return determinant of system matrix
  */
  template <unsigned dim>
  double scaledGaussElimination(CORE::LINALG::Matrix<dim, dim>& A,  ///< (in)    : system matrix
      CORE::LINALG::Matrix<dim, 1>& b,                              ///< (in)    : right-hand-side
      CORE::LINALG::Matrix<dim, 1>& x                               ///< (out)   : solution vector
  );


}  // namespace CORE::LINALG

BACI_NAMESPACE_CLOSE

#endif  // LINALG_GAUSS_H
