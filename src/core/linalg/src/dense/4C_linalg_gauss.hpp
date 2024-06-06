/*----------------------------------------------------------------------*/
/*! \file

\brief Gauss elimination for small nxn systems

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_LINALG_GAUSS_HPP
#define FOUR_C_LINALG_GAUSS_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
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
  valtype gaussElimination(Core::LinAlg::Matrix<dim, dim, valtype>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<dim, 1, valtype>& b,  ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<dim, 1, valtype>& x   ///< (out)   : solution vector
  );


  /*!
    \brief computes a Gaussian elimination for a linear system of equations after infnorm scaling

    \tparam dim      (in)    : dimension of the matrix
    \return determinant of system matrix
  */
  template <unsigned dim>
  double scaledGaussElimination(Core::LinAlg::Matrix<dim, dim>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<dim, 1>& b,                              ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<dim, 1>& x                               ///< (out)   : solution vector
  );


}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
