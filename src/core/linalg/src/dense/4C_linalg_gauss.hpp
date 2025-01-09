// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

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
    \tparam valtype  (in)    : type of values in the matrix - standard floating point or long
    precision
    \return determinant of system matrix
  */
  template <bool do_piv, unsigned dim, typename Valtype>
  Valtype gauss_elimination(
      Core::LinAlg::Matrix<dim, dim, Valtype>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<dim, 1, Valtype>& b,    ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<dim, 1, Valtype>& x     ///< (out)   : solution vector
  );


  /*!
    \brief computes a Gaussian elimination for a linear system of equations after infnorm scaling

    \tparam dim      (in)    : dimension of the matrix
    \return determinant of system matrix
  */
  template <unsigned dim>
  double scaled_gauss_elimination(Core::LinAlg::Matrix<dim, dim>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<dim, 1>& b,                                ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<dim, 1>& x                                 ///< (out)   : solution vector
  );


}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
