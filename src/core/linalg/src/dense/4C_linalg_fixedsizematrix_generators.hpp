/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of generators of different special fixed size matrices

\level 0
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINALG_FIXEDSIZEMATRIX_GENERATORS_HPP
#define FOUR_C_LINALG_FIXEDSIZEMATRIX_GENERATORS_HPP


#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  template <unsigned int size, typename value_type = double>
  LinAlg::Matrix<size, size, value_type> DiagonalMatrix(const value_type diagonal_value)
  {
    LinAlg::Matrix<size, size, value_type> diag_matrix(true);
    for (unsigned int i = 0; i < size; ++i) diag_matrix(i, i) = diagonal_value;
    return diag_matrix;
  }

  template <unsigned int size, typename value_type = double>
  LinAlg::Matrix<size, size, value_type> IdentityMatrix()
  {
    return DiagonalMatrix<size, value_type>(1);
  }
}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
