/*---------------------------------------------------------------------*/
/*! \file

\brief Conversion of LINALG matrix with CLN type to double and vice versa

\level 3


*----------------------------------------------------------------------*/

#ifndef BACI_UTILS_CLN_MATRIX_CONVERSION_HPP
#define BACI_UTILS_CLN_MATRIX_CONVERSION_HPP

#include "baci_config.hpp"

#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_utils_exceptions.hpp"

#include <cstddef>
#include <sstream>
#include <unordered_map>
#include <vector>

BACI_NAMESPACE_OPEN

namespace CORE::CLN
{
  //! convert LINALG matrix with CLN values to LINALG matrix with double values
  template <unsigned int num_row, unsigned int num_col>
  void ConvClnDouble(const LINALG::Matrix<num_row, num_col, ClnWrapper>& in,
      LINALG::Matrix<num_row, num_col, double>& out)
  {
    for (unsigned int idx = 0; idx < in.numRows() * in.numCols(); ++idx)
    {
      out.A()[idx] = cln::double_approx(in.A()[idx].Value());
    }
  }

  //! convert LINALG matrix with double values to LINALG matrix with CLN values
  template <unsigned int num_row, unsigned int num_col>
  void ConvDoubleCLN(const LINALG::Matrix<num_row, num_col, double>& in,
      LINALG::Matrix<num_row, num_col, ClnWrapper>& out, const int precision = 20)
  {
    for (unsigned int idx = 0; idx < in.numRows() * in.numCols(); ++idx)
    {
      ClnWrapper clnnum;
      // zeros do not convert properly to CLN (loss of precision)
      if ((in.A()[idx] == 0.0))
      {
        // returning the cached value from the ClnWrapper cln table
        clnnum = 0.0;
      }
      else
        clnnum = cln::cl_float(in.A()[idx], cln::float_format(precision));
      out.A()[idx] = clnnum;
    }
  }

  //! convert LINALG matrix with CLN values to a new LINALG matrix with CLN values with different
  //! precision
  template <unsigned int num_row, unsigned int num_col>
  void UpdatePresicion(const LINALG::Matrix<num_row, num_col, ClnWrapper>& in,
      LINALG::Matrix<num_row, num_col, ClnWrapper>& out, const int precision = 20)
  {
    for (unsigned int idx = 0; idx < in.numRows() * in.numCols(); ++idx)
    {
      ClnWrapper clnnum;
      // zeros do not convert properly to CLN (loss of precision)
      if ((in.A()[idx] == 0.0))
      {
        // returning the cached value from the ClnWrapper cln table
        clnnum = 0.0;
      }
      else
        clnnum = cln::cl_float(in.A()[idx].Value(), cln::float_format(precision));
      out.A()[idx] = clnnum;
    }
  }

}  // namespace CORE::CLN

BACI_NAMESPACE_CLOSE

#endif
