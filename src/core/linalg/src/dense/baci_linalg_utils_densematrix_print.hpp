/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of dense matrix printing methods for namespace CORE::LINALG

\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_LINALG_UTILS_DENSEMATRIX_PRINT_HPP
#define FOUR_C_LINALG_UTILS_DENSEMATRIX_PRINT_HPP

#include "baci_config.hpp"

#include "baci_linalg_serialdensematrix.hpp"

#include <string>

BACI_NAMESPACE_OPEN

namespace CORE::LINALG
{
  /*!
  \brief Print content of SerialDenseMatrix in Matlab format

  @param[in] filename Name of output file
  @param[in] A Matrix to be printed
  @param[in] newfile Flage to force printing to a new file (instead of appeding to an existing file)
  */
  void PrintSerialDenseMatrixInMatlabFormat(
      std::string filename, const CORE::LINALG::SerialDenseMatrix& A, const bool newfile = true);

}  // namespace CORE::LINALG

BACI_NAMESPACE_CLOSE

#endif
