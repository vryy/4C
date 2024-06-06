/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of dense matrix printing methods for namespace Core::LinAlg

\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_LINALG_UTILS_DENSEMATRIX_PRINT_HPP
#define FOUR_C_LINALG_UTILS_DENSEMATRIX_PRINT_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensematrix.hpp"

#include <string>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  /*!
  \brief Print content of SerialDenseMatrix in Matlab format

  @param[in] filename Name of output file
  @param[in] A Matrix to be printed
  @param[in] newfile Flage to force printing to a new file (instead of appeding to an existing file)
  */
  void PrintSerialDenseMatrixInMatlabFormat(
      std::string filename, const Core::LinAlg::SerialDenseMatrix& A, const bool newfile = true);

}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
