/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of dense matrix printing methods for namespace CORE::LINALG

\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_LINALG_UTILS_DENSEMATRIX_PRINT_H
#define BACI_LINALG_UTILS_DENSEMATRIX_PRINT_H

#include <string>

#include <Epetra_SerialDenseMatrix.h>

namespace CORE::LINALG
{
  /*!
  \brief Print content of SerialDenseMatrix in Matlab format

  @param[in] filename Name of output file
  @param[in] A Matrix to be printed
  @param[in] newfile Flage to force printing to a new file (instead of appeding to an existing file)
  */
  void PrintSerialDenseMatrixInMatlabFormat(
      std::string filename, const Epetra_SerialDenseMatrix& A, const bool newfile = true);

}  // namespace CORE::LINALG

#endif  // LINALG_UTILS_DENSEMATRIX_PRINT_H
