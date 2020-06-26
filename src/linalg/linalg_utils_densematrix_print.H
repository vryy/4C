/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of dense matrix printing methods for namespace LINALG

\level 0
\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/

#ifndef SRC_LINALG_LINALG_UTILS_DENSEMATRIX_PRINT_H_
#define SRC_LINALG_LINALG_UTILS_DENSEMATRIX_PRINT_H_

#include <string>

#include "Epetra_SerialDenseMatrix.h"

namespace LINALG
{
  /*!
  \brief Print content of SerialDenseMatrix in Matlab format

  @param[in] filename Name of output file
  @param[in] A Matrix to be printed
  @param[in] newfile Flage to force printing to a new file (instead of appeding to an existing file)
  */
  void PrintSerialDenseMatrixInMatlabFormat(
      std::string filename, const Epetra_SerialDenseMatrix& A, const bool newfile = true);

}  // namespace LINALG

#endif /* SRC_LINALG_LINALG_UTILS_DENSEMATRIX_PRINT_H_ */
