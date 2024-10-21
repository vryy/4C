// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_utils_densematrix_print.hpp"

#include <fstream>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::LinAlg::print_serial_dense_matrix_in_matlab_format(
    std::string filename, const Core::LinAlg::SerialDenseMatrix& A, const bool newfile)
{
  std::ofstream os;

  // open file for writing
  if (newfile)
    os.open(filename.c_str(), std::fstream::trunc);
  else
    os.open(filename.c_str(), std::fstream::ate | std::fstream::app);

  const int NumMyRows = A.numRows();
  const int NumMyColumns = A.numCols();

  for (int i = 0; i < NumMyRows; i++)
  {
    for (int j = 0; j < NumMyColumns; j++)
    {
      os << std::setw(10) << i + 1;  // increase index by one for matlab
      os << std::setw(10) << j + 1;  // increase index by one for matlab
      os << std::setw(30) << std::setprecision(16) << std::scientific << A(i, j);
      os << std::endl;
    }
  }

  os << std::flush;

  // close file
  os.close();

  // just to be sure
  if (os.is_open()) os.close();
}

FOUR_C_NAMESPACE_CLOSE
