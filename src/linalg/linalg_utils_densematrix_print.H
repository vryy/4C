/*----------------------------------------------------------------------*/
/*! \file
* \file linalg_utils_densematrix_print.H

\brief A collection of dense matrix printing methods for namespace LINALG

\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*-----------------------------------------------------------------------*/

#ifndef SRC_LINALG_LINALG_UTILS_DENSEMATRIX_PRINT_H_
#define SRC_LINALG_LINALG_UTILS_DENSEMATRIX_PRINT_H_

#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "../drt_lib/drt_exporter.H"
#include "linalg_blocksparsematrix.H"

namespace LINALG
{
  // forward declaration
  // class Solver;

  /*!
            \brief Print sparsity pattern of a matrix to postscript file

            creates a file Epetra::Matrix.ps in current directory where the exact
            name of file depends on the exact type of class.

            \note works in parallel and serial!
            */
  void PrintSparsityToPostscript(const Epetra_RowMatrix& A);

  /*!
   \brief Print content of a sparse matrix in Matlab format

   The sparse matrix is printed in Matlab sparse format to
   a file with specified file name.

   \note - Works in parallel and serial!
   \note - When reading this in Matlab, entries are located according to
   their global row and column ID, even if corresponding DOFs are condensed
   from global system.
   */
  void PrintMatrixInMatlabFormat(
      std::string fname, const Epetra_CrsMatrix& A, const bool newfile = true);

  /*!
   \brief Print content of SerialDenseMatrix in Matlab format

   */
  void PrintSerialDenseMatrixInMatlabFormat(
      std::string fname, const Epetra_SerialDenseMatrix& A, const bool newfile = true);

  /*!
   \brief Print content of a block sparse matrix in Matlab format

   The block sparse matrix is printed in Matlab sparse format to
   a file with specified file name.

   \note - Works in parallel and serial!
   \note - When reading this in Matlab, entries are located according to
   their global row and column ID, even if corresponding DOFs are condensed
   from global system.
   */
  void PrintBlockMatrixInMatlabFormat(std::string fname, const BlockSparseMatrixBase& A);

  /*!
   \brief Print content of a vector in Matlab format

   The Vector V is printed in Matlab sparse format to
   a file with specified file name.

   */
  void PrintVectorInMatlabFormat(
      std::string fname, const Epetra_Vector& V, const bool newfile = true);

  /*!
   \brief Print content of a map in Matlab format

   The map is printed in Matlab numbering (+1) to
   a file with specified file name.

   */
  void PrintMapInMatlabFormat(std::string fname, const Epetra_Map& map, const bool newfile = true);

}  // namespace LINALG

#endif /* SRC_LINALG_LINALG_UTILS_DENSEMATRIX_PRINT_H_ */
