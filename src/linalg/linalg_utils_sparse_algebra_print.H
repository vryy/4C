/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of sparse matrix printing methods for namespace LINALG

\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef SRC_LINALG_LINALG_UTILS_SPARSE_ALGEBRA_PRINT_H_
#define SRC_LINALG_LINALG_UTILS_SPARSE_ALGEBRA_PRINT_H_

#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "../drt_lib/drt_exporter.H"

// Forward declarations
namespace LINALG
{
  class BlockSparseMatrixBase;
}  // namespace LINALG

namespace LINALG
{
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

#ifndef TRILINOS_Q1_2015
  /*!
  \brief Write matrix to ASCII-file via Xpetra::IO

  Write using \c Xpetra::IO<>::Write() to be compatible with \c Xpetra::IO<>::Read().

  \note This is a utility function for debugging, e.g. to dump matrices and read them in some
  other code.

  @param[in] filename Name of output file
  @param[in] matrix Matrix to be written
  */
  void WriteEpetraCrsMatrixAsXpetra(
      const std::string& filename, Teuchos::RCP<Epetra_CrsMatrix> matrix);

  /*!
  \brief Write MultiVector to ASCII-file via Xpetra::IO

  Write using \c Xpetra::IO<>::Write() to be compatible with \c Xpetra::IO<>::Read().

  \note This is a utility function for debugging, e.g. to dump MultiVectors and read them in some
  other code.

  @param[in] filename Name of output file
  @param[in] vec MultiVector to be written
  */
  void WriteEpetraMultiVectorAsXpetra(
      const std::string& filename, Teuchos::RCP<Epetra_MultiVector> vec);
#endif

}  // namespace LINALG

#endif /* SRC_LINALG_LINALG_UTILS_SPARSE_ALGEBRA_PRINT_H_ */
