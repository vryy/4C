/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of sparse matrix printing methods for namespace Core::LinAlg

\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_LINALG_UTILS_SPARSE_ALGEBRA_PRINT_HPP
#define FOUR_C_LINALG_UTILS_SPARSE_ALGEBRA_PRINT_HPP

#include "4C_config.hpp"

#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCPDecl.hpp>

FOUR_C_NAMESPACE_OPEN

// Forward declarations
namespace Core::LinAlg
{
  class BlockSparseMatrixBase;
}  // namespace Core::LinAlg

namespace Core::LinAlg
{
  //! Print content of @p sparsematrix in Matlab format to file @p filename. Create new file or
  //! overwrite exisiting one if @p newfile is true
  void PrintMatrixInMatlabFormat(
      const std::string& filename, const Epetra_CrsMatrix& sparsematrix, const bool newfile = true);

  //! Print content of @p blockmatrix in Matlab format to file @p filename
  void PrintBlockMatrixInMatlabFormat(
      const std::string& filename, const BlockSparseMatrixBase& blockmatrix);

  //! Print content of @p vector in Matlab format to file @p filename. Create new file or overwrite
  //! exisiting one if @p newfile is true
  void PrintVectorInMatlabFormat(
      const std::string& filename, const Epetra_Vector& vector, const bool newfile = true);

  //! Print content of @p map in Matlab format to file @p filename. Create new file or overwrite
  //! exisiting one if @p newfile is true
  void PrintMapInMatlabFormat(
      const std::string& filename, const Epetra_Map& map, const bool newfile = true);

}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
