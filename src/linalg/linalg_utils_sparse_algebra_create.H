/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of algebraic creation methods for namespace LINALG

\level 0
\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/
#ifndef SRC_LINALG_LINALG_UTILS_SPARSE_ALGEBRA_CREATE_H_
#define SRC_LINALG_LINALG_UTILS_SPARSE_ALGEBRA_CREATE_H_

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_exporter.H"
#include "linalg_blocksparsematrix.H"

namespace LINALG
{
  // forward declaration
  // class Solver;

  /*!
   \brief Create a new Epetra_CrsMatrix and return RefcountPtr to it

   \param rowmap (in): row map of matrix
   \param npr (in): estimated number of entries per row.
   (need not be exact, better should be too big rather then too small)
   */
  Teuchos::RCP<Epetra_CrsMatrix> CreateMatrix(const Epetra_Map& rowmap, const int npr);

  /*!
   \brief Create a new Epetra_Vector and return RefcountPtr to it

   \param rowmap (in): row map of vector
   \param init (in): initializa vector to zero upon construction
   */
  Teuchos::RCP<Epetra_Vector> CreateVector(const Epetra_BlockMap& rowmap, const bool init = true);

  /*!
   \brief Create a new Epetra_MultiVector and return RefcountPtr to it

   \param rowmap (in): row map of vector
   \param rowmap (in): number of vectors
   \param init (in): initializa vector to zero upon construction
   */
  Teuchos::RCP<Epetra_MultiVector> CreateMultiVector(
      const Epetra_BlockMap& rowmap, const int numrows, const bool init = true);

  /*!
   \brief Create an Epetra_Map from a set of gids

   This is one of the basic operations that is needed every so often.

   \param gids The local gids of this map
   \param comm The map's communicator
   */
  Teuchos::RCP<Epetra_Map> CreateMap(const std::set<int>& gids, const Epetra_Comm& comm);

  /*!
   \brief Create an Epetra_Map from a vector of gids

   This is one of the basic operations that is needed every so often.

   \param gids The local gids of this map
   \param comm The map's communicator
   */
  Teuchos::RCP<Epetra_Map> CreateMap(const std::vector<int>& gids, const Epetra_Comm& comm);

}  // namespace LINALG

#endif /* SRC_LINALG_LINALG_UTILS_SPARSE_ALGEBRA_CREATE_H_ */
