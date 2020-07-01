/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of algebraic creation methods for namespace LINALG

\level 0
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

  /*!
      \brief Creates MultiMapExtractor to split dofs at certain position

      We assume that each node possesses ndim dofs of one field and (optionally) remaining dofs of
      another field. The dof row map is thus split into two.

      The ndim dofs are assigned to map 0 (the other map) and the remaining dofs to map 1 (the
      condition map).

      \param dis : (in) Discretization
      \param ndim : (in) dimensions of map 0
      \param extractor : (out) ready made map splitter

      \author u.kue
      \date 02/08
     */
  void CreateMapExtractorFromDiscretization(
      const DRT::DiscretizationInterface& dis, int ndim, LINALG::MultiMapExtractor& extractor);

  /*!
    \brief Creates MapExtractor to split dofs at certain position

    We assume that each node possesses ndim dofs of one field and (optionally) remaining dofs of
    another field. The dof row map is thus split into two.

    The ndim dofs are assigned to map 0 (the other map) and the remaining dofs to map 1 (the
    condition map).

    Nodal dofs and total map come from dofset.

    \param dis : (in) Discretization
    \param dofset : (in) Degree of freedom set
    \param ndim : (in) dimensions of map 0
    \param extractor : (out) ready made map splitter

    \author u.kue
    \date 02/08
   */
  void CreateMapExtractorFromDiscretization(const DRT::Discretization& dis,
      const DRT::DofSetInterface& dofset, int ndim, LINALG::MapExtractor& extractor);


  /*!
    \brief Creates MultiMapExtractor to split dofs at certain position

    We assume that each node possesses ndim_field1 + ndim_field2 dofs. The dof row map is thus
    split into two.

    The ndim_field1 dofs are assigned to map 0 (the other map) and the ndim_field2 to map 1 (the
    condition map).

    \param dis : (in) Discretization
    \param ndim_field1 : (in) dimensions of map 0
    \param ndim_field2 : (in) dimensions of map 1
    \param extractor : (out) ready made map splitter

    \author schott
    \date 12/11
   */
  void CreateMapExtractorFromDiscretization(const DRT::Discretization& dis, int ndim_field1,
      int ndim_field2, LINALG::MultiMapExtractor& extractor);

}  // namespace LINALG

#endif /* SRC_LINALG_LINALG_UTILS_SPARSE_ALGEBRA_CREATE_H_ */
