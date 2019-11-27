/*----------------------------------------------------------------------*/
/*! \file
 * \file linalg_utils_densematrix_manipulation.H

\brief A collection of dense matrix manipulation methods for namespace LINALG

\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*-----------------------------------------------------------------------*/

#ifndef SRC_LINALG_LINALG_UTILS_DENSEMATRIX_MANIPULATION_H_
#define SRC_LINALG_LINALG_UTILS_DENSEMATRIX_MANIPULATION_H_

#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_exporter.H"
#include "linalg_blocksparsematrix.H"

// forward declarations
class Epetra_IntVector;

namespace IO
{
  class DiscretizationWriter;
}

namespace LINALG
{
  // forward declaration
  // class Solver;

  /*!
       \brief Export a vector to a different map

       Values of source are copied to target where maps don't have to match.
       Prerequisite: Either the map of source OR the map of target has to be unique
       (will be tested)
       \warning When source is overlapping (and therefore target is unique), values
       in the overlapping region are inserted into the target on a first come
       first serve basis, meaning they should be equal in the source to
       be deterministic
       \param source (in) : source vector values are taken from
       \param target (out): target vector values will be inserted in
       */
  void Export(const Epetra_MultiVector& source, Epetra_MultiVector& target);

  /*!
   \brief Export a vector to a different map

   Values of source are copied to target where maps don't have to match.
   Prerequisite: Either the map of source OR the map of target has to be unique
   (will be tested)
   \warning When source is overlapping (and therefore target is unique), values
   in the overlapping region are inserted into the target on a first come
   first serve basis, meaning they should be equal in the source to
   be deterministic
   \param source (in) : source vector values are taken from
   \param target (out): target vector values will be inserted in
   */
  void Export(const Epetra_IntVector& source, Epetra_IntVector& target);

  /*
   \brief Insert a list of elements in a given global row of the matrix with
   out-of-bound and error check

   \param A (in):           Considered matrix
   \param GlobalRow(in):    Row number (in global coordinates) to put elements.
   \param NumEntries(in):   Number of entries.
   \param Values(in):       Values to enter.
   \param Indices(in):      Global column indices corresponding to values.
   */
  void InsertGlobalValues(Teuchos::RCP<Epetra_CrsMatrix> mat, int GlobalRow, int NumEntries,
      double* Values, int* Indices);

  /*! \brief Extract a partial Epetra_Vector from a given source vector
   *         on each proc without communication
   *
   *  This methods uses a given partial map to create the partial target vector.
   *
   *  \param source     (in) : source vector ( read-only )
   *  \param target_map (in) : map of the new target vector ( read-only )
   *
   *  \return the extracted partial Epetra_Vector as RCP
   *
   *  \author hiermeier \date 03/17 */
  Teuchos::RCP<Epetra_Vector> ExtractMyVector(
      const Epetra_Vector& source, const Epetra_Map& target_map);

  /*! \brief Extract a partial Eptra_Vector from a given source vector
   *         on each proc without communication
   *
   *  \param source (in) : source vector ( read-only )
   *  \param target (out): this target vector is going to be filled
   *
   *  \author hiermeier \date 03/17 */
  void ExtractMyVector(const Epetra_Vector& source, Epetra_Vector& target);

  /*! \brief Extract a partial Epetra_Vector from a given source vector
   *         and scale the entries on each proc without communication
   *
   *  \param scalar_source (in) : scale the target entries by this factor
   *  \param source        (in) : source vector ( read-only )
   *  \param scalar_target (in) : scale the source entries by this factor
   *  \param target        (out): this target vector is going to be filled with
   *                              the combined scaled values of the source and
   *                              target vectors
   *
   *  \author hiermeier \date 03/17 */
  void ExtractMyVector(double scalar_source, const Epetra_Vector& source, double scalar_target,
      Epetra_Vector& target);

  /*!
   \brief Call FillComplete on a Epetra_CrsMatrix (for square matrices only!)
   */
  void Complete(Epetra_CrsMatrix& A);

  /*!
   \brief Call FillComplete on a Epetra_CrsMatrix (for square matrices only!)

   This is the Teuchos::RCP wrapper of the above method.
   */
  inline void Complete(Teuchos::RCP<Epetra_CrsMatrix> A)
  {
    LINALG::Complete(*A);
    return;
  }

  /*!
   \brief Call FillComplete on a Epetra_CrsMatrix (for rectangular and square matrices)
   */
  void Complete(Epetra_CrsMatrix& A, const Epetra_Map& domainmap, const Epetra_Map& rangemap);

  /*!
   \brief Call FillComplete on a Epetra_CrsMatrix (for rectangular and square matrices)

   This is the Teuchos::RCP wrapper of the above method.
   */
  inline void Complete(
      Teuchos::RCP<Epetra_CrsMatrix> A, const Epetra_Map& domainmap, const Epetra_Map& rangemap)
  {
    LINALG::Complete(*A, domainmap, rangemap);
    return;
  }

  /*!
   \brief split a matrix into a 2x2 block system where the rowmap of one of the blocks is given

   Splits a given matrix into a 2x2 block system where the rowmap of one of the blocks is given
   on input. Blocks A11 and A22 are assumed to be square.
   All values on entry have to be Teuchos::null except the given rowmap and matrix A.
   Note that either A11rowmap or A22rowmap or both have to be nonzero. In case
   both rowmaps are supplied they have to be an exact and nonoverlapping split of A->RowMap().
   Matrix blocks are FillComplete() on exit.

   \param A         : Matrix A on input
   \param A11rowmap : rowmap of A11 or null
   \param A22rowmap : rowmap of A22 or null
   \param A11       : on exit matrix block A11
   \param A12       : on exit matrix block A12
   \param A21       : on exit matrix block A21
   \param A22       : on exit matrix block A22
   */
  bool SplitMatrix2x2(Teuchos::RCP<Epetra_CrsMatrix> A, Teuchos::RCP<Epetra_Map>& A11rowmap,
      Teuchos::RCP<Epetra_Map>& A22rowmap, Teuchos::RCP<Epetra_CrsMatrix>& A11,
      Teuchos::RCP<Epetra_CrsMatrix>& A12, Teuchos::RCP<Epetra_CrsMatrix>& A21,
      Teuchos::RCP<Epetra_CrsMatrix>& A22);

  /*!
   \brief split a matrix into a 2x2 block system

   Splits a given matrix into a 2x2 block system. All values on entry have to be
   Teuchos::null except the given rowmap(s) / domainmap(s) and matrix A.
   Note that either A11rowmap or A22rowmap or both have to be nonzero!
   Note that either A11domainmap or A22domainmap or both have to be nonzero!
   In case both rowmaps / domainmaps are supplied they have to be an exact and
   nonoverlapping split of A->RowMap() / A->DomainMap().
   Matrix blocks are FillComplete() on exit.

   \param A            : Matrix A on input
   \param A11rowmap    : rowmap of A11 or null
   \param A22rowmap    : rowmap of A22 or null
   \param A11domainmap : domainmap of A11 or null
   \param A22domainmap : domainmap of A22 or null
   \param A11          : on exit matrix block A11
   \param A12          : on exit matrix block A12
   \param A21          : on exit matrix block A21
   \param A22          : on exit matrix block A22
   */
  bool SplitMatrix2x2(Teuchos::RCP<LINALG::SparseMatrix> A, Teuchos::RCP<Epetra_Map>& A11rowmap,
      Teuchos::RCP<Epetra_Map>& A22rowmap, Teuchos::RCP<Epetra_Map>& A11domainmap,
      Teuchos::RCP<Epetra_Map>& A22domainmap, Teuchos::RCP<LINALG::SparseMatrix>& A11,
      Teuchos::RCP<LINALG::SparseMatrix>& A12, Teuchos::RCP<LINALG::SparseMatrix>& A21,
      Teuchos::RCP<LINALG::SparseMatrix>& A22);

  /** \brief Insert a diagonal row vector into a unfilled SparseMatrix
   *         on each proc without communication
   *
   *  \param mat (out) : Unfilled matrix
   *  \param diag (in) : Given diagonal (row-layout)
   *
   *  Return 0, if successful. If the given matrix is already filled, the method
   *  returns -1. In this case you should use ReplaceDiagonalValues(), instead.
   *
   *  \author hiermeier \date 03/17 */
  int InsertMyRowDiagonalIntoUnfilledMatrix(LINALG::SparseMatrix& mat, const Epetra_Vector& diag);

  /*!
   \brief split a rowmap of matrix A

   splits A->RowMap() into 2 maps, where one is given on input and the other map
   is created as complementary map. The complementary map is returned.

   \param Amap      : Map to split on input
   \param Agiven    : on entry submap that is given and part of Amap
   \return the remainder map of Amap that is not overlapping with Agiven
   */
  Teuchos::RCP<Epetra_Map> SplitMap(const Epetra_Map& Amap, const Epetra_Map& Agiven);

  /*!
   \brief merges two given Epetra_Maps

   merges input map1 and input map2, both of which have to be unique,
   but may be overlapping, to a new map and returns Teuchos::RCP to it.

   \param map1         : one map to be merged
   \param map2         : the other map to be merged
   \param allowoverlap : when set to false, an error is thrown if the result
   map is overlapping (default = true, overlap allowed)
   \return the (sorted) merged map of input maps map1 and map2
   */
  Teuchos::RCP<Epetra_Map> MergeMap(
      const Epetra_Map& map1, const Epetra_Map& map2, bool overlap = true);

  /*!
   \brief find the intersection set of two given Epetra_Maps

   Find the insection set of input map1 and input map2.

   \param map1         : first map
   \param map2         : second map
   \return the (sorted) intersection map of input maps map1 and map2
   */
  Teuchos::RCP<Epetra_Map> IntersectMap(const Epetra_Map& map1, const Epetra_Map& map2);


  /*!
   \brief merges two given Epetra_Maps

   merges input map1 and input map2 (given as Teuchos::RCP), both of which
   have to be unique, but may be overlapping, to a new map and returns
   Teuchos::RCP to it. The case that one or both input Teuchos::RCPs are null is
   detected and handled appropriately.

   \param map1         : one map to be merged
   \param map2         : the other map to be merged
   \param allowoverlap : when set to false, an error is thrown if the result
   map is overlapping (default = true, overlap allowed)
   \return the (sorted) merged map of input maps map1 and map2
   */
  Teuchos::RCP<Epetra_Map> MergeMap(const Teuchos::RCP<const Epetra_Map>& map1,
      const Teuchos::RCP<const Epetra_Map>& map2, bool overlap = true);

  /*! \brief Extract a locally overlapping sub-map from the source map
   *
   *  Find all gids in the source map which are also part of the reference map.
   *  This function works only locally on each proc. No communication is
   *  performed.
   *
   *  \param(in) src_map  source map which shall be considered
   *  \param(in) ref_map  reference map
   *  \return overlapping sub-map
   *
   *  \author hiermeier \date 03/18 */
  Teuchos::RCP<Epetra_Map> ExtractMyOverlappingSubMap(
      const Epetra_BlockMap& src_map, const Epetra_BlockMap& ref_map);

  /*!
     \brief split a vector into 2 non-overlapping pieces (Teuchos::RCP version)

     \param xmap    : map of vector to be split
     \param x       : vector to be split
     \param x1map   : map of first vector to be extracted
     \param x1      : first vector to be extracted
     \param x2map   : map of second vector to be extracted
     \param x2      : second vector to be extracted

     */
  bool SplitVector(const Epetra_Map& xmap, const Epetra_Vector& x, Teuchos::RCP<Epetra_Map>& x1map,
      Teuchos::RCP<Epetra_Vector>& x1, Teuchos::RCP<Epetra_Map>& x2map,
      Teuchos::RCP<Epetra_Vector>& x2);

  /*!
   \brief split a vector into 2 non-overlapping pieces (Teuchos::RCP version)

   \param xmap    : map of vector to be split
   \param x       : vector to be split
   \param x1map   : map of first vector to be extracted
   \param x1      : first vector to be extracted
   \param x2map   : map of second vector to be extracted
   \param x2      : second vector to be extracted

   */
  bool SplitVector(const Epetra_Map& xmap, const Epetra_Vector& x,
      Teuchos::RCP<const Epetra_Map>& x1map, Teuchos::RCP<Epetra_Vector>& x1,
      Teuchos::RCP<const Epetra_Map>& x2map, Teuchos::RCP<Epetra_Vector>& x2);

  /*! \brief compute a set of dof maps from a provided set of node maps
   *
   *  \param(in) dofset_id  considered dof set id. In the default case this value
   *                        can be set to zero.
   *  \param(in) node_maps  set of provided node maps
   *  \param(in) discret    underlying discretization
   *  \param(out) dof_maps  set of corresponding dof maps
   *
   *  \author hiermeier \date 03/18 */
  void ComputeDofMapsFromNodeMaps(const int dofset_id,
      const std::vector<Teuchos::RCP<Epetra_Map>>& node_maps,
      const DRT::DiscretizationInterface& discret, std::vector<Teuchos::RCP<Epetra_Map>>& dof_maps);

  /*! \brief Compute a dof map corresponding to the given node map
   *
   *  \param(in) dofset_id  considered dof set id. In the default case this value
   *                        can be set to zero.
   *  \param(in) node_map   provided node map
   *  \param(in) discret    underlying discretization
   *  \return corresponding dof map
   *
   *  \author hiermeier \date 03/18 */
  Teuchos::RCP<Epetra_Map> ComputeDofMapFromNodeMap(
      const int dofset_id, const Epetra_Map& node_map, const DRT::DiscretizationInterface& discret);

}  // namespace LINALG

#endif /* SRC_LINALG_LINALG_UTILS_DENSEMATRIX_MANIPULATION_H_ */
