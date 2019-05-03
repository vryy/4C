/*!----------------------------------------------------------------------
\file fsi_matrixtransform.H

\brief Utilities for matrix transformations in FSI

\maintainer Martin Kronbichler

\level 1
*/
/*---------------------------------------------------------------------*/

#ifndef FSI_MATRIXTRANSFORM_H
#define FSI_MATRIXTRANSFORM_H

#include <map>

#include <Teuchos_RCP.hpp>
#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include <Epetra_CrsMatrix.h>

// forward declarations
namespace ADAPTER
{
  class Coupling;
  class CouplingConverter;
}  // namespace ADAPTER

namespace LINALG
{
  class SparseMatrix;
}

namespace FSI
{
  namespace UTILS
  {
    /// extract submatrix of the src map and transform it to a new col map
    /*!

    Overlapping monolithic FSI adds matrices from different fields at the
    interface. These matrices belong to different row maps. Thus adding them
    requires moving one of them to a new row map. The relations between these
    maps are managed by ADAPTER::Coupling objects. In a parallel setting there
    is a master and a slave side (in case of matrix transformations we use
    source and destination abstraction via ADAPTER::CouplingConverter). The
    parallel distribution of both is arbitrary. And in addition there are
    permuted master and slave maps, that match the respective other side. So
    the row map transformation requires a parallel redistribution followed by
    a row map exchange.

    This operator does not utilize Epetra_CrsMatrix::ReplaceRowMap() but
    copies the temporary matrix. This is required both to preserve the
    internal Epetra_CrsMatrix from the destination LINALG::SparseMatrix and
    because the destination matrix row map may be much larger than the source
    matrix row map.

    The operator is meant to be usable on its own and operate on both row and
    column transformations (if the respective converter is given).

    An additional feature of this class is that it can assign matrix blocks
    from one field to block matrix slots on another field. As opposed to
    MatrixColTransform, this method extracts a logical block from the input
    matrix without any split call.

    \note All matrix transformation operators work with filled and unfilled
    destination matrices. The source matrix is never changed. The destination
    matrix is not reallocated and its filled state is not explicitly
    changed. There is a LINALG::SparseMatrix::Zero() call if addmatrix==false
    and this can reset the filled state if the matrix graph is not preserved
    by the LINALG::SparseMatrix object.

    \warning This class encapsulates a very specific operation needed for
    monolithic FSI block matrix assemble.

    \sa MatrixColTransform, MatrixRowTransform, MatrixRowColTransform
    \author kronbichler
    \date 11/15
     */
    class MatrixLogicalSplitAndTransform
    {
     public:
      /// construct
      MatrixLogicalSplitAndTransform() : havegidmap_(false) {}

      /// transformation operation
      /*!
        The call operator to be used for a matrix data copy between \p src and
        \p dst matrices.

        \param src            (i) source matrix
        \param logical_range_map  (i) sub-map that defines a logical block matrix row within src
        \param logical_domain_map (i) sub-map that defines a logical block matrix column within src
        \param scale          (i) scaling factor to be applied
        \param row_converter  (i) src to dst abstraction on ADAPTER::Coupling, if NULL: no row
        transform is done \param col_converter  (i) src to dst abstraction on ADAPTER::Coupling, if
        NULL: no column transform is done \param dst          (i/o) destination matrix \param
        exactmatch     (i) do not drop any source values if true \param addmatrix      (i) remove
        current dst values if false
       */
      bool operator()(const LINALG::SparseMatrix& src, const Epetra_Map& logical_range_map,
          const Epetra_Map& logical_domain_map, const double scale,
          const ::ADAPTER::CouplingConverter* row_converter,
          const ::ADAPTER::CouplingConverter* col_converter, LINALG::SparseMatrix& dst,
          bool exactmatch = true, bool addmatrix = false);

     private:
      /// setup column map matching between source and destination gids
      /*!
        Internal method.
       */
      void SetupGidMap(const Epetra_Map& rowmap, const Epetra_Map& colmap,
          const ::ADAPTER::CouplingConverter* converter, const Epetra_Comm& comm);

      /// copy values from source to destination matrix
      /*!
        Internal method.
       */
      void InternalAdd(Teuchos::RCP<Epetra_CrsMatrix> esrc, const Epetra_Map& logical_range_map,
          const Epetra_Map& logical_domain_map, const Epetra_Map& matching_dst_rows,
          Teuchos::RCP<Epetra_CrsMatrix> edst, bool exactmatch, double scale);

      /// fast method that adds into filled matrices
      /*!
        Internal method called by InternalAdd.
       */
      void AddIntoFilled(Teuchos::RCP<Epetra_CrsMatrix> esrc, const Epetra_Map& logical_range_map,
          const Epetra_Map& logical_domain_map, const Epetra_Vector& selector,
          const Epetra_Map& matching_dst_rows, Teuchos::RCP<Epetra_CrsMatrix> edst, bool exactmatch,
          double scale);

      /// slow method that adds into unfilled matrices
      /*!
        Internal method called by InternalAdd.
       */
      void AddIntoUnfilled(Teuchos::RCP<Epetra_CrsMatrix> esrc, const Epetra_Map& logical_range_map,
          const Epetra_Map& logical_domain_map, const Epetra_Vector& selector,
          const Epetra_Map& matching_dst_rows, Teuchos::RCP<Epetra_CrsMatrix> edst, bool exactmatch,
          double scale);

      /// source and destination gid matching
      std::map<int, int> gidmap_;

      /// setup done flag
      bool havegidmap_;

      /// localized version of gidmap_
      std::vector<int> lidvector_;

      /// exporter to communicate matrix to new row map
      Teuchos::RCP<Epetra_Export> exporter_;
    };



    /// communicate matrix to new row map
    /*!

      Overlapping monolithic FSI adds matrices from different fields at the
      interface. These matrices belong to different row maps. Thus adding them
      requires moving one of them to a new row map. The relations between these
      maps are managed by ADAPTER::Coupling objects. In a parallel setting there
      is a master and a slave side (in case of matrix transformations we use
      source and destination abstraction via ADAPTER::CouplingConverter). The
      parallel distribution of both is arbitrary. And in addition there are
      permuted master and slave maps, that match the respective other side. So
      the row map transformation requires a parallel redistribution followed by
      a row map exchange.

      This operator does not utilize Epetra_CrsMatrix::ReplaceRowMap() but
      copies the temporary matrix. This is required both to preserve the
      internal Epetra_CrsMatrix from the destination LINALG::SparseMatrix and
      because the destination matrix row map may be much larger than the source
      matrix row map.

      The operator is meant to be usable on its own or as part of the composed
      MatrixRowColTransform operator.

      \note The implementation is done by MatrixLogicalSplitAndTransform

      \note All matrix transformation operators work with filled and unfilled
      destination matrices. The source matrix is never changed. The destination
      matrix is not reallocated and its filled state is not explicitly
      changed. There is a LINALG::SparseMatrix::Zero() call if addmatrix==false
      and this can reset the filled state if the matrix graph is not preserved
      by the LINALG::SparseMatrix object.

      \warning This class encapsulates a very specific operation needed for
      monolithic FSI block matrix assemble.

      \sa MatrixLogicalSplitAndTransform, MatrixColTransform, MatrixRowColTransform
      \author u.kue
      \date 05/08
     */
    class MatrixRowTransform
    {
     public:
      /// transformation operation
      /*!
        The call operator to be used for a matrix data copy between \p src and
        \p dst matrices.

        \param src       (i) source matrix
        \param scale     (i) scaling factor to be applied
        \param converter (i) src to dst abstraction on ADAPTER::Coupling
        \param dst     (i/o) destination matrix
        \param addmatrix (i) remove current dst values if false
       */
      bool operator()(const LINALG::SparseMatrix& src, double scale,
          const ::ADAPTER::CouplingConverter& converter, LINALG::SparseMatrix& dst,
          bool addmatrix = false);

     private:
      /// object that does the actual work
      MatrixLogicalSplitAndTransform transformer;
    };


    /// communicate matrix to new col map
    /*!

      Overlapping monolithic FSI needs to assign matrix blocks from one field to
      block matrix slots belonging to another field. For some matrix blocks the
      row map stays the same but the column map changes.

      A special point here is that the source matrix column might include more
      values than the respective destination matrix column. This happens with
      fluid matrices that include velocity and pressure values whereas the
      corresponding ale matrices just contain displacement values. In such a
      case it is possible to advice the transformation to drop the pressure
      values (\p exactmatch=false ). By default additional values raise a
      \p dserror .

      \note The implementation is done by MatrixLogicalSplitAndTransform

      \note All matrix transformation operators work with filled and unfilled
      destination matrices. The source matrix is never changed. The destination
      matrix is not reallocated and its filled state is not explicitly
      changed. There is a LINALG::SparseMatrix::Zero() call if addmatrix==false
      and this can reset the filled state if the matrix graph is not preserved
      by the LINALG::SparseMatrix object.

      \warning This class encapsulates a very specific operation needed for monolithic
      FSI block matrix assemble.

      \sa MatrixLogicalSplitAndTransform, MatrixRowTransform, MatrixRowColTransform
      \author u.kue
      \date 05/08
     */
    class MatrixColTransform
    {
     public:
      /// transformation operation
      /*!
        The call operator to be used for a matrix data copy between \p src and
        \p dst matrices.

        \param rowmap    (i) row map of full source block matrix
        \param colmap    (i) col map of full source block matrix
        \param src       (i) source matrix
        \param scale     (i) scaling factor to be applied
        \param converter (i) src to dst abstraction on ADAPTER::Coupling
        \param dst     (i/o) destination matrix
        \param exactmatch (i) do not drop any source values if true
        \param addmatrix (i) remove current dst values if false
       */
      bool operator()(const Epetra_Map& rowmap, const Epetra_Map& colmap,
          const LINALG::SparseMatrix& src, double scale,
          const ::ADAPTER::CouplingConverter& converter, LINALG::SparseMatrix& dst,
          bool exactmatch = true, bool addmatrix = false);

     private:
      /// object that does the actual work
      MatrixLogicalSplitAndTransform transformer;
    };


    /// communicate matrix to new row map and col map
    /*!

      A combined row and column map exchange between source and destination
      matrix.

      \note The implementation is done by MatrixLogicalSplitAndTransform

      \note All matrix transformation operators work with filled and unfilled
      destination matrices. The source matrix is never changed. The destination
      matrix is not reallocated and its filled state is not explicitly
      changed. There is a LINALG::SparseMatrix::Zero() call if addmatrix==false
      and this can reset the filled state if the matrix graph is not preserved
      by the LINALG::SparseMatrix object.

      \warning This class encapsulates a very specific operation needed for monolithic
      FSI block matrix assemble.

      \sa MatrixLogicalSplitAndTransform, MatrixRowTransform, MatrixColTransform
      \author u.kue
      \date 05/08
     */
    class MatrixRowColTransform
    {
     public:
      /// transformation operation
      /*!
        The call operator to be used for a matrix data copy between \p src and
        \p dst matrices.

        \param src          (i) source matrix
        \param scale        (i) scaling factor to be applied
        \param rowconverter (i) src to dst abstraction on ADAPTER::Coupling
        \param colconverter (i) src to dst abstraction on ADAPTER::Coupling
        \param dst        (i/o) destination matrix
        \param exactmatch   (i) do not drop any source values if true
        \param addmatrix    (i) remove current dst values if false
       */
      bool operator()(const LINALG::SparseMatrix& src, double scale,
          const ::ADAPTER::CouplingConverter& rowconverter,
          const ::ADAPTER::CouplingConverter& colconverter, LINALG::SparseMatrix& dst,
          bool exactmatch = true, bool addmatrix = false);

     private:
      /// object that does the actual work
      MatrixLogicalSplitAndTransform transformer;
    };

  }  // namespace UTILS

}  // namespace FSI

#endif
