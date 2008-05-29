
#ifdef CCADISCRET

#ifndef FSI_MATRIXTRANSFORM_H
#define FSI_MATRIXTRANSFORM_H

#include "../drt_lib/linalg_systemmatrix.H"
#include "../drt_adapter/adapter_coupling.H"

#include <map>

#include <Teuchos_RCP.hpp>
#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include <Epetra_CrsMatrix.h>


namespace FSI
{

  /// base
  /*!

    This class encapsulates a very specific operation.

    \author u.kue
    \date 05/08
   */
  class MatrixTransform
  {
  protected:

    void InsertValues(Teuchos::RCP<Epetra_CrsMatrix> edst,
                      const Epetra_Map& dstrowmap,
                      const Epetra_Map& dstmap,
                      int row,
                      int NumEntries,
                      double *Values,
                      int *Indices);
  };

  /// communicate matrix to new row map
  /*!

    This class encapsulates a very specific operation.

    \author u.kue
    \date 05/08
   */
  class MatrixRowTransform : public MatrixTransform
  {
  public:

    /// transformation operation
    bool operator()(const LINALG::SparseMatrix& src,
                    double scale,
                    const ADAPTER::Coupling::Converter& converter,
                    LINALG::SparseMatrix& dst);

    /// redistribute source matrix according to permutation map
    Teuchos::RCP<Epetra_CrsMatrix> Redistribute(const LINALG::SparseMatrix& src,
                                                const Epetra_Map& permsrcmap);

    /// copy values from source to destination matrix
    void MatrixInsert(Teuchos::RCP<Epetra_CrsMatrix> permsrc,
                      const Epetra_Map& dstmap,
                      Teuchos::RCP<Epetra_CrsMatrix> edst);

  private:
    Teuchos::RCP<Epetra_Export> exporter_;
  };


  /// communicate matrix to new col map
  /*!

    This class encapsulates a very specific operation.

    \author u.kue
    \date 05/08
   */
  class MatrixColTransform : public MatrixTransform
  {
  public:

    /// construct
    MatrixColTransform() : havegidmap_(false) {}

    /// transformation operation
    bool operator()(const LINALG::BlockSparseMatrixBase& fullsrc,
                    const LINALG::SparseMatrix& src,
                    double scale,
                    const ADAPTER::Coupling::Converter& converter,
                    LINALG::SparseMatrix& dst,
                    bool exactmatch=true);

  private:
    std::map<int,int> gidmap_;
    bool havegidmap_;
  };
}

#endif
#endif
