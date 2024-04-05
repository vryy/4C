/*----------------------------------------------------------------------*/
/*! \file

\brief Specification of convection downwind numbering of a matrix

\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINALG_DOWNWINDMATRIX_HPP
#define FOUR_C_LINALG_DOWNWINDMATRIX_HPP

// Trilinos includes
#include "baci_config.hpp"

#include "baci_linalg_mapextractor.hpp"
#include "baci_linalg_sparsematrix.hpp"
#include "baci_utils_exceptions.hpp"

#include <Epetra_CrsMatrix.h>
#include <Epetra_IntVector.h>
#include <EpetraExt_Reindex_CrsMatrix.h>
#include <EpetraExt_Reindex_LinearProblem.h>
#include <EpetraExt_Reindex_MultiVector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

namespace CORE::LINALG
{
  /*!

  \brief Creates a convection downwind numbering of a matrix
  */
  class DownwindMatrix
  {
   public:
    /*!
    \brief Standard Constructor

    \param A (in) : Matrix to create downwind numbering for
    \param nv (in): number of velocity degrees of freedom per node, may not be zero
    \param np (in): number of pressure degrees of freedom per node, may be zero
    \param tau (in): threshold parameter, somewhere between 1.0 and about 4.0
    \param outlevel (in): level of output to screen (0 no output, >0 output)
    */
    explicit DownwindMatrix(Teuchos::RCP<Epetra_CrsMatrix> A, const int nv, const int np,
        const double tau, const int outlevel);

    /*!
    \brief Destructor
    */
    virtual ~DownwindMatrix() = default;
    /*!
    \brief reindex a given matrix
    */
    Teuchos::RCP<Epetra_CrsMatrix> Permute(Epetra_CrsMatrix* Aorig)
    {
      // shallow reindex the matrix
      EpetraExt::CrsMatrix_Reindex reindexer(*ndofrowmap_);
      Epetra_CrsMatrix& reA = reindexer(*Aorig);
      reA.FillComplete();
      Teuchos::RCP<Epetra_CrsMatrix> dwA =
          Teuchos::rcp(new Epetra_CrsMatrix(::Copy, DownwindRowMap(), reA.MaxNumEntries()));
      dwA->Export(reA, *sexporter_, Insert);
      if (!dwA->Filled()) dwA->FillComplete();
      return dwA;
    }

    /*!
    \brief reindex a given vector
    */
    Teuchos::RCP<Epetra_MultiVector> Permute(Epetra_MultiVector* xorig)
    {
      Epetra_MultiVector* dwx = new Epetra_Vector(DownwindRowMap(), false);
      Permute(xorig, dwx);
      return Teuchos::rcp(dwx);
    }

    /*!
    \brief reindex a given vector
    */
    void Permute(Epetra_MultiVector* xorig, Epetra_MultiVector* dwx)
    {
      // shallow reindex the vector
      EpetraExt::MultiVector_Reindex reindexer(*ndofrowmap_);
      Epetra_MultiVector& rex = reindexer(*xorig);
      dwx->Export(rex, *sexporter_, Insert);
      return;
    }

    /*!
    \brief undo reindexing of a given vector
    */
    void InvPermute(Epetra_MultiVector* xdw, Epetra_MultiVector* xorig)
    {
      EpetraExt::MultiVector_Reindex reindex(*ndofrowmap_);
      Epetra_MultiVector& rex = reindex(*xorig);
      rex.Export(*xdw, *rexporter_, Insert);
      return;
    }

    /*!
    \brief Get downwinded row map of problem
    */
    inline const Epetra_Map& DownwindRowMap() { return *sndofrowmap_; }

   private:
    // don't want copy-ctor and = operator
    DownwindMatrix(DownwindMatrix& old);
    DownwindMatrix operator=(const DownwindMatrix& old);

    /// setup phase of downwinding
    void Setup(const Epetra_CrsMatrix& A);

    /// sets an index in the Bey&Wittum method
    void SetF(
        const int i, int& nf, Epetra_IntVector& index, const Epetra_CrsMatrix& graph, int rec);

    /// sets an index in the Hackbusch method
    void SetL(
        const int i, int& nl, Epetra_IntVector& index, const Epetra_CrsMatrix& graph, int rec);

    /// Do downwinding according to Bey & Wittum
    void DownwindBeyWittum(const Epetra_CrsMatrix& nnodegraph, Epetra_IntVector& index,
        const Epetra_IntVector& oninflow);

    /// Do downwinding according to Hackbusch
    void DownwindHackbusch(const Epetra_CrsMatrix& nnodegraph, Epetra_IntVector& index,
        const Epetra_IntVector& oninflow);

    /// test whether k is successor of i
    inline bool IsSuccessor(const int k, const int gi, const Epetra_CrsMatrix& graph)
    {
      int numentries;
      int* indices;
      double* values;
      graph.ExtractMyRowView(k, numentries, values, indices);
      for (int j = 0; j < numentries; ++j)
        if (graph.ColMap().GID(indices[j]) == gi) return true;
      return false;
    }

    const int outlevel_;  // level of output
    const int nv_;        // number velocity dofs
    const int np_;        // number pressure dofs
    const int bs_;        // total nodal block size (nv_ + np_)
    const double tau_;    // strong connections cutoff ratio

    Teuchos::RCP<Epetra_Map> ndofrowmap_;  // result new row map of system
    Teuchos::RCP<Epetra_Map>
        sndofrowmap_;  // result new row map of system with dofs sorted in ascending order

    Teuchos::RCP<Epetra_Export> sexporter_;  // export from ndofrowmap_ to sndofrowmap_
    Teuchos::RCP<Epetra_Export> rexporter_;  // export the other direction

  };  // class  DownwindMatrix
}  // namespace CORE::LINALG


BACI_NAMESPACE_CLOSE

#endif  // LINALG_DOWNWINDMATRIX_H
