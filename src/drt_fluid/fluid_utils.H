/*!----------------------------------------------------------------------
\file fluid_utils.H
\brief utility functions for fluid problems

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef FLUID_UTILS_H
#define FLUID_UTILS_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "../drt_lib/drt_discret.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_lib/linalg_mapextractor.H"
#include "../drt_lib/linalg_systemmatrix.H"

/*!
 * \brief namespace for fluid related stuff
 */
namespace FLUID_UTILS
{

  /// velocity-pressure block matrix split strategy
  class VelPressSplitStrategy
  {
  public:

    /// construct with a block matrix base
    explicit VelPressSplitStrategy(LINALG::BlockSparseMatrixBase& mat)
      : mat_(mat), numdim_(-1) {}

    /// find row block to a given row gid
    int RowBlock(int lrow, int rgid)
    {
      if ((lrow%(numdim_+1))<numdim_)
        return 0;
      return 1;
    }

    /// find column block to a given column gid
    int ColBlock(int rblock, int lcol, int cgid)
    {
      if ((lcol%(numdim_+1))<numdim_)
        return 0;
      return 1;
    }

    /// assemble into the given block
    void Assemble(int eid,
                  int myrank,
                  const Epetra_SerialDenseMatrix& Aele,
                  const std::vector<int>& lmrow,
                  const std::vector<int>& lmrowowner,
                  const std::vector<int>& lmcol)
    {
      const int lrowdim = (int)lmrow.size();
      const int lcoldim = (int)lmcol.size();

      // loop rows of local matrix
      for (int lrow=0; lrow<lrowdim; ++lrow)
      {

        // check ownership of row
        if (lmrowowner[lrow]!=myrank)
          continue;

        int rgid = lmrow[lrow];
        int rblock = RowBlock(lrow,rgid);

        for (int lcol=0; lcol<lcoldim; ++lcol)
        {
          double val = Aele(lrow,lcol);
          int cgid = lmcol[lcol];
          int cblock = ColBlock(rblock,lcol,cgid);

          LINALG::SparseMatrix& matrix = mat_.Matrix(rblock,cblock);
          matrix.Assemble(val,rgid,cgid);
        }
      }
    }

    /// assemble into the given block
    void Assemble(double val, int rgid, int cgid)
    {
      int rblock = RowBlock(0,rgid);
      int cblock = ColBlock(rblock,0,cgid);
      LINALG::SparseMatrix& matrix = mat_.Matrix(rblock,cblock);
      matrix.Assemble(val,rgid,cgid);
    }

    /// assemble the remaining ghost entries
    void Complete() {}

    /// set number of velocity dofs
    void SetNumdim(int numdim) { numdim_ = numdim; }

  private:

    /// my block matrix base
    LINALG::BlockSparseMatrixBase& mat_;

    /// number of velocity dofs
    int numdim_;
  };


  /// (FSI) interface block matrix split strategy
  class InterfaceSplitStrategy : public LINALG::DefaultBlockMatrixStrategy
  {
  public:
    explicit InterfaceSplitStrategy(LINALG::BlockSparseMatrixBase& mat)
      : LINALG::DefaultBlockMatrixStrategy(mat)
    {}

    /// assemble into the given block
    void Assemble(int eid,
                  int myrank,
                  const Epetra_SerialDenseMatrix& Aele,
                  const std::vector<int>& lmrow,
                  const std::vector<int>& lmrowowner,
                  const std::vector<int>& lmcol)
    {
      if (condelements_->find(eid)!=condelements_->end())
      {
        // if we have an element with conditioned nodes, we have to do the
        // default assembling
        LINALG::DefaultBlockMatrixStrategy::Assemble(eid,myrank,Aele,lmrow,lmrowowner,lmcol);
      }
      else
      {
        // if there are no conditioned nodes we can simply assemble to the
        // internal matrix
        LINALG::SparseMatrix& matrix = Mat().Matrix(0,0);
        matrix.Assemble(eid,Aele,lmrow,lmrowowner,lmcol);
      }
    }

    void Assemble(double val, int rgid, int cgid)
    {
      // forward single value assembling
      LINALG::DefaultBlockMatrixStrategy::Assemble(val,rgid,cgid);
    }

    void SetCondElements(Teuchos::RCP<std::set<int> > condelements)
    { condelements_ = condelements; }

  private:
    Teuchos::RCP<std::set<int> > condelements_;
  };

/// Split dof row map of discretization into velocity and pressure maps
/*!

  We assume that each node possesses ndim velocity dofs and (optionally)
  one pressure dof. The dof row map is thus split in two.

  The velocity dofs are assigned to map 0 (the other map) and the pressure
  dofs are assigned to map 1 (the condition map).

  \param dis : (in) Discretization
  \param ndim : (in) dimensions of problem (2 or 3)
  \param extractor : (out) ready made map splitter

  \author u.kue
  \date 02/08
 */
void SetupFluidSplit(const DRT::Discretization& dis,
                     int ndim,
                     LINALG::MapExtractor& extractor);

/// Split dof row map of discretization into velocity and pressure maps
/*!

  We assume that each node possesses ndim velocity dofs and (optionally)
  one pressure dof. The dof row map is thus split in two.

  The velocity dofs are assigned to map 0 (the other map) and the pressure
  dofs are assigned to map 1 (the condition map).

  \author a.ger
  \date 03/08
 */
void SetupXFluidSplit(
        const DRT::Discretization&   dis,       ///< fluid discretization
        const RCP<XFEM::DofManager>  dofman,    ///< information about the dofs per node
        LINALG::MapExtractor&        extractor  ///< ready made map splitter
        );

}

#endif  // #ifndef FLUID_UTILS_H
#endif  // #ifdef CCADISCRET
