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
#include "../drt_lib/linalg_blocksparsematrix.H"

namespace FLD
{

/*!
 * \brief namespace for fluid related stuff
 */
namespace UTILS
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

  \param dis : (in) Discretization
  \param ndim : (in) dimensions of problem (2 or 3)
  \param extractor : (out) ready made map splitter

  \author u.kue
  \date 02/08
 */
void SetupFluidSplit(const DRT::Discretization& dis,
                     const DRT::DofSet& dofset,
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

/*!

\brief calculate lift&drag forces and angular momenta

Lift and drag forces are based upon the right hand side true-residual entities
of the corresponding nodes. The contribution of the end node of a line is entirely
added to a present L&D force.

Idea of this routine:

create

map< label, set<DRT::Node*> >

which is a set of nodes to each L&D Id
nodal forces of all the nodes within one set are added to one L&D force

Notice: Angular moments obtained from lift&drag forces currently refere to the
        initial configuration, i.e. are built with the coordinates X of a particular
        node irrespective of its current position.

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>


  \param dis          (in ) the discretisation (node distribution, conditions)
  \param force        (in ) the true residual (the nodal forces)
  \param params       (in ) parameterlist including the lifdrag flag
  \param liftdragvals (out) the computed values for lift and drag in an array 
  

  \author chfoe 
  \date 11/07
 */
void LiftDrag(
  const DRT::Discretization&     dis         , //! fluid discretization
  const Epetra_Vector&           force       , //! vector of nodalforces
  const ParameterList&           params      , //! parameters 
  RCP<map<int,vector<double> > > liftdragvals  //! computed ld values
  );

/*!
  \brief integrate mass flow over surfaces

  for each condition Id, compute the flow through the surface


  \return a map, where for each condition Id, we get the flowrate.
          positive and negative signs indicate net inflow and outflow


  \author Axel Gerstenberger
  \date 10/08
 */
std::map<int,double> ComputeSurfaceFlowrates(
    DRT::Discretization&       dis  ,      ///< the discretisation (node distribution, conditions)
    const RCP<Epetra_Vector>   velnp       ///< solution vector with velocities and pressure 
                                           ///< (only velocities are used)
    );

} // namespace UTILS
} // namespace FLD

#endif  // #ifndef FLUID_UTILS_H
#endif  // #ifdef CCADISCRET
