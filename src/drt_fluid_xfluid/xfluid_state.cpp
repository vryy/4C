/*!----------------------------------------------------------------------
\file xfluid_state.cpp
\brief State class for (in)stationary XFEM fluid problems

Attention:
These classes are still prototypes and have to be completed.

<pre>
Maintainer:  Raffaela Kruse
             kruse@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15240
</pre>

*----------------------------------------------------------------------*/

#include "xfluid_state.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_krylov_projector.H"

#include "../drt_xfem/xfem_fluidwizard.H"
#include "../drt_xfem/xfem_fluiddofset.H"
#include "../drt_xfem/xfem_neumann.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../drt_fluid/fluid_utils.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"


/*----------------------------------------------------------------------*
 |  Constructor for XFluidState                             kruse 08/14 |
 *----------------------------------------------------------------------*/
FLD::XFluidState::XFluidState(
  Teuchos::RCP<const Epetra_Map> & xfluiddofrowmap) :
  xfluiddofrowmap_(xfluiddofrowmap)
{
  InitSystemMatrix();
  InitStateVectors();
}

/*----------------------------------------------------------------------*
 |  Constructor for XFluidState                             kruse 08/14 |
 *----------------------------------------------------------------------*/
FLD::XFluidState::XFluidState(
  Teuchos::RCP<const Epetra_Map> & xfluiddofrowmap,
  const Epetra_Map &               slavedofrowmap,
  const Epetra_Map &               slavedofcolmap) :
  xfluiddofrowmap_(xfluiddofrowmap)
{
  InitSystemMatrix();
  InitStateVectors();
  InitCouplingMatricesAndRhs(slavedofrowmap,slavedofcolmap);
}

/*----------------------------------------------------------------------*
 |  Initialize (coupling) matrices & rhs-vectors
 |                                                         kruse 08/14
 *----------------------------------------------------------------------*/
void FLD::XFluidState::InitSystemMatrix()
{
  // create an EpetraFECrs matrix that does communication for non-local rows and columns
  // * this enables to do the evaluate loop over just row elements instead of col elements
  // * time consuming assemble for cut elements is done only once on a unique row processor
  // REMARK: call the SparseMatrix: * explicitdirichlet = true (is used in ApplyDirichlet, false uses Epetra memory based operations
  //                                                            that are not compatible with the FE-Matrix)
  //                                * savegraph = true/false: To save the graph (pattern for non-zero entries) leads to a speedup in the assembly of the matrix
  //                                    for subsequent assemblies.
  //                                    However, do not use the savegraph-option when the matrix graph can change during the usage of
  //                                    this object of SparseMatrix or use the Reset()-function instead of Zero().
  //                                    For XFEM-problems, the matrix graph changes between timesteps, however,
  //                                    then a new state class and sparsematrix is created, otherwise a Reset()-function has to be called instead
  //                                    of the Zero()-function. We are using the save-graph option.
  // * the estimate of the number of nonzero entries is adapted to hex8 elements with 8 adjacent elements around a node
  //   + edge-based couplings component-wise v_x->u_x, v_y->u_y, v_z->u_z, q->p
  //   number of non-zeros (for hex8 elements): 108+54 = 162
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*xfluiddofrowmap_,162,true,true,LINALG::SparseMatrix::FE_MATRIX));
}

void FLD::XFluidState::InitCouplingMatricesAndRhs(const Epetra_Map & slavedofrowmap,
                                                  const Epetra_Map & slavedofcolmap)
{
  // savegraph flag set to true, as there is no change in the matrix graph expected for the lifetime
  // of this state container
  C_xs_  = Teuchos::rcp(new LINALG::SparseMatrix(*xfluiddofrowmap_,0,true,true,LINALG::SparseMatrix::FE_MATRIX));
  C_sx_  = Teuchos::rcp(new LINALG::SparseMatrix(slavedofrowmap,0,true,true,LINALG::SparseMatrix::FE_MATRIX));
  C_ss_  = Teuchos::rcp(new LINALG::SparseMatrix(slavedofrowmap,0,true,true,LINALG::SparseMatrix::FE_MATRIX));
  rhC_s_ = LINALG::CreateVector(slavedofrowmap,true);
  rhC_s_col_= LINALG::CreateVector(slavedofcolmap,true);
}

/*----------------------------------------------------------------------*
 |  Initialize state vectors                                kruse 08/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidState::InitStateVectors()
{
  // Vectors passed to the element
  // -----------------------------
  // velocity/pressure at time n+1, n and n-1
  velnp_ = LINALG::CreateVector(*xfluiddofrowmap_,true);
  veln_  = LINALG::CreateVector(*xfluiddofrowmap_,true);
  velnm_ = LINALG::CreateVector(*xfluiddofrowmap_,true);

  // velocity/pressure at time n+alpha_F
  velaf_ = LINALG::CreateVector(*xfluiddofrowmap_,true);


  // acceleration/(scalar time derivative) at time n+1 and n
  accnp_ = LINALG::CreateVector(*xfluiddofrowmap_,true);
  accn_  = LINALG::CreateVector(*xfluiddofrowmap_,true);

  // acceleration/(scalar time derivative) at time n+alpha_M/(n+alpha_M/n)
  accam_ = LINALG::CreateVector(*xfluiddofrowmap_,true);

  // scalar at time n+alpha_F/n+1 and n+alpha_M/n
  // (only required for low-Mach-number case)
  // ... this is a dummy to avoid errors
  scaaf_ = LINALG::CreateVector(*xfluiddofrowmap_,true);
  scaam_ = LINALG::CreateVector(*xfluiddofrowmap_,true);

  // history vector
  hist_ = LINALG::CreateVector(*xfluiddofrowmap_,true);

  // the vector containing body and surface forces
  neumann_loads_= LINALG::CreateVector(*xfluiddofrowmap_,true);

  // rhs: standard (stabilized) residual vector (rhs for the incremental form)
  residual_     = LINALG::CreateVector(*xfluiddofrowmap_,true);
  trueresidual_ = LINALG::CreateVector(*xfluiddofrowmap_,true);

  // right hand side vector for linearized solution;
  rhs_ = LINALG::CreateVector(*xfluiddofrowmap_,true);

  // nonlinear iteration increment vector
  incvel_ = LINALG::CreateVector(*xfluiddofrowmap_,true);

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_   = LINALG::CreateVector(*xfluiddofrowmap_,true);

  // Todo: ALE-displacements and grid-velocities for the case of an
  // XFEM-ALE-fluid are not initialized, as this option is not
  // implemented yet!

  //if(alefluid_)...
  dserror("ALE displacments!!!");
}

/*----------------------------------------------------------------------*
 |  Set dirichlet- and velocity/pressure-map extractor      kruse 08/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidState::SetupMapExtractors(const Teuchos::RCP<DRT::Discretization> & xfluiddiscret,
                                          const double & time)
{
  // create dirichlet map extractor
  Teuchos::ParameterList eleparams;
  // other parameters needed by the elements
  eleparams.set("total time",time);
  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  xfluiddiscret->EvaluateDirichlet(eleparams, zeros_, Teuchos::null, Teuchos::null,
                                   Teuchos::null, dbcmaps_);

  zeros_->PutScalar(0.0);

  // create vel-pres splitter
  const int numdim = DRT::Problem::Instance()->NDim();
  velpressplitter_ = Teuchos::rcp( new LINALG::MapExtractor());
  FLD::UTILS::SetupFluidSplit(*xfluiddiscret, numdim, 1, *velpressplitter_);
}

/// setup map extractors for dirichlet maps & vel/pres-splitter
void FLD::XFluidState::SetupKSPMapExtractor(
  const Teuchos::RCP<DRT::Discretization> & xfluiddiscret
  )
{
  // create map of nodes involved in Krylov projection
  kspsplitter_ = Teuchos::rcp(new FLD::UTILS::KSPMapExtractor);
  kspsplitter_->Setup(*(xfluiddiscret));
}


/*----------------------------------------------------------------------*
 |  Set xfem fluid wizard                                   kruse 08/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidStateCreator::SetWizard(
  const Teuchos::RCP<DRT::Discretization> & discret)
{
  wizard_ = Teuchos::rcp( new XFEM::FluidWizardLevelSet(*discret) );

}
/*----------------------------------------------------------------------*
 |  Set xfem fluid wizard                                   kruse 08/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidStateCreator::SetWizard(
  const Teuchos::RCP<DRT::Discretization> & discret,
  const Teuchos::RCP<DRT::Discretization> & boundarydiscret)
{
  // Initialize the mesh wizard
  wizard_ = Teuchos::rcp( new XFEM::FluidWizardMesh(*discret, *boundarydiscret) );
}

/*----------------------------------------------------------------------*
 |  Perform the cut and fill state container                kruse 08/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<FLD::XFluidState> FLD::XFluidStateCreator::Create(
  const Teuchos::RCP<DRT::Discretization>&  discret,
  const Teuchos::RCP<DRT::Discretization>&  boundarydiscret,
  const Epetra_Vector &                     idispcol,
  Teuchos::ParameterList &                  solver_params,
  const int                                 step,
  const double &                            time,
  bool                                      coupling
)
{
  if (wizard_ == Teuchos::null)
    dserror("Uninitialized mesh wizard. Cannot create state.");

  //if (alefluid_)...
    dserror("FLD::XFluidStateCreator::Create: Add Ale displacements to Cut()!!! - at the moment Teuchos::null!");

  //--------------------------------------------------------------------------------------
  // the XFEM::FluidWizardMesh is created based on the xfluid-discretization and the boundary discretization
  // the FluidWizardMesh creates also a cut-object of type GEO::CutWizardMesh which performs the "CUT"
  wizard_->Cut( false,                         // include_inner
                idispcol,                      // interface displacements
                VolumeCellGaussPointBy_,       // how to create volume cell Gauss points?
                BoundCellGaussPointBy_,        // how to create boundary cell Gauss points?
                true,                          // parallel cut framework
                gmsh_cut_out_,                 // gmsh output for cut library
                Teuchos::null,                 // no ale displacements
                true                           // find point positions
                );

  //--------------------------------------------------------------------------------------
  // set the new dofset after cut
  const int maxNumMyReservedDofsperNode = (maxnumdofsets_)*4;
  dofset_ = wizard_->DofSet(maxNumMyReservedDofsperNode);

  const int restart = DRT::Problem::Instance()->Restart();

  if ((step < 1) or restart)
    minnumdofsets_ = discret->DofRowMap()->MinAllGID();

  dofset_->MinGID(minnumdofsets_); // set the minimal GID of xfem dis
  discret->ReplaceDofSet( dofset_, true );

  discret->FillComplete();

  //print all dofsets
  discret->GetDofSetProxy()->PrintAllDofsets(discret->Comm());

  //--------------------------------------------------------------------------------------
  // recompute nullspace based on new number of dofs per node
  // REMARK: this has to be done after replacing the discret_' dofset (via xfluid.discret_->ReplaceDofSet)
  discret->ComputeNullSpaceIfNecessary(solver_params,true);

  Teuchos::RCP<const Epetra_Map> xfluiddofrowmap = Teuchos::rcp(
      new Epetra_Map(*discret->DofRowMap()));

  Teuchos::RCP<FLD::XFluidState> state;
  if (coupling)
  {
    state = Teuchos::rcp(
      new FLD::XFluidState(xfluiddofrowmap,*boundarydiscret->DofRowMap(),*boundarydiscret->DofColMap()));
  }
  else
  {
    state = Teuchos::rcp(
          new FLD::XFluidState(xfluiddofrowmap));
  }

  state->SetupMapExtractors(discret,time);

  return state;
}

///*----------------------------------------------------------------------*
// |  Constructor for XFluidMeshState                         kruse 08/14 |
// *----------------------------------------------------------------------*/
//FLD::XFluidMeshState::XFluidMeshState(Teuchos::RCP<const Epetra_Map> &  xfluiddofrowmap,
//                                      const Epetra_Map &                slavedofrowmap,
//                                      const Epetra_Map &                slavedofcolmap
//                                      ) : XFluidState(xfluiddofrowmap)
//{
//  InitCouplingMatricesAndRhs(slavedofrowmap,slavedofcolmap);
//}
//
///*----------------------------------------------------------------------*
// |  Initialize state vectors                                kruse 08/14 |
// *----------------------------------------------------------------------*/
//void FLD::XFluidMeshState::InitStateVectors()
//{
//  FLD::XFluidState::InitStateVectors();
//}
//
///*----------------------------------------------------------------------*
// |  Initialize state vectors                                kruse 08/14 |
// *----------------------------------------------------------------------*/
//void FLD::XFluidLevelSetState::InitStateVectors()
//{
//  FLD::XFluidState::InitStateVectors();
//}
