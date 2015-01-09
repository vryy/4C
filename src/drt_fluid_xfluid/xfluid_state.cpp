/*!----------------------------------------------------------------------
\file xfluid_state.cpp
\brief State class for (in)stationary XFEM fluid problems

<pre>
Maintainer:  Raffaela Kruse and Benedikt Schott
             [kruse,schott]@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15240
</pre>

*----------------------------------------------------------------------*/

#include "xfluid_state.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_utils.H"

#include "../drt_cut/cut_cutwizard.H"

#include "../drt_xfem/xfem_dofset.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../drt_fluid/fluid_utils.H"

/*----------------------------------------------------------------------*
 |  Constructor for XFluidState                             kruse 08/14 |
 *----------------------------------------------------------------------*/
FLD::XFluidState::XFluidState(
  Teuchos::RCP<GEO::CutWizard> & wizard,
  Teuchos::RCP<XFEM::XFEMDofSet> &  dofset,
  Teuchos::RCP<const Epetra_Map> & xfluiddofrowmap) :
  xfluiddofrowmap_(xfluiddofrowmap),
  dofset_(dofset),
  wizard_(wizard)
{
  InitSystemMatrix();
  InitStateVectors();
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


/*----------------------------------------------------------------------*
 |  Initialize state vectors                               schott 12/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidState::ZeroCouplingMatricesAndRhs()
{
  // zero all coupling matrices and rhs vectors
  C_xs_->Zero();
  C_sx_->Zero();
  C_ss_->Zero();
  rhC_s_->Scale(0.0);
  rhC_s_col_->Scale(0.0);
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
  xfluiddiscret->EvaluateDirichlet(eleparams, zeros_, Teuchos::null, Teuchos::null, Teuchos::null, dbcmaps_);

  zeros_->PutScalar(0.0);

  // create vel-pres splitter
  const int numdim = DRT::Problem::Instance()->NDim();
  velpressplitter_ = Teuchos::rcp( new LINALG::MapExtractor());
  FLD::UTILS::SetupFluidSplit(*xfluiddiscret, numdim, 1, *velpressplitter_);
}


