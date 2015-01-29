/*!----------------------------------------------------------------------
\file xfluidfluidnew.cpp
\brief Control routine for fluid-fluid (in)stationary solvers with XFEM

ATTENTION! Class is still a prototype. Does not provide full (and correct)
functionality yet!

<pre>
Maintainer:  Raffaela Kruse
             kruse@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
*----------------------------------------------------------------------*/

#include "xfluidfluidnew.H"
#include "xfluid_state_creator.H"

#include "xfluidfluidresulttest.H"

#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid/fluid_utils.H"

#include "../drt_lib/drt_discret_faces.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_dofset_transparent_independent.H"
#include "../drt_lib/drt_utils_parallel.H"

#include "../drt_io/io_pstream.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

#include "../drt_fluid_ele/fluid_ele_parameter_xfem.H"

#include "../drt_xfem/xfem_edgestab.H"
#include "../drt_xfem/xfem_dofset.H"
#include "../drt_xfem/xfluidfluid_timeInt.H"


FLD::XFluidFluidNew::XFluidFluidNew(
  const Teuchos::RCP<FLD::FluidImplicitTimeInt> & embedded_fluid,     ///< embedded fluid
  const Teuchos::RCP<DRT::Discretization>&        xfluiddis,          ///< background fluid discretization
  const Teuchos::RCP<LINALG::Solver>&             solver,             ///< fluid solver
  const Teuchos::RCP<Teuchos::ParameterList>&     params,             ///< xfluid params
  const Teuchos::RCP<IO::DiscretizationWriter>&   output,             ///< discretization writer for paraview output
  bool                                            ale_xfluid,         ///< background (XFEM) fluid in ALE-formulation
  bool                                            ale_fluid           ///< embedded fluid in ALE-formulation
) : XFluid(
  xfluiddis,
  embedded_fluid->Discretization(),
  solver,
  params,
  output,
  ale_xfluid),
  embedded_fluid_(embedded_fluid),
  alefluid_(ale_fluid)
{
  // Write background fluid mesh
 // output_->WriteMesh(0,0.0);
}

FLD::XFluidFluidNew::~XFluidFluidNew()
{
}

void FLD::XFluidFluidNew::Init()
{
  // set parameters specific for fluid-fluid coupling
  SetXFluidFluidParams();

  //DRT::UTILS::PrintParallelDistribution(*boundarydis_);

  embedded_fluid_->Init();

  // base class init
  XFluid::Init();

  if(meshcoupl_dis_.size() != 1) dserror("we expect exact one mesh coupling discretization for Xfluidfluid at the moment!");

  soliddis_ = meshcoupl_dis_[0];

  // make the dofset of boundarydis be a subset of the embedded dis
  Teuchos::RCP<Epetra_Map> newcolnodemap = DRT::UTILS::ComputeNodeColMap(
      soliddis_,boundarydis_);
  soliddis_->Redistribute(*(embedded_fluid_->Discretization()->NodeRowMap()), *newcolnodemap);
  Teuchos::RCP<DRT::DofSet> newdofset=Teuchos::rcp(new DRT::TransparentIndependentDofSet(
      soliddis_,false));
  boundarydis_->ReplaceDofSet(newdofset); // do not call this with true!!
  boundarydis_->FillComplete();

  INPAR::FLUID::CalcError calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_,"calculate error");
  if ((coupling_strategy_ != INPAR::XFEM::Xfluid_Sided_Coupling) or (calcerr != INPAR::FLUID::no_error_calculation))
  {
    PrepareEmbeddedDistribution();
    CreateBoundaryEmbeddedMap();
  }

  //-------------------------------------------------------------------

  if (alefluid_)
    dispnpoldstate_ = Teuchos::rcp(new Epetra_Vector(*embedded_fluid_->Dispnp()));

  //--------------------------------------------------
  // Create XFluidFluid State
  //-----------------------------------------------
  const int restart = DRT::Problem::Instance()->Restart();

  if (restart)
  {
    embedded_fluid_->ReadRestart(restart);
  }

  state_ = this->GetNewState();

  //----------------------------------------------------------------
  // create auxiliary discretization of outer embedded fluid elements
  // Nitsche-Eigenvalue-Problem and for application of EOS
  // pressure term to outer embedded element
  //----------------------------------------------------------------
  // Todo: Shift to XFEM::MeshCouplingFluidFluid
  if (nitsche_evp_ || xff_eos_pres_emb_layer_)
    CreateEmbeddedBoundaryDiscretization();

  // non-stationary fluid-fluid interface requires
  // proper time integration approach
  if (alefluid_)
    xfluidfluid_timeint_ =  Teuchos::rcp(new XFEM::XFluidFluidTimeIntegration(
        XFluid::discret_,
        embedded_fluid_->Discretization(),
        state_->Wizard(), step_,
        xfem_timeintapproach_,*params_));

}

void FLD::XFluidFluidNew::SetDirichletNeumannBC()
{
  XFluid::SetDirichletNeumannBC();
//
//  Teuchos::ParameterList eleparams;
//  // other parameters needed by the elements
//  eleparams.set("total time",time_);
//
//  soliddis_->ClearState();
//  soliddis_->SetState("velaf",embedded_fluid_->Velnp());
//  //don't call this with the mapextractor. Otherwise the Mapextractor will
//  //be built again.
//  soliddis_->EvaluateDirichlet(eleparams,embedded_fluid_->Velnp(),Teuchos::null,Teuchos::null,Teuchos::null);
//  soliddis_->ClearState();
}

void FLD::XFluidFluidNew::SetXFluidFluidParams()
{
  Teuchos::ParameterList&   params_xf_stab = params_->sublist("XFLUID DYNAMIC/STABILIZATION");

  hybrid_lm_l2_proj_     = DRT::INPUT::IntegralValue<INPAR::XFEM::Hybrid_LM_L2_Proj>(params_xf_stab, "HYBRID_LM_L2_PROJ");

  xff_conv_stab_scaling_ = DRT::INPUT::IntegralValue<INPAR::XFEM::XFF_ConvStabScaling>(params_xf_stab,"XFF_CONV_STAB_SCALING");

  // Flag, whether any face-based terms are active
  eval_eos_ = edge_based_ || ghost_penalty_;

  // additional eos pressure stabilization on the elements of the embedded discretization,
  // that contribute to the interface
  xff_eos_pres_emb_layer_ = DRT::INPUT::IntegralValue<bool>(params_xf_stab,"XFF_EOS_PRES_EMB_LAYER");

  // whether an eigenvalue problem has to be solved to estimate Nitsche's parameter
  nitsche_evp_ = (DRT::INPUT::IntegralValue<INPAR::XFEM::ViscStab_TraceEstimate>(params_xf_stab,"VISC_STAB_TRACE_ESTIMATE")
                  == INPAR::XFEM::ViscStab_TraceEstimate_eigenvalue);

  // get general XFEM/XFFSI specific parameters
  monolithic_approach_  = DRT::INPUT::IntegralValue<INPAR::XFEM::Monolithic_xffsi_Approach>(params_->sublist("XFLUID DYNAMIC/GENERAL"),"MONOLITHIC_XFFSI_APPROACH");
  xfem_timeintapproach_ = DRT::INPUT::IntegralValue<INPAR::XFEM::XFluidFluidTimeInt>(params_->sublist("XFLUID DYNAMIC/GENERAL"),"XFLUIDFLUID_TIMEINT");

  // get information about active shape derivatives
  active_shapederivatives_ = alefluid_ && params_->get<bool>("shape derivatives");
}

void FLD::XFluidFluidNew::TimeUpdate()
{
  XFluid::Update();
  embedded_fluid_->Update();
}

void FLD::XFluidFluidNew::CutAndSetStateVectors()
{
  if (!alefluid_)
   return;

  // save the old state
  staten_ = state_;

  const int restart = DRT::Problem::Instance()->Restart();
  // if restart
  if (restart and ((restart+1) == step_))
  {
    xfluidfluid_timeint_->CreateBgNodeMapsForRestart(staten_->Wizard());
  }
}

Teuchos::RCP<LINALG::BlockSparseMatrixBase> FLD::XFluidFluidNew::BlockSystemMatrix(
    Teuchos::RCP<Epetra_Map> innermap,
    Teuchos::RCP<Epetra_Map> condmap)
{
  //Map of fluid FSI DOFs: condmap
  //Map of inner fluid DOFs: innermap

  //Get the fluid-fluid system matrix as sparse matrix
  Teuchos::RCP<LINALG::SparseMatrix> sparsesysmat = SystemMatrix();

  //F_{II}, F_{I\Gamma}, F_{\GammaI}, F_{\Gamma\Gamma}
  Teuchos::RCP<LINALG::SparseMatrix> fii, fig, fgi, fgg;
  // Split sparse system matrix into blocks according to the given maps
  LINALG::SplitMatrix2x2(sparsesysmat,innermap,condmap,innermap,condmap,fii,fig,fgi,fgg);
  // create a new block matrix out of the 4 blocks
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockmat = LINALG::BlockMatrix2x2(*fii,*fig,*fgi,*fgg);

  if ( blockmat == Teuchos::null )
    dserror("Creation of fluid-fluid block matrix failed.");

  return blockmat;
}

Teuchos::RCP<const Epetra_Map> FLD::XFluidFluidNew::PressureRowMap()
{
  return state_->xffluidvelpressplitter_->CondMap();
}

Teuchos::RCP<const Epetra_Map> FLD::XFluidFluidNew::VelocityRowMap()
{
  return state_->xffluidvelpressplitter_->OtherMap();
}

Teuchos::RCP<DRT::ResultTest> FLD::XFluidFluidNew::CreateFieldTest()
{
  return Teuchos::rcp(new FLD::XFluidFluidResultTest(*this));
}

Teuchos::RCP<FLD::XFluidFluidState> FLD::XFluidFluidNew::GetNewState()
{
  if (alefluid_)
  {
    LINALG::Export(*dispnp_,*idispnp_);
  }

  Teuchos::RCP<FLD::XFluidFluidState> state = state_creator_->Create(
    xdiscret_,
    embedded_fluid_->Discretization(),
    Teuchos::null, //!< col vector holding background ALE displacements for backdis
    solver_->Params(),
    step_,
    time_);

  // increment vector for merged background & embedded fluid
  // (not the classical Newton increment but the difference to
  // the value at the last time step)
  stepinc_ = LINALG::CreateVector(*state->xffluiddofrowmap_,true);

  // build a merged map from fluid-fluid dbc-maps
  state->CreateMergedDBCMapExtractor(embedded_fluid_->GetDBCMapExtractor());

  // create object for edgebased stabilization
  if (eval_eos_)
    edgestab_ =  Teuchos::rcp(new XFEM::XFEM_EdgeStab(state->Wizard(), discret_, include_inner_));

  return state;
}

void FLD::XFluidFluidNew::CreateState()
{
  // new cut for this time step
  state_ = this->GetNewState();
  XFluid::state_ = state_;

  if (!alefluid_)
    return;

  // map of background-fluid's standard and enriched node-ids and
  // their dof-gids for new cut
  xfluidfluid_timeint_->SaveBgNodeMapsAndCreateNew(state_->Wizard());
}

void FLD::XFluidFluidNew::AssembleMatAndRHS(
    int itnum                           ///< iteration number
)
{
  // evaluate elements
  soliddis_->SetState("velaf",embedded_fluid_->Velnp());
  embedded_fluid_->PrepareSolve();
  XFluid::AssembleMatAndRHS(itnum);

  // merge the residuals and velnp each into one large Epetra_Vector for the composed system
  state_->xffluidsplitter_->InsertXFluidVector(state_->velnp_, state_->xffluidvelnp_);
  state_->xffluidsplitter_->InsertFluidVector(embedded_fluid_->Velnp(), state_->xffluidvelnp_);

  state_->xffluidsplitter_->InsertXFluidVector(state_->residual_, state_->xffluidresidual_);
  state_->xffluidsplitter_->InsertFluidVector(embedded_fluid_->Residual(), state_->xffluidresidual_);

  // build a merged dbc map extractor from both discretizations
  state_->CreateMergedDBCMapExtractor(embedded_fluid_->GetDBCMapExtractor());
}

void FLD::XFluidFluidNew::PrepareShapeDerivatives(
  const Teuchos::RCP<const LINALG::MultiMapExtractor> fsiextractor,
  const Teuchos::RCP<std::set<int> > condelements)
{
  if (! active_shapederivatives_)
    return;

  // here we initialize the shapederivates
  // REMARK: the shape derivatives matrix results from linearization w.r.t. ALE-displacements
  // and therefore solely knows ALE-dof - here we use "extended shapederivatives" including
  // background fluid entries, that are set to zero
  Teuchos::RCP<LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy> > mat =
    Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>(*fsiextractor,*fsiextractor,108,false,true));
  mat->SetCondElements(condelements);
  extended_shapederivatives_ = mat;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool FLD::XFluidFluidNew::ConvergenceCheck(
  int          itnum,
  int          itemax,
  const double ittol)
{

  bool stopnonliniter = XFluid::ConvergenceCheck(itnum,itemax,ittol);

  //--------- Apply Dirichlet boundary conditions to system of equations
  //          residual displacements are supposed to be zero at
  //          boundary conditions

  const int mc_idx=0;
  Teuchos::RCP<XFluidState::CouplingState> & coup_state = state_->coup_state_[mc_idx];

  // adding rhC_s_ (coupling contribution) to residual of embedded fluid
  for (int iter=0; iter<coup_state->rhC_s_->MyLength();++iter)
  {
    Teuchos::RCP<Epetra_Vector> emb_residual = embedded_fluid_->Residual();
    const int rhsdgid = coup_state->rhC_s_->Map().GID(iter);
    if (coup_state->rhC_s_->Map().MyGID(rhsdgid) == false) dserror("rhsd_ should be on all processors");
    if (emb_residual->Map().MyGID(rhsdgid))
      (*emb_residual)[emb_residual->Map().LID(rhsdgid)] =
          (*emb_residual)[emb_residual->Map().LID(rhsdgid)] +
          (*coup_state->rhC_s_)[coup_state->rhC_s_->Map().LID(rhsdgid)];
    else dserror("Interface dof %d does not belong to embedded discretization!",rhsdgid);
  }

  // assemble subsequent system matrices into one
  state_->xffluidsysmat_->Zero();
  state_->xffluidsysmat_->Add(*state_->sysmat_,false,1.0,0.0);
  state_->xffluidsysmat_->Add(*embedded_fluid_->SystemMatrix(),false,1.0,1.0);
  state_->xffluidsysmat_->Add(*coup_state->C_xs_,false,1.0,1.0);
  state_->xffluidsysmat_->Add(*coup_state->C_sx_,false,1.0,1.0);
  state_->xffluidsysmat_->Add(*coup_state->C_ss_,false,1.0,1.0);
  state_->xffluidsysmat_->Complete();

  return stopnonliniter;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidNew::UpdateByIncrement()
{
  state_->xffluidvelnp_->Update(1.0,*state_->xffluidincvel_,1.0);

  // extract velnp_
  state_->velnp_ = state_->xffluidsplitter_->ExtractXFluidVector(state_->xffluidvelnp_);
  embedded_fluid_->WriteAccessVelnp() = state_->xffluidsplitter_->ExtractFluidVector(state_->xffluidvelnp_);

  // extract residual
  state_->residual_ = state_->xffluidsplitter_->ExtractXFluidVector(state_->xffluidresidual_);
  embedded_fluid_->Residual() = state_->xffluidsplitter_->ExtractFluidVector(state_->xffluidresidual_);

  // Update the fluid material velocity along the interface (ivelnp_), source (in): state_.velnp_
  LINALG::Export(*(embedded_fluid_->Velnp()),*(ivelnp_));
  boundarydis_->SetState("ivelnp",ivelnp_);
}
