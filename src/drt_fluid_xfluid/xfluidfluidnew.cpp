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

#include "xfluidresulttest.H"

#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid/fluid_utils.H"

#include "../drt_lib/drt_discret_faces.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_dofset_transparent_independent.H"
#include "../drt_lib/drt_utils_parallel.H"

#include "../drt_io/io_pstream.H"
#include "../drt_io/io.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

#include "../drt_fluid_ele/fluid_ele_parameter_xfem.H"

#include "../drt_xfem/xfem_edgestab.H"
#include "../drt_xfem/xfem_dofset.H"
#include "../drt_xfem/xfem_condition_manager.H"
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
  xfluiddis->Writer(),
  ale_xfluid),
  embedded_fluid_(embedded_fluid),
  ale_embfluid_(ale_fluid)
{
  xfluiddis->Writer()->WriteMesh(0,0.0);
}

FLD::XFluidFluidNew::~XFluidFluidNew()
{
}

void FLD::XFluidFluidNew::Init()
{
  // set parameters specific for fluid-fluid coupling
  SetXFluidFluidParams();

  // initialize embedded fluid
  embedded_fluid_->Init();

  // base class init
  XFluid::Init();

  if (meshcoupl_dis_.size() != 1) dserror("we expect exactly one mesh coupling discretization for Xfluidfluid at the moment!");

  if (DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_,"calculate error") != INPAR::FLUID::no_error_calculation)
  {
    mc_xff_->RedistributeForErrorCalculation();
  }

  // recreate internal faces of DiscretizationFaces (as the distribution of the embedded
  // discretization may have changed)
  bool redistributed =
      (DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_,"calculate error") != INPAR::FLUID::no_error_calculation
      || mc_xff_->GetAveragingStrategy() == INPAR::XFEM::Embedded_Sided
      || mc_xff_->GetAveragingStrategy() == INPAR::XFEM::Mean);
  if (redistributed)
    embedded_fluid_->CreateFacesExtension();

  //--------------------------------------------------
  // Create XFluidFluid State
  //-----------------------------------------------
  const int restart = DRT::Problem::Instance()->Restart();

  if (restart)
  {
    embedded_fluid_->ReadRestart(restart);
  }

  if (ale_embfluid_)
    dispnpoldstate_ = Teuchos::rcp(new Epetra_Vector(*embedded_fluid_->Dispnp()));
}

void FLD::XFluidFluidNew::SetXFluidFluidParams()
{
  Teuchos::ParameterList&   params_xf_stab = params_->sublist("XFLUID DYNAMIC/STABILIZATION");

  hybrid_lm_l2_proj_     = DRT::INPUT::IntegralValue<INPAR::XFEM::Hybrid_LM_L2_Proj>(params_xf_stab, "HYBRID_LM_L2_PROJ");

  xff_conv_stab_scaling_ = DRT::INPUT::IntegralValue<INPAR::XFEM::XFF_ConvStabScaling>(params_xf_stab,"XFF_CONV_STAB_SCALING");

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
  active_shapederivatives_ = ale_embfluid_ && params_->get<bool>("shape derivatives");
}

void FLD::XFluidFluidNew::SetDirichletNeumannBC()
{
  XFluid::SetDirichletNeumannBC();
  embedded_fluid_->SetDirichletNeumannBC();
}

void FLD::XFluidFluidNew::SetInterfaceFixed()
{
  mc_xff_->SetInterfaceFixed();
}

void FLD::XFluidFluidNew::SetInterfaceFree()
{
  mc_xff_->SetInterfaceFree();
}

void FLD::XFluidFluidNew::PrepareTimeStep()
{
  XFluid::PrepareTimeStep();
  embedded_fluid_->PrepareTimeStep();
}

void FLD::XFluidFluidNew::PrepareSolve()
{
  XFluid::PrepareSolve();

  // merge the velnp each into one large Epetra_Vector for the composed system
  xff_state_->xffluidsplitter_->InsertXFluidVector(xff_state_->velnp_, xff_state_->xffluidvelnp_);
  xff_state_->xffluidsplitter_->InsertFluidVector(embedded_fluid_->Velnp(), xff_state_->xffluidvelnp_);

  if (mc_xff_->GetAveragingStrategy() == INPAR::XFEM::Embedded_Sided
      and nitsche_evp_)
  {
    mc_xff_->EstimateNitscheTraceMaxEigenvalue(embedded_fluid_->Dispnp());
  }
}

void FLD::XFluidFluidNew::Evaluate(
  Teuchos::RCP<const Epetra_Vector> stepinc  ///< solution increment between time step n and n+1
)
{
  // split step increment
  Teuchos::RCP<Epetra_Vector> stepinc_xfluid  = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> stepinc_emb = Teuchos::null;
  if (stepinc!=Teuchos::null)
  {
    stepinc_xfluid = xff_state_->xffluidsplitter_->ExtractXFluidVector(stepinc);
    stepinc_emb = xff_state_->xffluidsplitter_->ExtractFluidVector(stepinc);

    // iteration increment
    xff_state_->xffluidincvel_->Update(1.0,*stepinc,-1.0,*stepinc_,0.0);
  }

  // call evaluate on both fluid fields
  embedded_fluid_->Evaluate(stepinc_emb);

  XFluid::Evaluate(stepinc_xfluid);

  if (active_shapederivatives_)
  {
    extended_shapederivatives_->Add(*embedded_fluid_->ShapeDerivatives(),false,1.0,0.0);
  }

  // step increment
  stepinc_->Update(1.0,*xff_state_->xffluidvelnp_,-1.0,*xff_state_->xffluidveln_,0.0);

}

void FLD::XFluidFluidNew::TimeUpdate()
{
  XFluid::TimeUpdate();
  embedded_fluid_->TimeUpdate();

  xff_state_->xffluidveln_->Update(1.0,*xff_state_->xffluidvelnp_,0.0);
  xff_state_->xffluidincvel_->PutScalar(0.0);
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
  return xff_state_->xffluidvelpressplitter_->CondMap();
}

Teuchos::RCP<const Epetra_Map> FLD::XFluidFluidNew::VelocityRowMap()
{
  return xff_state_->xffluidvelpressplitter_->OtherMap();
}

Teuchos::RCP<DRT::ResultTest> FLD::XFluidFluidNew::CreateFieldTest()
{
  return Teuchos::rcp(new FLD::XFluidResultTest(*this));
}

Teuchos::RCP<FLD::XFluidState> FLD::XFluidFluidNew::GetNewState()
{
  // further use type-cast pointer to MeshCouplingFluidFluid
  if (mc_xff_ == Teuchos::null)
  {
    const std::string cond_name("XFEMSurfFluidFluid");
    mc_xff_ = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFluidFluid>(condition_manager_->GetMeshCoupling(cond_name));

    if (mc_xff_ == Teuchos::null)
      dserror("Failed to cast to MeshCouplingFluidFluid");
  }
  if (ale_embfluid_)
  {
    LINALG::Export(*embedded_fluid_->Dispnp(),*mc_xff_->IDispnp());
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

  return state;
}

void FLD::XFluidFluidNew::CreateState()
{
  // new cut for this time step

  XFluid::CreateState();
  xff_state_ = Teuchos::rcp_dynamic_cast<FLD::XFluidFluidState>(XFluid::state_);

  if (xff_state_ == Teuchos::null)
    dserror("Failed to create an instance of XFluidFluidState.");

  if (!ale_embfluid_)
    return;

  if (xfluidfluid_timeint_ == Teuchos::null)
    xfluidfluid_timeint_ =  Teuchos::rcp(new XFEM::XFluidFluidTimeIntegration(
      XFluid::discret_, embedded_fluid_->Discretization(),
      xff_state_->Wizard(), step_,
      xfem_timeintapproach_,*params_));

  // map of background-fluid's standard and enriched node-ids and
  // their dof-gids for new cut
  samemaps_ = xfluidfluid_timeint_->SaveBgNodeMapsAndCreateNew(xff_state_->Wizard());
}

void FLD::XFluidFluidNew::AssembleMatAndRHS(
  int itnum                           ///< iteration number
)
{
  TEUCHOS_FUNC_TIME_MONITOR("FLD::XFluidFluidNew::AssembleMatAndRHS");

  if (mc_xff_->GetAveragingStrategy() ==
      INPAR::XFEM::Embedded_Sided)
  {
    mc_xff_->GetCouplingDis()->ClearState();
    // set velocity and displacement state for embedded fluid
    embedded_fluid_->SetStateTimInt();
    mc_xff_->GetCouplingDis()->SetState("veln",embedded_fluid_->Veln());
    if (ale_embfluid_)
      mc_xff_->GetCouplingDis()->SetState(
        "dispnp",embedded_fluid_->Dispnp());
  }

  // export interface velocities
  LINALG::Export(*(embedded_fluid_->Velnp()),*(mc_xff_->IVelnp()));
  LINALG::Export(*(embedded_fluid_->Veln()),*(mc_xff_->IVeln()));

  // evaluate elements of XFluid part
  XFluid::AssembleMatAndRHS(itnum);

  // insert XFluid residual to merged
  xff_state_->xffluidsplitter_->InsertXFluidVector(xff_state_->residual_,
    xff_state_->xffluidresidual_);

  // evaluate elements of embedded fluid
  embedded_fluid_->PrepareSolve();

  // add coupling contribution embedded residual
  const int mc_idx = 0;
  Teuchos::RCP<XFluidState::CouplingState> & coup_state = xff_state_->coup_state_[mc_idx];

  {
    Teuchos::RCP<Epetra_Vector> emb_residual = LINALG::CreateVector(*embedded_fluid_->DofRowMap(),true);
    emb_residual->Update(1.0,*embedded_fluid_->Residual(),0.0);

    // adding rhC_s_ (coupling contribution) to residual of embedded fluid
    for (int iter=0; iter<coup_state->rhC_s_->MyLength();++iter)
    {
      const int rhsgid = coup_state->rhC_s_->Map().GID(iter);
      if (coup_state->rhC_s_->Map().MyGID(rhsgid) == false) dserror("rhC_s_ should be on all processors");
      if (emb_residual->Map().MyGID(rhsgid))
        (*emb_residual)[emb_residual->Map().LID(rhsgid)] +=
            (*coup_state->rhC_s_)[coup_state->rhC_s_->Map().LID(rhsgid)];
      else dserror("Interface dof %d does not belong to embedded discretization!",rhsgid);
    }

    // add embedded part of merged residual
    xff_state_->xffluidsplitter_->InsertFluidVector(emb_residual,
        xff_state_->xffluidresidual_);
  }

  // assemble XFluid and embedded fluid system matrices into one
  xff_state_->xffluidsysmat_->Zero();
  xff_state_->xffluidsysmat_->Add(*xff_state_->sysmat_,false,1.0,0.0);
  xff_state_->xffluidsysmat_->Add(*embedded_fluid_->SystemMatrix(),false,1.0,1.0);
  xff_state_->xffluidsysmat_->Add(*coup_state->C_xs_,false,1.0,1.0);
  xff_state_->xffluidsysmat_->Add(*coup_state->C_sx_,false,1.0,1.0);
  xff_state_->xffluidsysmat_->Add(*coup_state->C_ss_,false,1.0,1.0);
  xff_state_->xffluidsysmat_->Complete();
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
void FLD::XFluidFluidNew::UpdateByIncrement()
{
  // Update merged
  XFluid::UpdateByIncrement();
  // update xfluid
  xff_state_->velnp_->Update(1.0,*xff_state_->xffluidsplitter_->ExtractXFluidVector(xff_state_->xffluidvelnp_),0.0);
//  // update embedded fluid
//  //embedded_fluid_->IterUpdate(xff_state_->xffluidsplitter_->ExtractFluidVector(xff_state_->xffluidincvel_));
  embedded_fluid_->WriteAccessVelnp()->Update(1.0,*xff_state_->xffluidsplitter_->ExtractFluidVector(xff_state_->xffluidvelnp_),0.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidNew::XTimint_StoreOldStateData(const bool firstcall_in_timestep)
{
  // save the old state
  xff_staten_ = xff_state_;

  const int restart = DRT::Problem::Instance()->Restart();

  // if restart
  if (restart and ((restart+1) == step_))
  {
    xfluidfluid_timeint_->CreateBgNodeMapsForRestart(xff_staten_->Wizard());
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidNew::XTimint_DoTimeStepTransfer(const bool screen_out)
{
  //---------------------------------------------------------------
  if(myrank_==0 and screen_out) IO::cout << "XFEM::TIMEINTEGRATION: ..." << IO::endl;

  //---------------------------------------------------------------
  if(timealgo_ !=  INPAR::FLUID::timeint_one_step_theta) dserror("check which vectors have to be reconstructed for non-OST scheme");

  switch (monolithic_approach_)
  {
  case INPAR::XFEM::XFFSI_Full_Newton:
    SetXFluidStateVectors(embedded_fluid_->Dispnp());
    break;
  case INPAR::XFEM::XFFSI_FixedALE_Interpolation:
    SetEmbFluidStateVectors(dispnpoldstate_);
    SetXFluidStateVectors(dispnpoldstate_);
    break;
  case INPAR::XFEM::XFFSI_FixedALE_Partitioned:
    SetXFluidStateVectors(dispnpoldstate_);
    break;
  default:
    dserror("Unknown monolithic XFFSI approach.");
    break;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidNew::Output()
{
  XFluid::Output();
  embedded_fluid_->Output();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidNew::SetXFluidStateVectors(Teuchos::RCP<const Epetra_Vector> disp)
{
  if (xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_FullProj or
      xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_KeepGhostValues or
     (xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_ProjIfMoved and (not samemaps_)) or
      xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_IncompProj)
  {
    // export the vectors to the column distribution map
    Teuchos::RCP<Epetra_Vector> velnpcol = LINALG::CreateVector(*embedded_fluid_->Discretization()->DofColMap(),true);
    Teuchos::RCP<Epetra_Vector> velncol  = LINALG::CreateVector(*embedded_fluid_->Discretization()->DofColMap(),true);
    Teuchos::RCP<Epetra_Vector> accnpcol = LINALG::CreateVector(*embedded_fluid_->Discretization()->DofColMap(),true);
    Teuchos::RCP<Epetra_Vector> velnmcol  = LINALG::CreateVector(*embedded_fluid_->Discretization()->DofColMap(),true);
    Teuchos::RCP<Epetra_Vector> accncol = LINALG::CreateVector(*embedded_fluid_->Discretization()->DofColMap(),true);
    Teuchos::RCP<Epetra_Vector> dispcol  = LINALG::CreateVector(*embedded_fluid_->Discretization()->DofColMap(),true);

    LINALG::Export(*embedded_fluid_->Velnp(),*velnpcol);
    LINALG::Export(*embedded_fluid_->Veln(), *velncol);
    LINALG::Export(*embedded_fluid_->Velnm(), *velnmcol);
    LINALG::Export(*embedded_fluid_->Accnp(),*accnpcol);
    LINALG::Export(*embedded_fluid_->Accn(),*accncol);
    LINALG::Export(*disp,  *dispcol);

    // Todo: we export far more than we need...

    // we have five state vectors, which need values from the last time step
    xfluidfluid_timeint_->SetNewBgStateVectors(xff_state_->velnp_,
                                               xff_state_->veln_,
                                               xff_state_->velnm_,
                                               xff_state_->accnp_,
                                               xff_state_->accn_,
                                               xff_staten_->velnp_,
                                               xff_staten_->veln_,
                                               xff_staten_->velnm_,
                                               xff_staten_->accnp_,
                                               xff_staten_->accn_,
                                               velnpcol,
                                               velncol,
                                               velnmcol,
                                               accnpcol,
                                               accncol,
                                               dispcol);

    // project and enforce incompressibility on projection patch
    if (xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_IncompProj)
    {
      xfluidfluid_timeint_->EnforceIncompAfterProjection(
          xff_state_->Wizard(),
          xff_staten_->Wizard(),
          xff_state_->velnp_,
          xff_state_->veln_,
          xff_state_->velnm_,
          xff_state_->dbcmaps_);

    }
  }
  // Note: if Xff_TimeInt_ProjIfMoved is chosen and the maps remain the same
  // (TODO: they remain the same just for one dofset)
  // the enriched values are not projected from the embedded fluid anymore.
  else if (xfem_timeintapproach_ == INPAR::XFEM::Xff_TimeInt_ProjIfMoved and samemaps_)
  {
    // we use the old velocity as start value
    xff_state_->velnp_->Update(1.0,*xff_staten_->velnp_,0.0);
    xff_state_->veln_->Update( 1.0,*xff_staten_->veln_, 0.0);
    xff_state_->velnm_->Update(1.0,*xff_staten_->velnm_,0.0);
    xff_state_->accn_->Update( 1.0,*xff_staten_->accn_, 0.0);
    xff_state_->accnp_->Update(1.0,*xff_staten_->accnp_,0.0);
  }

  // write interpolated state vectors into merged veln and velnp
  xff_state_->xffluidsplitter_->InsertFluidVector(embedded_fluid_->Veln(),xff_state_->xffluidveln_);
  xff_state_->xffluidsplitter_->InsertXFluidVector(xff_state_->veln_,xff_state_->xffluidveln_);
  xff_state_->xffluidsplitter_->InsertXFluidVector(xff_state_->velnp_, xff_state_->xffluidvelnp_);
  xff_state_->xffluidsplitter_->InsertFluidVector(embedded_fluid_->Velnp(), xff_state_->xffluidvelnp_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluidFluidNew::SetEmbFluidStateVectors(Teuchos::RCP<Epetra_Vector> disp)
{
  // get nodal velocities and pressure from previous intersection state

  // export the vectors to the column distribution map
  Teuchos::RCP<Epetra_Vector> embvelnpcol  = LINALG::CreateVector(*embedded_fluid_->Discretization()->DofColMap());
  Teuchos::RCP<Epetra_Vector> embvelncol   = LINALG::CreateVector(*embedded_fluid_->Discretization()->DofColMap());
  Teuchos::RCP<Epetra_Vector> embaccnpcol  = LINALG::CreateVector(*embedded_fluid_->Discretization()->DofColMap());
  Teuchos::RCP<Epetra_Vector> embvelnmcol   = LINALG::CreateVector(*embedded_fluid_->Discretization()->DofColMap());
  Teuchos::RCP<Epetra_Vector> embaccncol   = LINALG::CreateVector(*embedded_fluid_->Discretization()->DofColMap());
  Teuchos::RCP<Epetra_Vector> embdispcol   = LINALG::CreateVector(*embedded_fluid_->Discretization()->DofColMap());
  Teuchos::RCP<Epetra_Vector> embdispnpcol = LINALG::CreateVector(*embedded_fluid_->Discretization()->DofColMap());

  LINALG::Export(*embedded_fluid_->Velnp(), *embvelnpcol);
  LINALG::Export(*embedded_fluid_->Veln(),  *embvelncol);
  if (embedded_fluid_->Velnm() != Teuchos::null)
    LINALG::Export(*embedded_fluid_->Velnm(), *embaccnpcol);
  LINALG::Export(*embedded_fluid_->Accnp(), *embaccnpcol);
  LINALG::Export(*embedded_fluid_->Accn(), *embaccnpcol);
  LINALG::Export(*disp,                     *embdispcol);
  if (embedded_fluid_->Dispnp() != Teuchos::null)
    LINALG::Export(*embedded_fluid_->Dispnp(),*embdispnpcol);

  xfluidfluid_timeint_->SetNewEmbStateVectors(
    embedded_fluid_->WriteAccessVelnp(),
    Teuchos::null,
    Teuchos::null,
    Teuchos::null,
    Teuchos::null,
    embvelnpcol,
    embvelncol,
    embvelnmcol,
    embaccnpcol,
    embaccncol,
    xff_staten_->velnp_,
    xff_staten_->veln_,
    xff_staten_->velnm_,
    xff_staten_->accnp_,
    xff_staten_->accn_,
    embdispnpcol,
    embdispcol);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::XFluidFluidNew::UpdateMonolithicFluidSolution(const Teuchos::RCP<const Epetra_Map>& fsidofmap)
{
  Teuchos::RCP<Epetra_Vector> fixedfsizeros = LINALG::CreateVector(*fsidofmap,true);
  // manipulate the dbc map extractor
  Teuchos::RCP<const Epetra_Map> dbcmap = embedded_fluid_->GetDBCMapExtractor()->CondMap();
  std::vector<Teuchos::RCP<const Epetra_Map> > condmaps;
  condmaps.push_back(dbcmap);
  condmaps.push_back(fsidofmap);
  Teuchos::RCP<Epetra_Map> condmerged = LINALG::MultiMapExtractor::MergeMaps(condmaps);

  Teuchos::RCP<LINALG::MapExtractor> fsidbcmapex =
      Teuchos::rcp(new LINALG::MapExtractor(*(embedded_fluid_->DofRowMap()), condmerged));

  xff_state_->CreateMergedDBCMapExtractor(fsidbcmapex);
  XFluid::SetElementTimeParameter();
  Solve();
  xff_state_->CreateMergedDBCMapExtractor(embedded_fluid_->GetDBCMapExtractor());
}//FLD::XFluidFluid::UpdateMonolithicFluidSolution()
