/*----------------------------------------------------------------------*/
/*! \file

\brief Control routine for fluid-fluid (in)stationary solvers with XFEM

\level 2

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
 */
/*----------------------------------------------------------------------*/

#include "xfluidfluid.H"
#include "xfluid_state_creator.H"
#include "xfluidresulttest.H"

#include "../drt_lib/drt_discret_faces.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_dofset_transparent_independent.H"
#include "../drt_lib/drt_utils_parallel.H"

#include "../drt_io/io_pstream.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

#include "../drt_fluid/fluid_utils.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_ele/fluid_ele_interface.H"
#include "../drt_fluid_ele/fluid_ele_factory.H"

#include "../drt_xfem/xfem_edgestab.H"
#include "../drt_xfem/xfem_dofset.H"
#include "../drt_xfem/xfem_condition_manager.H"
#include "../drt_xfem/xfem_mesh_projector.H"
#include "../drt_xfem/xfluid_timeInt.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FLD::XFluidFluid::XFluidFluid(const Teuchos::RCP<FLD::FluidImplicitTimeInt>& embedded_fluid,
    const Teuchos::RCP<DRT::Discretization>& xfluiddis, const Teuchos::RCP<LINALG::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params, bool ale_xfluid, bool ale_fluid)
    : XFluid(xfluiddis, embedded_fluid->Discretization(), Teuchos::null, solver, params,
          xfluiddis->Writer(), ale_xfluid),
      embedded_fluid_(embedded_fluid),
      projector_(Teuchos::rcp(
          new XFEM::MeshProjector(embedded_fluid_->Discretization(), discret_, *params_))),
      ale_embfluid_(ale_fluid),
      cond_name_("XFEMSurfFluidFluid")
{
  xfluiddis->Writer()->WriteMesh(0, 0.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FLD::XFluidFluid::XFluidFluid(const Teuchos::RCP<FLD::FluidImplicitTimeInt>& embedded_fluid,
    const Teuchos::RCP<DRT::Discretization>& xfluiddis,
    const Teuchos::RCP<DRT::Discretization>& soliddis, const Teuchos::RCP<LINALG::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params, bool ale_xfluid, bool ale_fluid)
    : XFluid(xfluiddis, soliddis, Teuchos::null, solver, params, xfluiddis->Writer(), ale_xfluid),
      embedded_fluid_(embedded_fluid),
      ale_embfluid_(ale_fluid),
      cond_name_("XFEMSurfFluidFluid")
{
  xfluiddis->Writer()->WriteMesh(0, 0.0);
  meshcoupl_dis_.push_back(embedded_fluid->Discretization());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FLD::XFluidFluid::~XFluidFluid() {}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluidFluid::Init(bool createinitialstate)
{
  // initialize embedded fluid
  embedded_fluid_->Init();

  // base class init
  XFluid::Init(false);

  // set parameters specific for fluid-fluid coupling
  SetXFluidFluidParams();


  if (createinitialstate) CreateInitialState();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluidFluid::CreateInitialState()
{
  // base class CreateInitialState
  XFluid::CreateInitialState();

  if (DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_, "calculate error") !=
      INPAR::FLUID::no_error_calculation)
  {
    mc_xff_->RedistributeForErrorCalculation();
  }

  // recreate internal faces of DiscretizationFaces (as the distribution of the embedded
  // discretization may have changed)
  if (DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_, "calculate error") !=
          INPAR::FLUID::no_error_calculation ||
      mc_xff_->GetAveragingStrategy() == INPAR::XFEM::Embedded_Sided ||
      mc_xff_->GetAveragingStrategy() == INPAR::XFEM::Mean)
  {
    embedded_fluid_->CreateFacesExtension();
  }

  // create internal faces for embedded discretization afterwards, if not full EOS on embedded
  // domain
  {
    Teuchos::ParameterList* stabparams = &(params_->sublist("RESIDUAL-BASED STABILIZATION"));
    if (xff_eos_pres_emb_layer_ && DRT::INPUT::IntegralValue<INPAR::FLUID::StabType>(*stabparams,
                                       "STABTYPE") == INPAR::FLUID::stabtype_residualbased)
    {
      Teuchos::RCP<DRT::DiscretizationFaces> facediscret =
          Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(
              embedded_fluid_->Discretization(), true);
      facediscret->CreateInternalFacesExtension(true);
    }
  }

  //--------------------------------------------------
  // Create XFluidFluid State
  //-----------------------------------------------
  const int restart = DRT::Problem::Instance()->Restart();

  if (restart)
  {
    embedded_fluid_->ReadRestart(restart);
  }

  if (ale_embfluid_) dispnpoldstate_ = Teuchos::rcp(new Epetra_Vector(*embedded_fluid_->Dispnp()));

  return;
}

void FLD::XFluidFluid::UseBlockMatrix(bool splitmatrix)
{
  // TODO: is it reasonable to init Block Matrix with npr > 0 when just blocks are assigned later?
  // Think about memory in this context

  // should we shift this creation to xfluidfluidstate-class?

  if (splitmatrix)
    xff_state_->xffluidsysmat_ =
        Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
            *XFluidFluidMapExtractor(), *XFluidFluidMapExtractor(), 108, false, true));
}

void FLD::XFluidFluid::SetXFluidFluidParams()
{
  Teuchos::ParameterList& params_xf_stab = params_->sublist("XFLUID DYNAMIC/STABILIZATION");

  // additional eos pressure stabilization on the elements of the embedded discretization,
  // that contribute to the interface
  xff_eos_pres_emb_layer_ =
      DRT::INPUT::IntegralValue<bool>(params_xf_stab, "XFF_EOS_PRES_EMB_LAYER");

  // whether an eigenvalue problem has to be solved to estimate Nitsche's parameter
  nitsche_evp_ =
      (DRT::INPUT::IntegralValue<INPAR::XFEM::ViscStab_TraceEstimate>(params_xf_stab,
           "VISC_STAB_TRACE_ESTIMATE") == INPAR::XFEM::ViscStab_TraceEstimate_eigenvalue);

  // get general XFEM/XFFSI specific parameters
  monolithic_approach_ = DRT::INPUT::IntegralValue<INPAR::XFEM::Monolithic_xffsi_Approach>(
      params_->sublist("XFLUID DYNAMIC/GENERAL"), "MONOLITHIC_XFFSI_APPROACH");
  xfem_timeintapproach_ = DRT::INPUT::IntegralValue<INPAR::XFEM::XFluidFluidTimeInt>(
      params_->sublist("XFLUID DYNAMIC/GENERAL"), "XFLUIDFLUID_TIMEINT");

  // get information about active shape derivatives
  active_shapederivatives_ = ale_embfluid_ && params_->get<bool>("shape derivatives");
}

void FLD::XFluidFluid::SetInitialFlowField(
    const INPAR::FLUID::InitialField initfield, const int startfuncno)
{
  XFluid::SetInitialFlowField(initfield, startfuncno);
  embedded_fluid_->SetInitialFlowField(initfield, startfuncno);
}

void FLD::XFluidFluid::SetInterfaceFixed() { mc_xff_->SetInterfaceFixed(); }

void FLD::XFluidFluid::SetInterfaceFree() { mc_xff_->SetInterfaceFree(); }

void FLD::XFluidFluid::PrepareTimeStep()
{
  embedded_fluid_->PrepareTimeStep();
  XFluid::PrepareTimeStep();
}

Teuchos::RCP<const Epetra_Vector> FLD::XFluidFluid::InitialGuess()
{
  xff_state_->xffluidsplitter_->InsertFluidVector(
      embedded_fluid_->InitialGuess(), xff_state_->xffluidincvel_);
  xff_state_->xffluidsplitter_->InsertXFluidVector(
      XFluid::InitialGuess(), xff_state_->xffluidincvel_);
  return xff_state_->xffluidincvel_;
}

void FLD::XFluidFluid::PrepareXFEMSolve()
{
  XFluid::PrepareXFEMSolve();

  // merge the velnp each into one large Epetra_Vector for the composed system
  xff_state_->xffluidsplitter_->InsertXFluidVector(xff_state_->velnp_, xff_state_->xffluidvelnp_);
  xff_state_->xffluidsplitter_->InsertFluidVector(
      embedded_fluid_->Velnp(), xff_state_->xffluidvelnp_);

  xff_state_->xffluidsplitter_->InsertXFluidVector(xff_state_->veln_, xff_state_->xffluidveln_);
  xff_state_->xffluidsplitter_->InsertFluidVector(
      embedded_fluid_->Veln(), xff_state_->xffluidveln_);
}

void FLD::XFluidFluid::Evaluate(
    Teuchos::RCP<const Epetra_Vector> stepinc  ///< solution increment between time step n and n+1
)
{
  // split step increment
  Teuchos::RCP<Epetra_Vector> stepinc_xfluid = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> stepinc_emb = Teuchos::null;

  if (stepinc != Teuchos::null)
  {
    stepinc_xfluid = xff_state_->xffluidsplitter_->ExtractXFluidVector(stepinc);
    stepinc_emb = xff_state_->xffluidsplitter_->ExtractFluidVector(stepinc);

    // compute increment
    xff_state_->xffluidincvel_->Update(1.0, *stepinc, -1.0, *stepinc_, 0.0);

    xff_state_->xffluiddbcmaps_->InsertCondVector(
        xff_state_->xffluiddbcmaps_->ExtractCondVector(xff_state_->xffluidzeros_),
        xff_state_->xffluidincvel_);

    // update embedded fluid solution by increment
    embedded_fluid_->UpdateIterIncrementally(
        xff_state_->xffluidsplitter_->ExtractFluidVector(xff_state_->xffluidincvel_));
  }

  if (mc_xff_->GetAveragingStrategy() == INPAR::XFEM::Embedded_Sided and nitsche_evp_)
  {
    mc_xff_->ResetEvaluatedTraceEstimates();
    if (ale_embfluid_)
      embedded_fluid_->Discretization()->SetState("dispnp", embedded_fluid_->Dispnp());
  }

  // evaluation of background fluid (new cut for full Newton approach)
  XFluid::UpdateByIncrements(stepinc_xfluid);
  XFluid::Evaluate();

  // update step increment
  stepinc_->Update(1.0, *xff_state_->xffluidvelnp_, -1.0, *xff_state_->xffluidveln_, 0.0);

  if (active_shapederivatives_)
  {
    extended_shapederivatives_->Add(*embedded_fluid_->ShapeDerivatives(), false, 1.0, 0.0);
  }

  xff_state_->xffluidincvel_->PutScalar(0.0);

  xff_state_->trueresidual_ = Teuchos::rcp(new Epetra_Vector(*xff_state_->xffluidresidual_));
  xff_state_->trueresidual_->PutScalar(ResidualScaling());
}

void FLD::XFluidFluid::TimeUpdate()
{
  embedded_fluid_->TimeUpdate();
  XFluid::TimeUpdate();
  xff_state_->xffluidveln_->Update(1.0, *xff_state_->xffluidvelnp_, 0.0);
}

Teuchos::RCP<LINALG::BlockSparseMatrixBase> FLD::XFluidFluid::BlockSystemMatrix(
    Teuchos::RCP<Epetra_Map> innermap, Teuchos::RCP<Epetra_Map> condmap)
{
  // Map of fluid FSI DOFs: condmap
  // Map of inner fluid DOFs: innermap

  // Get the fluid-fluid system matrix as sparse matrix
  Teuchos::RCP<LINALG::SparseMatrix> sparsesysmat = SystemMatrix();

  // F_{II}, F_{I\Gamma}, F_{\GammaI}, F_{\Gamma\Gamma}
  Teuchos::RCP<LINALG::SparseMatrix> fii, fig, fgi, fgg;
  // Split sparse system matrix into blocks according to the given maps
  LINALG::SplitMatrix2x2(sparsesysmat, innermap, condmap, innermap, condmap, fii, fig, fgi, fgg);
  // create a new block matrix out of the 4 blocks
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockmat =
      LINALG::BlockMatrix2x2(*fii, *fig, *fgi, *fgg);

  if (blockmat == Teuchos::null) dserror("Creation of fluid-fluid block matrix failed.");

  return blockmat;
}

Teuchos::RCP<const Epetra_Map> FLD::XFluidFluid::PressureRowMap()
{
  return xff_state_->xffluidvelpressplitter_->CondMap();
}

Teuchos::RCP<const Epetra_Map> FLD::XFluidFluid::VelocityRowMap()
{
  return xff_state_->xffluidvelpressplitter_->OtherMap();
}

Teuchos::RCP<DRT::ResultTest> FLD::XFluidFluid::CreateFieldTest()
{
  return Teuchos::rcp(new FLD::XFluidResultTest(*this));
}

Teuchos::RCP<FLD::XFluidState> FLD::XFluidFluid::GetNewState()
{
  // further use type-cast pointer to MeshCouplingFluidFluid
  if (mc_xff_ == Teuchos::null)
  {
    mc_xff_ = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFluidFluid>(
        condition_manager_->GetMeshCoupling(cond_name_));

    if (mc_xff_ == Teuchos::null) dserror("Failed to cast to MeshCouplingFluidFluid");
  }

  if (ale_embfluid_)
  {
    mc_xff_->UpdateDisplacementIterationVectors();  // update last iteration interface displacements
    LINALG::Export(*embedded_fluid_->Dispnp(), *mc_xff_->IDispnp());
  }

  state_it_++;

  Teuchos::RCP<FLD::XFluidFluidState> state =
      state_creator_->Create(xdiscret_, embedded_fluid_->Discretization(),
          Teuchos::null,  //!< col vector holding background ALE displacements for backdis
          solver_->Params(), step_, time_);

  // increment vector for merged background & embedded fluid
  // (not the classical Newton increment but the difference to
  // the value at the last time step)
  stepinc_ = LINALG::CreateVector(*state->xffluiddofrowmap_, true);

  // build a merged map from fluid-fluid dbc-maps
  state->CreateMergedDBCMapExtractor(embedded_fluid_->GetDBCMapExtractor());

  return state;
}

void FLD::XFluidFluid::CreateState()
{
  // free the pointer to state_ object to enable to destroy the state_ object
  xff_state_ = Teuchos::null;

  // new cut for this time step
  XFluid::CreateState();
  xff_state_ = Teuchos::rcp_dynamic_cast<FLD::XFluidFluidState>(XFluid::state_);

  if (xff_state_ == Teuchos::null) dserror("Failed to create an instance of XFluidFluidState.");

  if (!ale_embfluid_) return;
}

void FLD::XFluidFluid::AssembleMatAndRHS(int itnum  ///< iteration number
)
{
  TEUCHOS_FUNC_TIME_MONITOR("FLD::XFluidFluid::AssembleMatAndRHS");

  // evaluate elements of embedded fluid
  embedded_fluid_->PrepareSolve();

  if (mc_xff_->GetAveragingStrategy() == INPAR::XFEM::Embedded_Sided)
  {
    mc_xff_->GetCouplingDis()->ClearState();
    // set velocity and displacement state for embedded fluid
    embedded_fluid_->SetStateTimInt();
    mc_xff_->GetCouplingDis()->SetState("veln", embedded_fluid_->Veln());

    if (ale_embfluid_)
    {
      mc_xff_->GetCouplingDis()->SetState("dispnp", embedded_fluid_->Dispnp());
    }
  }

  // export interface velocities
  // TODO: shift to mesh coupling class
  LINALG::Export(*(embedded_fluid_->Velnp()), *(mc_xff_->IVelnp()));
  LINALG::Export(*(embedded_fluid_->Veln()), *(mc_xff_->IVeln()));

  // evaluate elements of XFluid part
  XFluid::AssembleMatAndRHS(itnum);

  // insert XFluid residual to merged
  xff_state_->xffluidsplitter_->InsertXFluidVector(
      xff_state_->residual_, xff_state_->xffluidresidual_);

  // add coupling contribution to embedded residual
  const int mc_idx = condition_manager_->GetMeshCouplingIndex(cond_name_);
  Teuchos::RCP<XFluidState::CouplingState>& coup_state = xff_state_->coup_state_[mc_idx];

  {
    // adding rhC_s_ (coupling contribution) to residual of embedded fluid
    for (int iter = 0; iter < coup_state->rhC_s_->MyLength(); ++iter)
    {
      const int rhsgid = coup_state->rhC_s_->Map().GID(iter);
      if (coup_state->rhC_s_->Map().MyGID(rhsgid) == false)
        dserror("rhC_s_ should be on all processors");
      if (embedded_fluid_->Residual()->Map().MyGID(rhsgid))
        (*embedded_fluid_->Residual())[embedded_fluid_->Residual()->Map().LID(rhsgid)] +=
            (*coup_state->rhC_s_)[coup_state->rhC_s_->Map().LID(rhsgid)];
      else
        dserror("Interface dof %d does not belong to embedded discretization!", rhsgid);
    }
  }

  // add additional EOS-pressure stabilization to interface-contributing
  // layer of embedded fluid
  if (xff_eos_pres_emb_layer_)
  {
    AddEosPresStabToEmbLayer();
  }

  // add embedded part of merged residual
  xff_state_->xffluidsplitter_->InsertFluidVector(
      embedded_fluid_->Residual(), xff_state_->xffluidresidual_);

  // assemble XFluid and embedded fluid system matrices into one

  // TODO: when creation is shifted state-class, we can ask the state class for this
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat_block =
      Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(xff_state_->xffluidsysmat_, false);
  if (sysmat_block != Teuchos::null)
  {
    sysmat_block->Assign(1, 1, LINALG::View, *xff_state_->sysmat_);
    sysmat_block->Assign(1, 0, LINALG::View, *coup_state->C_xs_);
    sysmat_block->Assign(0, 1, LINALG::View, *coup_state->C_sx_);
    embedded_fluid_->SystemMatrix()->UnComplete();
    embedded_fluid_->SystemMatrix()->Add(*coup_state->C_ss_, false, 1.0, 1.0);
    Teuchos::RCP<LINALG::SparseMatrix> alesysmat_sparse =
        Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(embedded_fluid_->SystemMatrix());
    sysmat_block->Assign(0, 0, LINALG::View, *alesysmat_sparse);
  }
  else
  {
    // TODO introduce a ZeroFluidFluidSysmat in xfluidfluid state (use PutScalar if
    // explicitDirichlet = false)
    xff_state_->xffluidsysmat_->Zero();
    xff_state_->xffluidsysmat_->Add(*xff_state_->sysmat_, false, 1.0, 0.0);
    xff_state_->xffluidsysmat_->Add(*embedded_fluid_->SystemMatrix(), false, 1.0, 1.0);
    xff_state_->xffluidsysmat_->Add(*coup_state->C_xs_, false, 1.0, 1.0);
    xff_state_->xffluidsysmat_->Add(*coup_state->C_sx_, false, 1.0, 1.0);
    xff_state_->xffluidsysmat_->Add(*coup_state->C_ss_, false, 1.0, 1.0);
  }

  xff_state_->xffluidsysmat_->Complete();
}

void FLD::XFluidFluid::PrepareShapeDerivatives(
    const Teuchos::RCP<const LINALG::MultiMapExtractor> fsiextractor,
    const Teuchos::RCP<std::set<int>> condelements)
{
  if (!active_shapederivatives_) return;

  // here we initialize the shapederivates
  // REMARK: the shape derivatives matrix results from linearization w.r.t. ALE-displacements
  // and therefore solely knows ALE-dof - here we use "extended shapederivatives" including
  // background fluid entries, that are set to zero
  Teuchos::RCP<LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>> mat =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<FLD::UTILS::InterfaceSplitStrategy>(
          *fsiextractor, *fsiextractor, 108, false, true));
  mat->SetCondElements(condelements);
  extended_shapederivatives_ = mat;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluidFluid::UpdateByIncrement()
{
  // Update merged
  XFluid::UpdateByIncrement();
  // update xfluid
  xff_state_->velnp_->Update(
      1.0, *xff_state_->xffluidsplitter_->ExtractXFluidVector(xff_state_->xffluidvelnp_), 0.0);
  // update embedded fluid
  // embedded_fluid_->IterUpdate(xff_state_->xffluidsplitter_->ExtractFluidVector(xff_state_->xffluidincvel_));
  embedded_fluid_->WriteAccessVelnp()->Update(
      1.0, *xff_state_->xffluidsplitter_->ExtractFluidVector(xff_state_->xffluidvelnp_), 0.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluidFluid::AddEosPresStabToEmbLayer()
{
  if (ale_embfluid_)
    embedded_fluid_->Discretization()->SetState("gridv", embedded_fluid_->GridVel());

  Teuchos::ParameterList faceparams;

  const Teuchos::RCP<DRT::DiscretizationFaces> xdiscret =
      Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(embedded_fluid_->Discretization(), true);

  // set additional faceparams according to ghost-penalty terms due to Nitsche's method
  faceparams.set("ghost_penalty_reconstruct", false);

  //------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> residual_col = LINALG::CreateVector(*xdiscret->DofColMap(), true);

  //------------------------------------------------------------
  const Epetra_Map* rmap = NULL;

  // TODO: do not create a new matrix all the time, why not creating an epetraFE matrix in
  // fluidimplicit directly?
  Teuchos::RCP<Epetra_FECrsMatrix> sysmat_FE;

  rmap = &(embedded_fluid_->SystemMatrix()->OperatorRangeMap());
  sysmat_FE = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy, *rmap, 256, false));

  // TODO: think about the dirichlet and savegraph flags when ApplyDirichlet or Zero is called
  Teuchos::RCP<LINALG::SparseMatrix> sysmat_linalg =
      Teuchos::rcp(new LINALG::SparseMatrix(Teuchos::rcp_static_cast<Epetra_CrsMatrix>(sysmat_FE),
          LINALG::View, true, true, LINALG::SparseMatrix::FE_MATRIX));

  //------------------------------------------------------------
  // loop over row faces

  const int numrowintfaces = xdiscret->NumMyRowFaces();

  for (int i = 0; i < numrowintfaces; ++i)
  {
    DRT::Element* actface = xdiscret->lRowFace(i);
    DRT::ELEMENTS::FluidIntFace* ele = dynamic_cast<DRT::ELEMENTS::FluidIntFace*>(actface);
    if (ele == NULL) dserror("expect FluidIntFace element");
    edgestab_->EvaluateEdgeStabBoundaryGP(faceparams, xdiscret,
        mc_xff_->GetAuxiliaryDiscretization(), ele, sysmat_linalg, residual_col);
  }

  //------------------------------------------------------------
  sysmat_linalg->Complete();
  embedded_fluid_->SystemMatrix()->UnComplete();

  (embedded_fluid_->SystemMatrix())->Add(*sysmat_linalg, false, 1.0, 1.0);
  embedded_fluid_->SystemMatrix()->Complete();
  //------------------------------------------------------------
  // need to export residual_col to embedded fluid residual
  {
    Epetra_Vector res_tmp(embedded_fluid_->Residual()->Map(), true);
    Epetra_Export exporter(residual_col->Map(), res_tmp.Map());
    int err = res_tmp.Export(*residual_col, exporter, Add);
    if (err) dserror("Export using exporter returned err=%d", err);
    embedded_fluid_->Residual()->Update(1.0, res_tmp, 1.0);
  }

  mc_xff_->GetCouplingDis()->ClearState();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool FLD::XFluidFluid::XTimint_ProjectFromEmbeddedDiscretization(
    const Teuchos::RCP<XFEM::XFluidTimeInt>& xfluid_timeint,  ///< xfluid time integration class
    std::vector<Teuchos::RCP<Epetra_Vector>>& newRowStateVectors,  ///< vectors to be reconstructed
    Teuchos::RCP<const Epetra_Vector> target_dispnp,  ///< displacement col - vector timestep n+1
    const bool screen_out                             ///< screen output?
)
{
  std::vector<Teuchos::RCP<const Epetra_Vector>> oldStateVectors;

  Teuchos::RCP<const Epetra_Vector> velncol = DRT::UTILS::GetColVersionOfRowVector(
      embedded_fluid_->Discretization(), embedded_fluid_->Veln());
  Teuchos::RCP<const Epetra_Vector> accncol = DRT::UTILS::GetColVersionOfRowVector(
      embedded_fluid_->Discretization(), embedded_fluid_->Accn());
  oldStateVectors.push_back(velncol);
  oldStateVectors.push_back(accncol);

  // get set of node-ids, that demand projection from embedded discretization
  std::map<int, std::set<int>>& projection_nodeToDof =
      xfluid_timeint->Get_NodeToDofMap_For_Reconstr(INPAR::XFEM::Xf_TimeInt_by_PROJ_from_DIS);

  Teuchos::RCP<const Epetra_Vector> disp =
      DRT::UTILS::GetColVersionOfRowVector(embedded_fluid_->Discretization(), dispnpoldstate_);
  projector_->SetSourcePositionVector(disp);
  projector_->SetSourceStateVectors(oldStateVectors);

  projector_->Project(projection_nodeToDof, newRowStateVectors, target_dispnp);

  int numfailed = 0;
  int my_numfailed = projection_nodeToDof.size();
  discret_->Comm().SumAll(&my_numfailed, &numfailed, 1);

  return numfailed == 0;

  // projector_->GmshOutput(step_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool FLD::XFluidFluid::XTimint_DoIncrementStepTransfer(
    const bool screen_out, const bool firstcall_in_timestep)
{
  // use increment step transfer if :
  // - at least one XFEM interface is moving (the fluid-fluid interface or any other)
  // - first call in time step to initialize velnp
  if (mc_xff_->HasMovingInterface() || meshcoupl_dis_.size() > 1 || firstcall_in_timestep)
  {
    return XFluid::XTimint_DoIncrementStepTransfer(screen_out, firstcall_in_timestep);
  }
  else
    return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::XFluidFluid::Output()
{
  XFluid::Output();
  embedded_fluid_->Output();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::XFluidFluid::UpdateMonolithicFluidSolution(
    const Teuchos::RCP<const Epetra_Map>& fsidofmap)
{
  // manipulate the dbc map extractor
  Teuchos::RCP<const Epetra_Map> dbcmap = embedded_fluid_->GetDBCMapExtractor()->CondMap();
  std::vector<Teuchos::RCP<const Epetra_Map>> condmaps;
  condmaps.push_back(dbcmap);
  condmaps.push_back(fsidofmap);
  Teuchos::RCP<Epetra_Map> condmerged = LINALG::MultiMapExtractor::MergeMaps(condmaps);

  Teuchos::RCP<LINALG::MapExtractor> fsidbcmapex =
      Teuchos::rcp(new LINALG::MapExtractor(*(embedded_fluid_->DofRowMap()), condmerged));

  // DBC map-extractor containing FSI-dof
  xff_state_->CreateMergedDBCMapExtractor(fsidbcmapex);
  Solve();
  xff_state_->CreateMergedDBCMapExtractor(embedded_fluid_->GetDBCMapExtractor());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::XFluidFluid::InterpolateEmbeddedStateVectors()
{
  XFEM::MeshProjector embedded_projector(
      embedded_fluid_->Discretization(), embedded_fluid_->Discretization(), *params_);
  std::vector<Teuchos::RCP<Epetra_Vector>> newRowStateVectors;

  newRowStateVectors.push_back(embedded_fluid_->WriteAccessVelnp());
  newRowStateVectors.push_back(embedded_fluid_->WriteAccessAccnp());

  std::vector<Teuchos::RCP<const Epetra_Vector>> oldStateVectors;

  Teuchos::RCP<const Epetra_Vector> velncol = DRT::UTILS::GetColVersionOfRowVector(
      embedded_fluid_->Discretization(), embedded_fluid_->Velnp());
  Teuchos::RCP<const Epetra_Vector> accncol = DRT::UTILS::GetColVersionOfRowVector(
      embedded_fluid_->Discretization(), embedded_fluid_->Accnp());
  oldStateVectors.push_back(velncol);
  oldStateVectors.push_back(accncol);

  Teuchos::RCP<const Epetra_Vector> srcdisp =
      DRT::UTILS::GetColVersionOfRowVector(embedded_fluid_->Discretization(), dispnpoldstate_);
  embedded_projector.SetSourcePositionVector(srcdisp);
  embedded_projector.SetSourceStateVectors(oldStateVectors);

  Teuchos::RCP<const Epetra_Vector> tardisp = DRT::UTILS::GetColVersionOfRowVector(
      embedded_fluid_->Discretization(), embedded_fluid_->Dispnp());

  embedded_projector.ProjectInFullTargetDiscretization(newRowStateVectors, tardisp);

  // embedded_projector.GmshOutput(step_,embedded_fluid_->Dispnp());

  embedded_fluid_->SetOldPartOfRighthandside();
  embedded_fluid_->SetDirichletNeumannBC();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<std::vector<double>> FLD::XFluidFluid::EvaluateErrorComparedToAnalyticalSol()
{
  // this function provides a general implementation for calculating error norms between computed
  // solutions and an analytical solution which is implemented or given by a function in the input
  // file

  INPAR::FLUID::CalcError calcerr =
      DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_, "calculate error");

  if (calcerr == INPAR::FLUID::no_error_calculation) return Teuchos::null;
  // set the time to evaluate errors
  //

  // define the norms that have to be computed

  //-------------------------------------------------------------------------------------------------------------------
  // domain error norms w.r.t incompressible Navier-Stokes equations
  //
  //--------------------------------------
  // background domain
  //--------------------------------------
  //
  // standard domain errors
  // 1.   || u - u_h ||_L2(Omega)              =   standard L2-norm for velocity
  // 2.   || grad( u - u_h ) ||_L2(Omega)      =   standard H1-seminorm for velocity
  // 3.   || u - u_h ||_H1(Omega)              =   standard H1-norm for velocity
  //                                           =   sqrt( || u - u_h ||^2_L2(Omega) + || grad( u -
  //                                           u_h ) ||^2_L2(Omega) )
  // 4.   || p - p_h ||_L2(Omega)              =   standard L2-norm for pressure
  //
  // viscosity-scaled domain errors
  // 5.   || nu^(+1/2) grad( u - u_h ) ||_L2(Omega)      =   visc-scaled H1-seminorm for velocity
  //                                                     =   nu^(+1/2) * || grad( u - u_h )
  //                                                     ||_L2(Omega) (for homogeneous visc)
  // 6.   || nu^(-1/2) (p - p_h) ||_L2(Omega)            =   visc-scaled L2-norm for pressure
  //                                                     =   nu^(-1/2) * || p - p_h ||_L2(Omega)
  //                                                     (for homogeneous visc)
  // 7.   || sigma^(+1/2) ( u - u_h ) ||_L2(Omega)       =   sigma-scaled L2-norm for velocity
  //                                                     =   sigma^(+1/2) * || u - u_h ||_L2(Omega)
  //                                                     (for homogeneous sigma)
  // 8.   || Phi^(+1/2) (p - p_h) ||_L2(Omega)           =   Phi-scaled L2-norm for pressure
  //                                                     =   Phi^(+1/2) * || p - p_h ||_L2(Omega)
  //                                                     (for homogeneous Phi)
  // with Phi^{-1} = sigma*CP^2 + |beta|*CP + nu + (|beta|*CP/sqrt(sigma*CP^2 + nu))^2, see
  // Massing,Schott,Wall Oseen paper
  //
  // 9. functional G=sin(x)( u,x - u,x exact ) (Sudhakar)
  //
  //
  //--------------------------------------
  // embedded domain
  //--------------------------------------
  //-------------------------------------------------------------------------------------------------------------------
  // domain error norms w.r.t incompressible Navier-Stokes/ Oseen equations
  //
  // standard domain errors
  // 1.   || u - u_h ||_L2(Omega)              =   standard L2-norm for velocity
  // 2.   || grad( u - u_h ) ||_L2(Omega)      =   standard H1-seminorm for velocity
  // 3.   || u - u_h ||_H1(Omega)              =   standard H1-norm for velocity
  //                                           =   sqrt( || u - u_h ||^2_L2(Omega) + || grad( u -
  //                                           u_h ) ||^2_L2(Omega) )
  // 4.   || p - p_h ||_L2(Omega)              =   standard L2-norm for pressure
  //
  // viscosity-scaled domain errors
  // 5.   || nu^(+1/2) grad( u - u_h ) ||_L2(Omega)      =   visc-scaled H1-seminorm for velocity
  //                                                     =   nu^(+1/2) * || grad( u - u_h )
  //                                                     ||_L2(Omega) (for homogeneous visc)
  // 6.   || nu^(-1/2) (p - p_h) ||_L2(Omega)            =   visc-scaled L2-norm for pressure
  //                                                     =   nu^(-1/2) * || p - p_h ||_L2(Omega)
  //                                                     (for homogeneous visc)
  // 7.   || sigma^(+1/2) ( u - u_h ) ||_L2(Omega)       =   sigma-scaled L2-norm for velocity
  //                                                     =   sigma^(+1/2) * || u - u_h ||_L2(Omega)
  //                                                     (for homogeneous sigma)
  // 8.   || Phi^(+1/2) (p - p_h) ||_L2(Omega)           =   Phi-scaled L2-norm for pressure
  //                                                     =   Phi^(+1/2) * || p - p_h ||_L2(Omega)
  //                                                     (for homogeneous Phi)
  // with Phi^{-1} = sigma*CP^2 + |beta|*CP + nu + (|beta|*CP/sqrt(sigma*CP^2 + nu))^2, see
  // Massing,Schott,Wall Oseen paper
  //
  // 9. functional G=sin(x)( u,x - u,x exact ) (Sudhakar)
  //-------------------------------------------------------------------------------------------------------------------
  // interface/boundary error norms at the XFEM-interface, boundary
  // w.r.t Nitsche's method to enforce interface/boundary conditions
  //
  // 1.   || nu^(+1/2) (u - u*) ||_H1/2(Gamma)             =  broken H1/2 Sobolev norm for
  // boundary/coupling condition
  // 2.   || nu^(+1/2) grad( u - u_h )*n ||_H-1/2(Gamma)   =  standard H-1/2 Sobolev norm for normal
  // flux (velocity part)
  // 3.   || nu^(-1/2) (p - p_h)*n ||_H-1/2(Gamma)         =  standard H-1/2 Sobolev norm for normal
  // flux (pressure part)
  // 4.   || (u*n)_inflow (u - u*) ||_L2(Gamma)            =  L^2 Sobolev norm for inflow
  // boundary/coupling condition
  // 5.   || (sigma*h+|u|+nu/h)^(+1/2) (u - u*)*n ||_L2(Gamma) =  L^2 Sobolev norm for mass
  // conservation coupling condition
  //
  //-------------------------------------------------------------------------------------------------------------------
  // errors introduced by stabilizations (edge-based fluid stabilizations and ghost-penalty
  // stabilizations)
  //
  // ...
  //-------------------------------------------------------------------------------------------------------------------

  // number of norms that have to be calculated
  const int num_dom_norms = 10;
  const int num_interf_norms = 8;
  const int num_stab_norms = 3;

  Epetra_SerialDenseVector cpu_dom_norms(num_dom_norms);
  Epetra_SerialDenseVector cpu_dom_norms_emb(num_dom_norms);
  Epetra_SerialDenseVector cpu_interf_norms(num_interf_norms);
  Epetra_SerialDenseVector cpu_stab_norms(num_stab_norms);

  Teuchos::RCP<Epetra_SerialDenseVector> glob_dom_norms_bg =
      Teuchos::rcp(new Epetra_SerialDenseVector(num_dom_norms));
  Teuchos::RCP<Epetra_SerialDenseVector> glob_dom_norms_emb =
      Teuchos::rcp(new Epetra_SerialDenseVector(num_dom_norms));
  Teuchos::RCP<Epetra_SerialDenseVector> glob_interf_norms =
      Teuchos::rcp(new Epetra_SerialDenseVector(num_interf_norms));
  Teuchos::RCP<Epetra_SerialDenseVector> glob_stab_norms =
      Teuchos::rcp(new Epetra_SerialDenseVector(num_stab_norms));

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("u and p at time n+1 (converged)", state_->velnp_);

  mc_xff_->GetCondDis()->ClearState();
  mc_xff_->GetCondDis()->SetState("velaf", embedded_fluid_->Velnp());
  // mc_xff_->GetCondDis()->SetState("dispnp", embedded_fluid_->Dispnp());

  mc_xff_->SetState();

  // evaluate domain error norms and interface/boundary error norms at XFEM-interface
  // loop row elements of background fluid
  XFluid::ComputeErrorNorms(glob_dom_norms_bg, glob_interf_norms, glob_stab_norms);

  //-----------------------------------------------
  // Embedded discretization
  //---------------------------------------------
  // set vector values needed by elements
  mc_xff_->GetCondDis()->ClearState();
  mc_xff_->GetCondDis()->SetState("u and p at time n+1 (converged)", embedded_fluid_->Velnp());

  // evaluate domain error norms and interface/boundary error norms at XFEM-interface
  // loop row elements
  const int numrowele_emb = mc_xff_->GetCondDis()->NumMyRowElements();
  for (int i = 0; i < numrowele_emb; ++i)
  {
    // local element-wise squared error norms
    Epetra_SerialDenseVector ele_dom_norms_emb(num_dom_norms);

    // pointer to current element
    DRT::Element* actele = mc_xff_->GetCondDis()->lRowElement(i);

    Teuchos::RCP<MAT::Material> mat = actele->Material();

    DRT::ELEMENTS::Fluid* ele = dynamic_cast<DRT::ELEMENTS::Fluid*>(actele);

    DRT::Element::LocationArray la(1);

    // get element location vector, dirichlet flags and ownerships
    actele->LocationVector(*mc_xff_->GetCondDis(), la, false);

    Epetra_SerialDenseMatrix elemat1;
    Epetra_SerialDenseMatrix elemat2;
    Epetra_SerialDenseVector elevec2;
    Epetra_SerialDenseVector elevec3;
    params_->set<int>("action", FLD::calc_fluid_error);

    DRT::ELEMENTS::FluidFactory::ProvideImplXFEM(actele->Shape(), "xfem")
        ->EvaluateService(ele, *params_, mat, *mc_xff_->GetCondDis(), la[0].lm_, elemat1, elemat2,
            ele_dom_norms_emb, elevec2, elevec3);

    // sum up (on each processor)
    cpu_dom_norms_emb += ele_dom_norms_emb;

  }  // end loop over embedded fluid elements

  //--------------------------------------------------------
  // reduce and sum over all procs

  for (int i = 0; i < num_dom_norms; ++i) (*glob_dom_norms_emb)(i) = 0.0;
  mc_xff_->GetCondDis()->Comm().SumAll(
      cpu_dom_norms_emb.Values(), glob_dom_norms_emb->Values(), num_dom_norms);

  // standard domain errors bg-dis
  double dom_bg_err_vel_L2 =
      0.0;  //  || u - u_b ||_L2(Omega)           =   standard L2-norm for velocity
  double dom_bg_err_vel_H1_semi =
      0.0;  //  || grad( u - u_b ) ||_L2(Omega)   =   standard H1-seminorm for velocity
  double dom_bg_err_vel_H1 =
      0.0;  //  || u - u_b ||_H1(Omega)           =   standard H1-norm for velocity
  double dom_bg_err_pre_L2 =
      0.0;  //  || p - p_b ||_L2(Omega)           =   standard L2-norm for pressure

  // viscosity-scaled domain errors
  double dom_bg_err_vel_H1_semi_nu_scaled =
      0.0;  //  || nu^(+1/2) grad( u - u_b ) ||_L2(Omega)  =   visc-scaled H1-seminorm for velocity
  double dom_bg_err_pre_L2_nu_scaled =
      0.0;  //  || nu^(-1/2) (p - p_b) ||_L2(Omega)        =   visc-scaled L2-norm for pressure
  double dom_bg_err_vel_L2_sigma_scaled =
      0.0;  //  || sigma^(+1/2) ( u - u_h ) ||_L2(Omega)       =   sigma-scaled L2-norm for velocity
  double dom_bg_err_pre_L2_Phi_scaled =
      0.0;  //  || Phi^(+1/2) (p - p_h) ||_L2(Omega)           =   Phi-scaled L2-norm for pressure

  // standard domain errors bg-dis
  double dom_emb_err_vel_L2 =
      0.0;  //  || u - u_e ||_L2(Omega)           =   standard L2-norm for velocity
  double dom_emb_err_vel_H1_semi =
      0.0;  //  || grad( u - u_e ) ||_L2(Omega)   =   standard H1-seminorm for velocity
  double dom_emb_err_vel_H1 =
      0.0;  //  || u - u_e ||_H1(Omega)           =   standard H1-norm for velocity
  double dom_emb_err_pre_L2 =
      0.0;  //  || p - p_e ||_L2(Omega)           =   standard L2-norm for pressure

  // viscosity-scaled domain errors
  double dom_emb_err_vel_H1_semi_nu_scaled =
      0.0;  //  || nu^(+1/2) grad( u - u_e ) ||_L2(Omega)  =   visc-scaled H1-seminorm for velocity
  double dom_emb_err_pre_L2_nu_scaled =
      0.0;  //  || nu^(-1/2) (p - p_e) ||_L2(Omega)        =   visc-scaled L2-norm for pressure
  double dom_emb_err_vel_L2_sigma_scaled =
      0.0;  //  || sigma^(+1/2) ( u - u_h ) ||_L2(Omega)       =   sigma-scaled L2-norm for velocity
  double dom_emb_err_pre_L2_Phi_scaled =
      0.0;  //  || Phi^(+1/2) (p - p_h) ||_L2(Omega)           =   Phi-scaled L2-norm for pressure

  // interface errors
  double interf_err_Honehalf = 0.0;  //  || nu^(+1/2) (u_b - u_e) ||_H1/2(Gamma)          =  broken
                                     //  H1/2 Sobolev norm for boundary/coupling condition
  double interf_err_Hmonehalf_u =
      0.0;  //  || nu^(+1/2) grad( u_b - u_e )*n ||_H-1/2(Gamma) =  broken H-1/2 Sobolev norm for
            //  normal flux (velocity part)
  double interf_err_Hmonehalf_p =
      0.0;  //  || nu^(-1/2) (p_b - p_e)*n ||_H-1/2(Gamma)         =  broken H-1/2 Sobolev norm for
            //  normal flux (pressure part)
  double interf_err_inflow = 0.0;     //  || (u*n)_inflow (u - u*) ||_L2(Gamma)            =  L^2
                                      //  Sobolev norm for inflow boundary/coupling condition
  double interf_err_mass_cons = 0.0;  //  || (sigma*h+|u|+nu/h)^(+1/2) (u - u*)*n ||_L2(Gamma) = L^2
                                      //  Sobolev norm for mass conservation coupling condition

  // sudhakar functional for testing integration
  double functional_bg = 0.0;
  double functional_emb = 0.0;

  dom_bg_err_vel_L2 = sqrt((*glob_dom_norms_bg)[0]);
  dom_bg_err_vel_H1_semi = sqrt((*glob_dom_norms_bg)[1]);
  dom_bg_err_vel_H1 = sqrt((*glob_dom_norms_bg)[2]);
  dom_bg_err_pre_L2 = sqrt((*glob_dom_norms_bg)[3]);

  dom_bg_err_vel_H1_semi_nu_scaled = sqrt((*glob_dom_norms_bg)[4]);
  dom_bg_err_pre_L2_nu_scaled = sqrt((*glob_dom_norms_bg)[5]);
  dom_bg_err_vel_L2_sigma_scaled = sqrt((*glob_dom_norms_bg)[6]);
  dom_bg_err_pre_L2_Phi_scaled = sqrt((*glob_dom_norms_bg)[7]);

  functional_bg = (*glob_dom_norms_bg)[8];

  dom_emb_err_vel_L2 = sqrt((*glob_dom_norms_emb)[0]);
  dom_emb_err_vel_H1_semi = sqrt((*glob_dom_norms_emb)[1]);
  dom_emb_err_vel_H1 = sqrt((*glob_dom_norms_emb)[2]);
  dom_emb_err_pre_L2 = sqrt((*glob_dom_norms_emb)[3]);

  dom_emb_err_vel_H1_semi_nu_scaled = sqrt((*glob_dom_norms_emb)[4]);
  dom_emb_err_pre_L2_nu_scaled = sqrt((*glob_dom_norms_emb)[5]);
  dom_emb_err_vel_L2_sigma_scaled = sqrt((*glob_dom_norms_emb)[6]);
  dom_emb_err_pre_L2_Phi_scaled = sqrt((*glob_dom_norms_emb)[7]);

  functional_emb = (*glob_dom_norms_emb)[8];


  interf_err_Honehalf = sqrt((*glob_interf_norms)[0]);
  interf_err_Hmonehalf_u = sqrt((*glob_interf_norms)[1]);
  interf_err_Hmonehalf_p = sqrt((*glob_interf_norms)[2]);
  interf_err_inflow = sqrt((*glob_interf_norms)[3]);
  interf_err_mass_cons = sqrt((*glob_interf_norms)[4]);

  if (myrank_ == 0)
  {
    {
      std::cout.precision(8);
      IO::cout << IO::endl
               << "---- error norm for analytical solution Nr. "
               << DRT::INPUT::get<INPAR::FLUID::CalcError>(*params_, "calculate error")
               << " ----------" << IO::endl;
      IO::cout << "-------------- domain error norms (background)------------" << IO::endl;
      IO::cout << "|| u - u_b ||_L2(Omega)                        =  " << dom_bg_err_vel_L2
               << IO::endl;
      IO::cout << "|| grad( u - u_b ) ||_L2(Omega)                =  " << dom_bg_err_vel_H1_semi
               << IO::endl;
      IO::cout << "|| u - u_b ||_H1(Omega)                        =  " << dom_bg_err_vel_H1
               << IO::endl;
      IO::cout << "|| p - p_b ||_L2(Omega)                        =  " << dom_bg_err_pre_L2
               << IO::endl;
      IO::cout << "-------------- domain error norms (embedded)  ------------" << IO::endl;
      IO::cout << "|| u - u_e ||_L2(Omega)                        =  " << dom_emb_err_vel_L2
               << IO::endl;
      IO::cout << "|| grad( u_ - u_h ) ||_L2(Omega)               =  " << dom_emb_err_vel_H1_semi
               << IO::endl;
      IO::cout << "|| u - u_e ||_H1(Omega)                        =  " << dom_emb_err_vel_H1
               << IO::endl;
      IO::cout << "|| p - p_e ||_L2(Omega)                        =  " << dom_emb_err_pre_L2
               << IO::endl;
      IO::cout << "----viscosity-scaled domain error norms (background)------" << IO::endl;
      IO::cout << "|| nu^(+1/2) grad( u - u_b ) ||_L2(Omega)      =  "
               << dom_bg_err_vel_H1_semi_nu_scaled << IO::endl;
      IO::cout << "|| nu^(-1/2) (p - p_b) ||_L2(Omega)            =  "
               << dom_bg_err_pre_L2_nu_scaled << IO::endl;
      IO::cout << "|| sigma^(+1/2) ( u - u_h ) ||_L2(Omega)       =  "
               << dom_bg_err_vel_L2_sigma_scaled << IO::endl;
      IO::cout << "|| Phi^(+1/2) (p - p_h) ||_L2(Omega)           =  "
               << dom_bg_err_pre_L2_Phi_scaled << IO::endl;
      IO::cout << "----viscosity-scaled domain error norms (embedded) ------" << IO::endl;
      IO::cout << "|| nu^(+1/2) grad( u - u_e ) ||_L2(Omega)      =  "
               << dom_emb_err_vel_H1_semi_nu_scaled << IO::endl;
      IO::cout << "|| nu^(-1/2) (p - p_e) ||_L2(Omega)            =  "
               << dom_emb_err_pre_L2_nu_scaled << IO::endl;
      IO::cout << "|| sigma^(+1/2) ( u - u_h ) ||_L2(Omega)       =  "
               << dom_emb_err_vel_L2_sigma_scaled << IO::endl;
      IO::cout << "|| Phi^(+1/2) (p - p_h) ||_L2(Omega)           =  "
               << dom_emb_err_pre_L2_Phi_scaled << IO::endl;
      IO::cout << "---------------------------------------------------------" << IO::endl;
      IO::cout << "-------------- interface/boundary error norms -----------" << IO::endl;
      IO::cout << "|| nu^(+1/2) (u_b - u_e) ||_H1/2(Gamma)            =  " << interf_err_Honehalf
               << IO::endl;
      IO::cout << "|| nu^(+1/2) grad( u_b - u_e )*n ||_H-1/2(Gamma)   =  " << interf_err_Hmonehalf_u
               << IO::endl;
      IO::cout << "|| nu^(-1/2) (p_b - p_e)*n ||_H-1/2(Gamma)         =  " << interf_err_Hmonehalf_p
               << IO::endl;
      IO::cout << "|| (u*n)_inflow (u_b - u_e) ||_L2(Gamma)           =  " << interf_err_inflow
               << IO::endl;
      IO::cout << "|| (sigma*h+|u|+nu/h)^(+1/2) (u_b - u_e)*n ||_L2(Gamma)  =  "
               << interf_err_mass_cons << IO::endl;
      IO::cout << "---------------------------------------------------------" << IO::endl;
      IO::cout << "-------------- Error on Functionals from solution  ------------" << IO::endl;
      IO::cout << " | sin(x) ( u,x - u,x exact ) | (background)        = " << functional_bg
               << IO::endl;
      IO::cout << " | sin(x) ( u,x - u,x exact ) | (embedded)          = " << functional_emb
               << IO::endl;
    }

    // append error of the last time step to the error file
    if ((step_ == stepmax_) or (time_ == maxtime_))  // write results to file
    {
      std::ostringstream temp;
      const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
      const std::string fname = simulation + ".xfem_abserror";

      std::ofstream f;
      f.open(fname.c_str(), std::fstream::ate | std::fstream::app);
      f << "#| " << simulation << "\n";
      f << "#| Step"
        << " | Time"
        << " | || u - u_b ||_L2(Omega)"
        << " | || grad( u - u_b ) ||_L2(Omega)"
        << " | || u - u_b ||_H1(Omega)"
        << " | || p - p_b ||_L2(Omega)"
        << " | || nu^(+1/2) grad( u - u_b ) ||_L2(Omega)"
        << " | || nu^(-1/2) ( p - p_b ) ||_L2(Omega)"
        << " | || sigma^(+1/2) ( u - u_h ) ||_L2(Omega)"
        << " | || Phi^(+1/2) (p - p_h) ||_L2(Omega)"
        << " | || u - u_e ||_L2(Omega)"
        << " | || grad( u - u_e ) ||_L2(Omega)"
        << " | || u - u_e ||_H1(Omega)"
        << " | || p - p_e ||_L2(Omega)"
        << " | || sigma^(+1/2) ( u - u_h ) ||_L2(Omega)"
        << " | || Phi^(+1/2) (p - p_h) ||_L2(Omega)"
        << " | || nu^(+1/2) grad( u - u_e ) ||_L2(Omega)"
        << " | || nu^(-1/2) ( p - p_e) ||_L2(Omega)"
        << " | || nu^(+1/2) (u_b - u_e) ||_H1/2(Gamma)"
        << " | || nu^(+1/2) grad( u_b - u_e )*n ||_H-1/2(Gamma)"
        << " | || nu^(-1/2) (p_b - p_e)*n |_H-1/2(Gamma)"
        << " | || (u*n)_inflow (u_b - u_e) ||_L2(Gamma)"
        << " | || (sigma*h+|u|+nu/h)^(+1/2) (u_b - u_e)*n ||_L2(Gamma)"
        << " |  | sin(x) ( u,x - u,x exact ) | (background)"
        << " |  | sin(x) ( u,x - u,x exact ) | (embedded)"
        << " |\n";
      f << step_ << " " << time_ << " " << dom_bg_err_vel_L2 << " " << dom_bg_err_vel_H1_semi << " "
        << dom_bg_err_vel_H1 << " " << dom_bg_err_pre_L2 << " " << dom_bg_err_vel_H1_semi_nu_scaled
        << " " << dom_bg_err_pre_L2_nu_scaled << " " << dom_bg_err_vel_L2_sigma_scaled << " "
        << dom_bg_err_pre_L2_Phi_scaled << " " << dom_emb_err_vel_L2 << " "
        << dom_emb_err_vel_H1_semi << " " << dom_emb_err_vel_H1 << " " << dom_emb_err_pre_L2 << " "
        << dom_emb_err_vel_H1_semi_nu_scaled << " " << dom_emb_err_pre_L2_nu_scaled << " "
        << dom_emb_err_vel_L2_sigma_scaled << " " << dom_emb_err_pre_L2_Phi_scaled << " "
        << interf_err_Honehalf << " " << interf_err_Hmonehalf_u << " " << interf_err_Hmonehalf_p
        << " " << interf_err_inflow << " " << interf_err_mass_cons << " " << functional_bg << " "
        << functional_emb << " "
        << "\n";
      f.flush();
      f.close();
    }
    std::ostringstream temp;
    const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
    const std::string fname = simulation + "_time.xfem_abserror";

    if (step_ == 1)
    {
      std::ofstream f;
      f.open(fname.c_str());
      f << "#| " << simulation << "\n";
      f << "#| Step"
        << " | Time"
        << " | || u - u_b ||_L2(Omega)"
        << " | || grad( u - u_b ) ||_L2(Omega)"
        << " | || u - u_b ||_H1(Omega)"
        << " | || p - p_b ||_L2(Omega)"
        << " | || nu^(+1/2) grad( u - u_b ) ||_L2(Omega)"
        << " | || nu^(-1/2) ( p - p_b ) ||_L2(Omega)"
        << " | || u - u_e ||_L2(Omega)"
        << " | || grad( u - u_e ) ||_L2(Omega)"
        << " | || u - u_e ||_H1(Omega)"
        << " | || p - p_e ||_L2(Omega)"
        << " | || nu^(+1/2) grad( u - u_e ) ||_L2(Omega)"
        << " | || nu^(-1/2) ( p - p_e) ||_L2(Omega)"
        << " | || nu^(+1/2) (u_b - u_e) ||_H1/2(Gamma)"
        << " | || nu^(+1/2) grad( u_b - u_e )*n ||_H-1/2(Gamma)"
        << " | || nu^(-1/2) (p_b - p_e)*n |_H-1/2(Gamma)"
        << " | || (u*n)_inflow (u_b - u_e) ||_L2(Gamma)"
        << " | || (sigma*h+|u|+nu/h)^(+1/2) (u_b - u_e)*n ||_L2(Gamma)"
        << " |  | sin(x) ( u,x - u,x exact ) | (background)"
        << " |  | sin(x) ( u,x - u,x exact ) | (embedded)"
        << " |\n";
      f << step_ << " " << time_ << " " << dom_bg_err_vel_L2 << " " << dom_bg_err_vel_H1_semi << " "
        << dom_bg_err_vel_H1 << " " << dom_bg_err_pre_L2 << " " << dom_bg_err_vel_H1_semi_nu_scaled
        << " " << dom_bg_err_pre_L2_nu_scaled << " " << dom_bg_err_vel_L2_sigma_scaled << " "
        << dom_bg_err_pre_L2_Phi_scaled << " " << dom_emb_err_vel_L2 << " "
        << dom_emb_err_vel_H1_semi << " " << dom_emb_err_vel_H1 << " " << dom_emb_err_pre_L2 << " "
        << dom_emb_err_vel_H1_semi_nu_scaled << " " << dom_emb_err_pre_L2_nu_scaled << " "
        << dom_emb_err_vel_L2_sigma_scaled << " " << dom_emb_err_pre_L2_Phi_scaled << " "
        << interf_err_Honehalf << " " << interf_err_Hmonehalf_u << " " << interf_err_Hmonehalf_p
        << " " << interf_err_inflow << " " << interf_err_mass_cons << " " << functional_bg << " "
        << functional_emb << " "
        << "\n";

      f.flush();
      f.close();
    }
    else
    {
      std::ofstream f;
      f.open(fname.c_str(), std::fstream::ate | std::fstream::app);
      f << step_ << " " << time_ << " " << dom_bg_err_vel_L2 << " " << dom_bg_err_vel_H1_semi << " "
        << dom_bg_err_vel_H1 << " " << dom_bg_err_pre_L2 << " " << dom_bg_err_vel_H1_semi_nu_scaled
        << " " << dom_bg_err_pre_L2_nu_scaled << " " << dom_bg_err_vel_L2_sigma_scaled << " "
        << dom_bg_err_pre_L2_Phi_scaled << " " << dom_emb_err_vel_L2 << " "
        << dom_emb_err_vel_H1_semi << " " << dom_emb_err_vel_H1 << " " << dom_emb_err_pre_L2 << " "
        << dom_emb_err_vel_H1_semi_nu_scaled << " " << dom_emb_err_pre_L2_nu_scaled << " "
        << dom_emb_err_vel_L2_sigma_scaled << " " << dom_emb_err_pre_L2_Phi_scaled << " "
        << interf_err_Honehalf << " " << interf_err_Hmonehalf_u << " " << interf_err_Hmonehalf_p
        << " " << interf_err_inflow << " " << interf_err_mass_cons << " " << functional_bg << " "
        << functional_emb << " "
        << "\n";

      f.flush();
      f.close();
    }
  }  // myrank = 0

  return Teuchos::null;
}
