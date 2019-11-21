/*----------------------------------------------------------------------------*/
/** \file

\brief Merge the mutual functionality of the XCONTACT state creators.

\maintainer Matthias Mayr

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "xfem_state_creator.H"
#include "xfield_state.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_discret_xfem.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_cut/cut_cutwizard.H"
#include "../drt_structure_xstructure/xstr_xstructure_structure_state.H"

#include "../drt_structure_new/str_timint_factory.H"

#include "../drt_xfem/xfem_dofset.H"
#include "../drt_xfem/xfem_condition_manager.H"



/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XFEM::StateCreator::StateCreator()
    : isinit_(false),
      issetup_(false),
      condition_manager_(Teuchos::null),
      nodal_dofset_strategy_(INPAR::CUT::NDS_Strategy_full),
      volume_cell_gauss_point_by_(),
      bound_cell_gauss_point_by_(),
      gmsh_cut_out_(false),
      numdof_(-1),
      maxnumdofsets_(-1),
      minnumdofsets_(-1),
      include_inner_(false)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::StateCreator::Init(const Teuchos::RCP<XFEM::ConditionManager>& condition_manager,
    const Teuchos::ParameterList& params_xfem, int numdof, int maxnumdofsets, int minnumdofsets,
    bool include_inner)
{
  issetup_ = false;

  condition_manager_ = condition_manager;

  nodal_dofset_strategy_ = DRT::INPUT::IntegralValue<INPAR::CUT::NodalDofSetStrategy>(
      params_xfem, "NODAL_DOFSET_STRATEGY");
  volume_cell_gauss_point_by_ =
      DRT::INPUT::IntegralValue<INPAR::CUT::VCellGaussPts>(params_xfem, "VOLUME_GAUSS_POINTS_BY");
  bound_cell_gauss_point_by_ =
      DRT::INPUT::IntegralValue<INPAR::CUT::BCellGaussPts>(params_xfem, "BOUNDARY_GAUSS_POINTS_BY");

  gmsh_cut_out_ = DRT::INPUT::IntegralValue<int>(params_xfem, "GMSH_CUT_OUT");

  numdof_ = numdof;
  maxnumdofsets_ = maxnumdofsets;
  minnumdofsets_ = minnumdofsets;

  include_inner_ = include_inner;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::StateCreator::Setup()
{
  CheckInit();

  /* Should be filled in one of the derived classes if necessary! This is
   * only a dummy. */

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<XFEM::XFieldState> XFEM::StateCreator::Create(
    const Teuchos::RCP<DRT::DiscretizationXFEM>& xfielddiscret,
    Teuchos::RCP<const Epetra_Vector> back_disp_col, Teuchos::ParameterList& solver_params,
    int step, double time, bool dosetup)
{
  CheckInitSetup();

  //----------------------------------------------------------------------
  // create the XFieldState object
  Teuchos::RCP<XFEM::XFieldState> state = CreateXFieldState();

  //----------------------------------------------------------------------
  // finish the XFieldState object (initialize and setup)
  CompleteState(
      state, xfielddiscret, Teuchos::null, back_disp_col, solver_params, step, time, dosetup);

  return state;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<XFEM::XFieldState> XFEM::StateCreator::Create(
    const Teuchos::RCP<DRT::DiscretizationInterface>& xfielddiscret,
    const Teuchos::RCP<DRT::DiscretizationInterface>& fielddiscret,
    Teuchos::RCP<const Epetra_Vector> back_disp_col, Teuchos::ParameterList& solver_params,
    int step, double time, bool dosetup)
{
  // try to cast things
  Teuchos::RCP<DRT::DiscretizationXFEM> xfielddiscret_ptr =
      Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(xfielddiscret, true);
  Teuchos::RCP<DRT::Discretization> fielddiscret_ptr =
      Teuchos::rcp_dynamic_cast<DRT::Discretization>(fielddiscret, true);
  // if success, do the actual function call
  return Create(
      xfielddiscret_ptr, fielddiscret_ptr, back_disp_col, solver_params, step, time, dosetup);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<XFEM::XFieldState> XFEM::StateCreator::Create(
    const Teuchos::RCP<DRT::DiscretizationXFEM>& xfielddiscret,
    const Teuchos::RCP<DRT::Discretization>& fielddiscret,
    Teuchos::RCP<const Epetra_Vector> back_disp_col, Teuchos::ParameterList& solver_params,
    int step, double time, bool dosetup)
{
  CheckInitSetup();

  //----------------------------------------------------------------------
  // create the XFieldFieldState object
  Teuchos::RCP<XFEM::XFieldState> state = CreateXFieldFieldState();

  //----------------------------------------------------------------------
  // finish the XFieldFieldState object (initialize and setup)
  CompleteState(
      state, xfielddiscret, fielddiscret, back_disp_col, solver_params, step, time, dosetup);

  return state;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::StateCreator::CompleteState(const Teuchos::RCP<XFEM::XFieldState>& state,
    const Teuchos::RCP<DRT::DiscretizationInterface>& xfielddiscret,
    const Teuchos::RCP<DRT::DiscretizationInterface>& fielddiscret,
    Teuchos::RCP<const Epetra_Vector> back_disp_col, Teuchos::ParameterList& solver_params,
    int step, double time, bool dosetup)
{
  // try to cast things
  Teuchos::RCP<DRT::DiscretizationXFEM> xfielddiscret_ptr =
      Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(xfielddiscret, true);

  Teuchos::RCP<DRT::Discretization> fielddiscret_ptr =
      Teuchos::rcp_dynamic_cast<DRT::Discretization>(fielddiscret, true);

  // if success, do the actual function call
  CompleteState(state, xfielddiscret_ptr, fielddiscret_ptr, back_disp_col, solver_params, step,
      time, dosetup);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::StateCreator::CompleteState(const Teuchos::RCP<XFEM::XFieldState>& state,
    const Teuchos::RCP<DRT::DiscretizationXFEM>& xfielddiscret,
    const Teuchos::RCP<DRT::Discretization>& fielddiscret,
    Teuchos::RCP<const Epetra_Vector> back_disp_col, Teuchos::ParameterList& solver_params,
    int step, double time, bool dosetup)
{
  if (state.is_null()) dserror("The state must exist at this point!");

  //----------------------------------------------------------------------
  // create new cut wizard & dof-set
  Teuchos::RCP<GEO::CutWizard> wizard = Teuchos::null;
  Teuchos::RCP<XFEM::XFEMDofSet> xdofset = Teuchos::null;

  CreateNewCutState(xdofset, wizard, xfielddiscret, back_disp_col, solver_params, step);

  //----------------------------------------------------------------------
  // Initialize the state object
  state->Init(condition_manager_, wizard, xdofset, xfielddiscret, fielddiscret);

  // setup the state object
  if (dosetup) state->Setup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum XFEM::StateStatus XFEM::StateCreator::CreateNewCutState(
    Teuchos::RCP<XFEM::XFEMDofSet>& xdofset, Teuchos::RCP<GEO::CutWizard>& wizard,
    const Teuchos::RCP<DRT::DiscretizationXFEM>& xdiscret,
    const Teuchos::RCP<const Epetra_Vector>& back_disp_col, Teuchos::ParameterList& solver_params,
    int step)
{
  CheckInitSetup();

  if (condition_manager_.is_null()) return state_unchanged;

  //----------------------------------------------------------------------
  // create new cut wizard
  CreateCutWizard(xdiscret, back_disp_col, wizard);

  //----------------------------------------------------------------------
  // performs the "CUT"
  wizard->Cut(include_inner_);

  //----------------------------------------------------------------------
  // set the new dofset after cut
  int max_num_my_reserved_dofs_per_node = (maxnumdofsets_)*numdof_;

  // create a new XFEM-dofset
  xdofset =
      Teuchos::rcp(new XFEM::XFEMDofSet(*wizard, max_num_my_reserved_dofs_per_node, *xdiscret));

  const int restart = DRT::Problem::Instance()->Restart();
  if ((step < 1) or restart) minnumdofsets_ = xdiscret->DofRowMap()->MinAllGID();

  // set the minimal GID of xfem dis
  xdofset->SetMinGID(minnumdofsets_);

  return FinishBackgroundDiscretization(xdofset, solver_params, *xdiscret);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::StateCreator::CreateCutWizard(const Teuchos::RCP<DRT::DiscretizationXFEM>& xdiscret,
    const Teuchos::RCP<const Epetra_Vector>& back_disp_col, Teuchos::RCP<GEO::CutWizard>& wizard)
{
  // new wizard using information about cutting sides from the condition_manager
  wizard = Teuchos::rcp(new GEO::CutWizard(xdiscret));

  // Set options for the cut wizard
  wizard->SetOptions(nodal_dofset_strategy_, volume_cell_gauss_point_by_,
      bound_cell_gauss_point_by_, gmsh_cut_out_, true, false, true);

  //----------------------------------------------------------------------
  // set state for all mesh cuttings
  AddCutterStates(wizard);

  //----------------------------------------------------------------------
  // set background state (background mesh displacements and level-set values)
  wizard->SetBackgroundState(back_disp_col, condition_manager_->GetLevelSetFieldCol(),
      condition_manager_->GetLevelSetCouplingGid());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum XFEM::StateStatus XFEM::StateCreator::FinishBackgroundDiscretization(
    const Teuchos::RCP<XFEM::XFEMDofSet>& xdofset, Teuchos::ParameterList& solver_params,
    DRT::DiscretizationXFEM& xdiscret)
{
  // field dofset has nds = 0
  xdiscret.ReplaceDofSet(0, xdofset, true);

  xdiscret.FillComplete(true, false, false);

  // print all dofsets
  xdiscret.GetDofSetProxy()->PrintAllDofsets(xdiscret.Comm());

  //----------------------------------------------------------------------
  // recompute nullspace based on new number of dofs per node
  /* REMARK: this has to be done after replacing the discret' dofset
   * (via discret_->ReplaceDofSet) */
  xdiscret.ComputeNullSpaceIfNecessary(solver_params, true);

  return state_changed;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::StateCreator::AddCutterStates(Teuchos::RCP<GEO::CutWizard>& wizard)
{
  // loop all mesh coupling objects
  for (int mc_idx = 0; mc_idx < condition_manager_->NumMeshCoupling(); ++mc_idx)
  {
    Teuchos::RCP<XFEM::MeshCoupling> mc_coupl = condition_manager_->GetMeshCoupling(mc_idx);

    wizard->AddCutterState(mc_idx, mc_coupl->GetCutterDis(), mc_coupl->GetCutterDispCol(),
        condition_manager_->GetMeshCouplingStartGID(mc_idx));
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<XFEM::XFieldState> XFEM::StateCreator::CreateXFieldState() const
{
  Teuchos::RCP<XFEM::XFieldState> state = Teuchos::null;
  ProblemType prb_type = DRT::Problem::Instance()->GetProblemType();
  switch (prb_type)
  {
    case prb_fluid_xfem:
      dserror("Not yet considered!");
      // state = Teuchos::rcp(new FLD::XFluidState());
      break;
    default:
      dserror("Unsupported problem type ( %d )!", prb_type);
      break;
  }
  return state;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<XFEM::XFieldState> XFEM::StateCreator::CreateXFieldFieldState() const
{
  Teuchos::RCP<XFEM::XFieldState> state = Teuchos::null;
  ProblemType prb_type = DRT::Problem::Instance()->GetProblemType();
  switch (prb_type)
  {
    case prb_fluid_xfem:
      dserror("Not yet considered!");
      //      state = Teuchos::rcp( new FLD::XFluidFluidState() );
      break;
    case prb_xcontact:
      state =
          Teuchos::rcp_dynamic_cast<XFEM::XFieldState>(STR::TIMINT::BuildDataGlobalState(), true);
      break;
    default:
      dserror("Unsupported problem type (%d)!", prb_type);
      break;
  }
  return state;
}
