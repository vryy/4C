/*----------------------------------------------------------------------------*/
/*!
\file ad_str_xcontact.cpp

\brief Interface between the xcontact algorithm and the structure algorithm

\maintainer Michael Hiermeier

\date Jun 17, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/


#include "ad_str_xcontact.H"

#include "../drt_contact_xcontact/xcontact_multi_discretization_wrapper.H"
#include "../drt_contact_xcontact/str_model_evaluator_xcontact.H"

#include "../drt_xfem/xfield_state.H"

#include "../drt_structure_new/str_timint_base.H"

#include "../drt_contact/contact_abstract_strategy.H"

#include "../drt_lib/drt_discret_xfem.H"

#include "../linalg/linalg_solver.H"

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_contact_xcontact/xcontact_cutwizard.H"
#include "../drt_contact_xcontact/xcontact_state_creator.H"

#include "../drt_io/io_pstream.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
ADAPTER::StructureXContact::StructureXContact(const Teuchos::RCP<ADAPTER::Structure>& structure)
    : ADAPTER::StructureWrapper(structure),
      isinit_(false),
      issetup_(false),
      structure_ptr_(Teuchos::rcp_dynamic_cast<StructureNew>(structure, true)),
      multi_discret_ptr_(Teuchos::rcp_dynamic_cast<XCONTACT::MultiDiscretizationWrapper>(
          structure_ptr_->DiscretizationInterface(), true)),
      p_xfem_general_ptr_(Teuchos::null),
      state_creator_ptr_(Teuchos::null),
      state_ptr_(Teuchos::null),
      num_dof_per_node_(-1),
      max_num_dof_sets_(-1),
      state_count_(-1),
      iscontact_(false)
{
  /* intentionally left blank */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureXContact::Init(
    const Teuchos::ParameterList& p_xfem_general, const int& num_dof_per_node)
{
  issetup_ = false;

  p_xfem_general_ptr_ = Teuchos::rcpFromRef<const Teuchos::ParameterList>(p_xfem_general);
  num_dof_per_node_ = num_dof_per_node;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureXContact::Setup()
{
  CheckInit();

  // get the maximum number of dof sets
  max_num_dof_sets_ = p_xfem_general_ptr_->get<int>("MAX_NUM_DOFSETS");

  // ---------------------------------------------------------------------
  // build, initialize and setup the state creator
  // ---------------------------------------------------------------------
  state_creator_ptr_ = Teuchos::rcp(new XCONTACT::StateCreator());
  state_creator_ptr_->Init(
      Teuchos::null, *p_xfem_general_ptr_, num_dof_per_node_, max_num_dof_sets_, 0, true);
  state_creator_ptr_->Setup();

  // ---------------------------------------------------------------------
  // create initial state
  // ---------------------------------------------------------------------
  /* The XFieldState class created here is indexed with 0. All further
   * first cuts of a new time-step have then index 1 and have to be reset
   * to 0 in PrepareTimeStep(). ToDo */
  state_count_ = 0;
  /* The GlobalState is half finished as well as the time integration
   * strategy. So we do not clear it here, but complete it instead. */
  CompleteInitialState();

  /* set the initial contact status to inactive. */
  XContactModel().SetContactStatus(iscontact_);

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureXContact::CompleteInitialState()
{
  state_ptr_ = Structure().XFieldState();

  /* complete the state object (be aware that the condition manager is
   * undefined and we do not create a new cut wizard.) */
  StateCreator().CompleteState(state_ptr_, MultiDiscret().DiscretPtr(XFEM::xstructure),
      MultiDiscret().DiscretPtr(XFEM::structure), Teuchos::null,
      Structure().LinearSolver()->Params(), Structure().GetTimeNp(), Structure().GetStepNp(), true);

  Structure().Setup();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureXContact::RecreateState()
{
  //---------------------------------------------------------------------------
  // get level-set values in nodal row layout
  STR::MODELEVALUATOR::XContact& xcontact = XContactModel();

  Teuchos::RCP<Epetra_Vector> levelset_values =
      Teuchos::rcp(new Epetra_Vector(xcontact.LevelSetValues()));
  levelset_values->ReplaceMap(xcontact.Strategy().SlRowNodes());

  //---------------------------------------------------------------------------
  // create the new state class ( vectors, matrices ... )
  StateCreator().Recreate(state_ptr_, structure_ptr_, MultiDiscret(), Teuchos::null,
      levelset_values, Structure().LinearSolver()->Params(), Structure().GetStepNp(),
      Structure().GetTimeNp(), true);

  //---------------------------------------------------------------------------
  // increase the state counter
  ++state_count_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureXContact::PrepareTimeStep()
{
  CheckInitSetup();
  // predict the new displacement field
  Structure().PrepareTimeStep();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool ADAPTER::StructureXContact::IsComingIntoContact()
{
  CheckInitSetup();

  // We are already in contact. Nothing to do.
  if (iscontact_) return false;

  XContactModel().EvaluateWeightedGap();
  const Epetra_Vector& wgap = XContactModel().GetWeightedGap();

  // check if we are in contact
  double min_gap_value = 0.0;
  wgap.MinValue(&min_gap_value);

  // update the iscontact flag
  iscontact_ = (min_gap_value <= 0.0);

  // return true, if we are coming into contact
  return iscontact_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& ADAPTER::StructureXContact::GetWeightedGap() const
{
  CheckInitSetup();
  return XContactModel().GetWeightedGap();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::STR::ConvergenceStatus ADAPTER::StructureXContact::Solve()
{
  CheckInitSetup();

  PreSolve();

  INPAR::STR::ConvergenceStatus status = INPAR::STR::conv_success;

  // ---------------------------------------------------------------------
  // Solve the direct contact problem
  // ---------------------------------------------------------------------
  if (iscontact_)
  {
    IO::cout << __LINE__ << " -- " << __PRETTY_FUNCTION__
             << ": Solve for the contact state is not yet implemented!" << IO::endl;
  }
  else
    status = Structure().Solve();

  // --- nothing more to do, if the contact state was not solved ---------
  if (not iscontact_) return status;

  // ---------------------------------------------------------------------
  // Solve m-times the sensibility problem
  // ---------------------------------------------------------------------

  // ---------------------------------------------------------------------
  // Solve the shape optimization problem
  // ---------------------------------------------------------------------

  return status;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureXContact::PreSolve()
{
  if (not iscontact_) return;

  RecreateState();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::StructureXContact::SetScaTraValuesInStructure_Np(const Epetra_Vector& phi_np)
{
  CheckInitSetup();
  STR::MODELEVALUATOR::XContact& xcontact = XContactModel();

  if (not iscontact_)
  {
    xcontact.SetLevelSetValuesPtr(Teuchos::null);
    return;
  }

  xcontact.SetLevelSetValuesPtr(MultiDiscret().ScaTra2Contact(phi_np, XFEM::xstructure));

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Comm& ADAPTER::StructureXContact::Comm() const { return MultiDiscret().Comm(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::XContact& ADAPTER::StructureXContact::XContactModel()
{
  try
  {
    STR::MODELEVALUATOR::XContact& evaluator = dynamic_cast<STR::MODELEVALUATOR::XContact&>(
        Structure().ModelEvaluator(INPAR::STR::model_contact));
    return evaluator;
  }
  catch (const std::bad_cast& e)
  {
    dserror(
        "Dynamic cast to STR::MODELEVALUATOR::XContact failed!\\"
        "( throw = \" %s \" )",
        e.what());
  }
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::MODELEVALUATOR::XContact& ADAPTER::StructureXContact::XContactModel() const
{
  try
  {
    const STR::MODELEVALUATOR::XContact& evaluator =
        dynamic_cast<const STR::MODELEVALUATOR::XContact&>(
            Structure().ModelEvaluator(INPAR::STR::model_contact));
    return evaluator;
  }
  catch (const std::bad_cast& e)
  {
    dserror("Dynamic cast to STR::MODELEVALUATOR::XContact failed!");
  }
  exit(EXIT_FAILURE);
}
