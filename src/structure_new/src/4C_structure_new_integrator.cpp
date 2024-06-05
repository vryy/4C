/*-----------------------------------------------------------*/
/*! \file

\brief Generic class for all explicit/implicit time integrators.


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_integrator.hpp"

#include "4C_global_data.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_solver_nonlin_nox_aux.hpp"
#include "4C_structure_new_dbc.hpp"
#include "4C_structure_new_model_evaluator.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_model_evaluator_structure.hpp"
#include "4C_structure_new_monitor_dbc.hpp"
#include "4C_structure_new_nox_nln_str_linearsystem.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_timint_noxinterface.hpp"
#include "4C_utils_epetra_exceptions.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_function.hpp"

#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::Integrator::Integrator()
    : isinit_(false),
      issetup_(false),
      modelevaluator_ptr_(Teuchos::null),
      eval_data_ptr_(Teuchos::null),
      sdyn_ptr_(Teuchos::null),
      gstate_ptr_(Teuchos::null),
      io_ptr_(Teuchos::null),
      dbc_ptr_(Teuchos::null),
      timint_ptr_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::Init(const Teuchos::RCP<STR::TIMINT::BaseDataSDyn>& sdyn_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataIO>& io_ptr, const Teuchos::RCP<STR::Dbc>& dbc_ptr,
    const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr)
{
  issetup_ = false;

  sdyn_ptr_ = sdyn_ptr;
  gstate_ptr_ = gstate_ptr;
  io_ptr_ = io_ptr;
  dbc_ptr_ = dbc_ptr;
  timint_ptr_ = timint_ptr;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::Setup()
{
  check_init();
  // ---------------------------------------------------------------------------
  // build model evaluator data container
  // ---------------------------------------------------------------------------
  eval_data_ptr_ = Teuchos::rcp(new STR::MODELEVALUATOR::Data());
  eval_data_ptr_->Init(timint_ptr_);
  eval_data_ptr_->Setup();

  // ---------------------------------------------------------------------------
  // build model evaluator
  // ---------------------------------------------------------------------------
  modelevaluator_ptr_ = Teuchos::rcp(new STR::ModelEvaluator());
  modelevaluator_ptr_->Init(
      eval_data_ptr_, sdyn_ptr_, gstate_ptr_, io_ptr_, Teuchos::rcp(this, false), timint_ptr_);
  modelevaluator_ptr_->Setup();

  // ---------------------------------------------------------------------------
  // build monitor for a tensile test
  // ---------------------------------------------------------------------------
  monitor_dbc_ptr_ = Teuchos::rcp(new STR::MonitorDbc);
  monitor_dbc_ptr_->Init(io_ptr_, *gstate_ptr_->GetDiscret(), *gstate_ptr_, *dbc_ptr_);
  monitor_dbc_ptr_->Setup();

  mt_energy_.Setup();

  // the issetup_ flag is not set here!!!
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::set_initial_displacement(
    const INPAR::STR::InitialDisp init, const int startfuncno)
{
  switch (init)
  {
    case INPAR::STR::initdisp_zero_disp:
    {
      global_state().GetDisN()->PutScalar(0.0);
      global_state().GetDisNp()->PutScalar(0.0);

      break;
    }
    case INPAR::STR::initdisp_disp_by_function:
    {
      const Epetra_Map* dofrowmap = global_state().GetDiscret()->dof_row_map();

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < global_state().GetDiscret()->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        const CORE::Nodes::Node* lnode = global_state().GetDiscret()->lRowNode(lnodeid);

        // the set of degrees of freedom associated with the node
        const std::vector<int> nodedofset = global_state().GetDiscret()->Dof(0, lnode);

        // loop nodal dofs
        for (unsigned int d = 0; d < nodedofset.size(); ++d)
        {
          const int dofgid = nodedofset[d];
          const int doflid = dofrowmap->LID(dofgid);

          // evaluate component k of spatial function
          const double initialval =
              GLOBAL::Problem::Instance()
                  ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(startfuncno - 1)
                  .Evaluate(lnode->X().data(), global_state().GetTimeN(), d);

          const int err = global_state().GetDisN()->ReplaceMyValues(1, &initialval, &doflid);
          if (err != 0) FOUR_C_THROW("dof not on proc");
        }
      }

      // initialize also the solution vector
      global_state().GetDisNp()->Update(1.0, *global_state().GetDisN(), 0.0);

      break;
    }
    default:
      FOUR_C_THROW("Unknown option for initial displacement: %d", init);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::check_init() const { FOUR_C_ASSERT(is_init(), "Call Init() first!"); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::check_init_setup() const
{
  FOUR_C_ASSERT(is_init() and is_setup(), "Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::ResetModelStates(const Epetra_Vector& x)
{
  check_init_setup();
  ModelEval().ResetStates(x);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::equilibrate_initial_state()
{
  check_init();

  // temporary right-hand-side
  Teuchos::RCP<Epetra_Vector> rhs_ptr =
      Teuchos::rcp(new Epetra_Vector(*global_state().DofRowMapView(), true));
  // wrap the rhs_ptr in a nox_epetra_Vector
  Teuchos::RCP<::NOX::Epetra::Vector> nox_rhs_ptr =
      Teuchos::rcp(new ::NOX::Epetra::Vector(rhs_ptr, ::NOX::Epetra::Vector::CreateView));

  // initialize a temporal structural stiffness matrix
  Teuchos::RCP<CORE::LINALG::SparseOperator> stiff_ptr =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*global_state().DofRowMapView(), 81, true, true));


  // overwrite initial state vectors with Dirichlet BCs
  // note that we get accelerations resulting from inhomogeneous Dirichlet conditions here
  const double& timen = (*global_state().GetMultiTime())[0];
  Teuchos::RCP<Epetra_Vector> disnp_ptr = global_state().GetDisNp();
  Teuchos::RCP<Epetra_Vector> velnp_ptr = global_state().GetVelNp();
  Teuchos::RCP<Epetra_Vector> accnp_ptr = global_state().GetAccNp();
  dbc().apply_dirichlet_bc(timen, disnp_ptr, velnp_ptr, accnp_ptr, false);


  // ---------------------------------------------------------------------------
  // compute the mass matrix
  // ---------------------------------------------------------------------------
  // set the evaluate parameters of the current base class
  reset_eval_params();
  // !!! evaluate the initial state !!!
  EvalData().SetTotalTime(gstate_ptr_->GetTimeN());

  // initialize the mass matrix and the Rayleigh damping matrix (optional)
  if (not ModelEval().initialize_inertia_and_damping(*disnp_ptr, *stiff_ptr))
    FOUR_C_THROW("initialize_inertia_and_damping failed!");

  /* If we are restarting the simulation, we do not have to calculate a
   * consistent acceleration, since we get it anyway from the restart file.
   * Furthermore, we keep the update routines untouched. Actually, the only
   * thing which is necessary is the calculation of the mass matrix. So we are
   * done at this point.                                    hiermeier 06/16 */
  /* However, if we need to restart the initial state, e.g. when starting dynamics
   * analysis from the static one, we should re-calculate the initial state
   * for a warm start-up */
  if (timint_ptr_->IsRestarting() && !timint_ptr_->is_restarting_initial_state())
  {
    return;
  }

  // build the entire initial right-hand-side
  if (not ModelEval().ApplyInitialForce(*disnp_ptr, *rhs_ptr))
    FOUR_C_THROW("ApplyInitialForce failed!");

  // add inertial and viscous contributions to rhs
  /* note: this needs to be done 'manually' here because in the RHS evaluation
   * routine of an ordinary time step, these contributions are scaled by weighting
   * factors inside the time integration scheme (e.g. alpha_f/m for GenAlpha) */
  rhs_ptr->Update(1.0, *global_state().GetFinertialNp(), 1.0);
  rhs_ptr->Update(1.0, *global_state().GetFviscoNp(), 1.0);

  /* Meier 2015: Here, we copy the mass matrix in the stiffness block in order to
   * not perform the Dirichlet conditions on the constant mass matrix later on.
   * This is necessary since we need the original mass matrix (without blanked
   * rows) on the Dirichlet DoFs in order to calculate correct reaction
   * forces.*/
  stiff_ptr->Add(*global_state().GetMassMatrix(), false, 1.0, 0.0);
  stiff_ptr->Complete();

  // treatment of elements with special element technology (e.g. pressure DOFs)
  global_state().apply_element_technology_to_acceleration_system(*stiff_ptr, *rhs_ptr);

  // ---------------------------------------------------------------------------
  // build a NOX::NLN::STR::LinearSystem
  // ---------------------------------------------------------------------------
  // get the structural linear solver
  std::map<enum NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>> str_linsolver;
  str_linsolver[NOX::NLN::sol_structure] =
      tim_int().GetDataSDyn().GetLinSolvers().at(INPAR::STR::model_structure);

  // copy the nox parameter-list
  Teuchos::ParameterList p_nox = tim_int().GetDataSDyn().GetNoxParams();
  NOX::NLN::AUX::set_printing_parameters(p_nox, global_state().GetComm());

  // create a copy of the initial displacement vector
  Teuchos::RCP<Epetra_Vector> soln_ptr =
      Teuchos::rcp(new Epetra_Vector(*global_state().DofRowMapView(), true));
  // wrap the soln_ptr in a nox_epetra_Vector
  Teuchos::RCP<::NOX::Epetra::Vector> nox_soln_ptr =
      Teuchos::rcp(new ::NOX::Epetra::Vector(soln_ptr, ::NOX::Epetra::Vector::CreateView));

  // Check if we are using a Newton direction
  std::string dir_str = p_nox.sublist("Direction").get<std::string>("Method");
  if (dir_str == "User Defined")
    dir_str = p_nox.sublist("Direction").get<std::string>("User Defined Method");
  if (dir_str != "Newton" and dir_str != "Modified Newton")
    FOUR_C_THROW(
        "The EquilibriateState predictor is currently only working for the "
        "direction-method \"Newton\".");

  // create the linear system
  // printing parameters
  Teuchos::ParameterList& p_print = p_nox.sublist("Printing", true);
  // linear solver parameters
  Teuchos::ParameterList& p_ls =
      p_nox.sublist("Direction", true).sublist("Newton", true).sublist("Linear Solver", true);

  Teuchos::RCP<NOX::NLN::LinearSystem> linsys_ptr = Teuchos::rcp(new NOX::NLN::STR::LinearSystem(
      p_print, p_ls, str_linsolver, Teuchos::null, Teuchos::null, stiff_ptr, *nox_soln_ptr));

  // (re)set the linear solver parameters
  p_ls.set<int>("Number of Nonlinear Iterations", 0);
  p_ls.set<int>("Current Time Step", global_state().GetStepNp());
  // ToDo Get the actual tolerance value
  p_ls.set<double>("Wanted Tolerance", 1.0e-6);
  // ---------------------------------------------------------------------------
  /* Meier 2015: Due to the Dirichlet conditions applied to the mass matrix, we
   * solely solve for the accelerations at non-Dirichlet DoFs while the
   * resulting accelerations at the Dirichlet-DoFs will be zero.
   * accelerations at DoFs with inhomogeneous Dirichlet conditions were already
   * added above. */
  // ---------------------------------------------------------------------------
  // solve the linear system
  if (stiff_ptr->NormInf() == 0.0) FOUR_C_THROW("You are about to invert a singular matrix!");

  linsys_ptr->applyJacobianInverse(p_ls, *nox_rhs_ptr, *nox_soln_ptr);
  nox_soln_ptr->scale(-1.0);

  // get the solution vector and add it into the acceleration vector
  accnp_ptr->Update(1.0, nox_soln_ptr->getEpetraVector(), 1.0);

  // re-build the entire initial right-hand-side with correct accelerations
  ModelEval().ApplyInitialForce(*disnp_ptr, *rhs_ptr);

  // call update routines to copy states from t_{n+1} to t_{n}
  // note that the time step is not incremented
  PreUpdate();
  UpdateStepState();
  UpdateStepElement();
  post_update();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::Integrator::current_state_is_equilibrium(const double& tol)
{
  check_init();

  // temporary right-hand-side
  Teuchos::RCP<Epetra_Vector> rhs_ptr =
      Teuchos::rcp(new Epetra_Vector(*global_state().DofRowMapView(), true));

  // overwrite initial state vectors with Dirichlet BCs
  const double& timen = (*global_state().GetMultiTime())[0];
  Teuchos::RCP<Epetra_Vector> disnp_ptr = global_state().GetDisNp();
  Teuchos::RCP<Epetra_Vector> velnp_ptr = global_state().GetVelNp();
  Teuchos::RCP<Epetra_Vector> accnp_ptr = global_state().GetAccNp();
  dbc().apply_dirichlet_bc(timen, disnp_ptr, velnp_ptr, accnp_ptr, false);

  // set the evaluate parameters of the current base class
  reset_eval_params();
  // !!! evaluate the initial state !!!
  EvalData().SetTotalTime(gstate_ptr_->GetTimeN());

  // build the entire right-hand-side
  ModelEval().ApplyInitialForce(*disnp_ptr, *rhs_ptr);

  // add viscous contributions to rhs
  rhs_ptr->Update(1.0, *global_state().GetFviscoNp(), 1.0);
  // add inertial contributions to rhs
  rhs_ptr->Update(1.0, *global_state().GetFinertialNp(), 1.0);

  double resnorm = 0.0;
  rhs_ptr->NormInf(&resnorm);

  return (resnorm < tol ? true : false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::determine_stress_strain()
{
  check_init_setup();
  ModelEval().determine_stress_strain();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::DetermineEnergy()
{
  check_init_setup();
  ModelEval().DetermineEnergy();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::Integrator::GetModelValue(const Epetra_Vector& x)
{
  FOUR_C_THROW(
      "This routine is not supported in the currently active time "
      "integration scheme.");
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::Integrator::get_total_mid_time_str_energy(const Epetra_Vector& x)
{
  check_init_setup();
  if (not mt_energy_.is_correctly_configured())
    FOUR_C_THROW(
        "You are trying to compute the mid-time energy in case of a non-static"
        " simulation, but you have not specified the desired energy averaging type."
        " Please add a meaningful MIDTIME_ENERGY_TYPE to the ---STRUCTURAL DYNAMIC"
        " section of your Input file.");

  Teuchos::RCP<const Epetra_Vector> disnp_ptr = global_state().ExtractDisplEntries(x);
  const Epetra_Vector& disnp = *disnp_ptr;

  set_state(disnp);

  Teuchos::RCP<const Epetra_Vector> velnp_ptr = global_state().GetVelNp();
  const Epetra_Vector& velnp = *velnp_ptr;

  EvalData().clear_values_for_all_energy_types();
  STR::MODELEVALUATOR::Structure& str_model =
      dynamic_cast<STR::MODELEVALUATOR::Structure&>(Evaluator(INPAR::STR::model_structure));

  Teuchos::RCP<const Epetra_Vector> dis_avg =
      mt_energy_.Average(disnp, *global_state().GetDisN(), GetIntParam());
  Teuchos::RCP<const Epetra_Vector> vel_avg =
      mt_energy_.Average(velnp, *global_state().GetVelN(), GetIntParam());

  str_model.DetermineEnergy(*dis_avg, vel_avg.get(), true);
  mt_energy_.int_energy_np_ = EvalData().GetEnergyData(STR::internal_energy);
  mt_energy_.kin_energy_np_ = EvalData().GetEnergyData(STR::kinetic_energy);
  global_state().GetFextNp()->Dot(*dis_avg, &mt_energy_.ext_energy_np_);

  CORE::IO::cout(CORE::IO::debug) << __LINE__ << " -- " << __PRETTY_FUNCTION__ << "\n";
  mt_energy_.Print(CORE::IO::cout.os(CORE::IO::debug));

  return mt_energy_.GetTotal();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::update_structural_energy()
{
  if (not mt_energy_.StoreEnergyN()) return;

  get_total_mid_time_str_energy(*global_state().GetDisNp());
  mt_energy_.CopyNpToN();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::determine_optional_quantity()
{
  check_init_setup();
  ModelEval().determine_optional_quantity();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::Integrator::determine_element_volumes(
    const Epetra_Vector& x, Teuchos::RCP<Epetra_Vector>& ele_vols)
{
  check_init_setup();
  STR::MODELEVALUATOR::Generic& model = Evaluator(INPAR::STR::model_structure);
  STR::MODELEVALUATOR::Structure& smodel = dynamic_cast<STR::MODELEVALUATOR::Structure&>(model);

  return smodel.determine_element_volumes(x, ele_vols);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::OutputStepState(CORE::IO::DiscretizationWriter& iowriter) const
{
  check_init_setup();
  ModelEval().OutputStepState(iowriter);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::MonitorDbc(CORE::IO::DiscretizationWriter& writer) const
{
  monitor_dbc_ptr_->Execute(writer);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::runtime_pre_output_step_state()
{
  check_init_setup();
  ModelEval().runtime_pre_output_step_state();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::runtime_output_step_state() const
{
  check_init_setup();
  ModelEval().runtime_output_step_state();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::PostOutput()
{
  check_init_setup();
  ModelEval().PostOutput();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::ResetStepState()
{
  check_init_setup();
  ModelEval().ResetStepState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::Integrator::get_condensed_update_norm(
    const enum NOX::NLN::StatusTest::QuantityType& qtype) const
{
  check_init_setup();

  double myupdatenorm = eval_data_ptr_->GetMyUpdateNorm(qtype);
  const enum ::NOX::Abstract::Vector::NormType normtype =
      eval_data_ptr_->get_update_norm_type(qtype);

  return get_condensed_global_norm(qtype, normtype, myupdatenorm);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::Integrator::get_condensed_previous_sol_norm(
    const enum NOX::NLN::StatusTest::QuantityType& qtype) const
{
  check_init_setup();

  double myprevsolnorm = eval_data_ptr_->get_my_previous_sol_norm(qtype);
  const enum ::NOX::Abstract::Vector::NormType normtype =
      eval_data_ptr_->get_update_norm_type(qtype);

  return get_condensed_global_norm(qtype, normtype, myprevsolnorm);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::Integrator::get_condensed_solution_update_rms(
    const enum NOX::NLN::StatusTest::QuantityType& qtype) const
{
  check_init_setup();
  // global relative mean square norm
  double grmsnorm = 0.0;
  // get proc data
  double myrmsnorm = eval_data_ptr_->GetMyRMSNorm(qtype);
  // get total dof number
  int gdofnumber = get_condensed_dof_number(qtype);
  // sum over all procs
  gstate_ptr_->GetComm().SumAll(&myrmsnorm, &grmsnorm, 1);

  // calculate the root mean square and return the value
  return (std::sqrt(grmsnorm) / static_cast<double>(gdofnumber));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::Integrator::get_condensed_dof_number(
    const enum NOX::NLN::StatusTest::QuantityType& qtype) const
{
  check_init_setup();
  // global dof number of the given quantity
  int gdofnumber = 0.0;
  int mydofnumber = eval_data_ptr_->GetMyDofNumber(qtype);
  gstate_ptr_->GetComm().SumAll(&mydofnumber, &gdofnumber, 1);
  return gdofnumber;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::Integrator::get_condensed_global_norm(
    const enum NOX::NLN::StatusTest::QuantityType& qtype,
    const enum ::NOX::Abstract::Vector::NormType& normtype, double& mynorm) const
{
  double gnorm = 0;

  switch (normtype)
  {
    case ::NOX::Abstract::Vector::OneNorm:
    {
      gstate_ptr_->GetComm().SumAll(&mynorm, &gnorm, 1);
      break;
    }
    case ::NOX::Abstract::Vector::TwoNorm:
    {
      gstate_ptr_->GetComm().SumAll(&mynorm, &gnorm, 1);
      gnorm = std::sqrt(gnorm);
      break;
    }
    case ::NOX::Abstract::Vector::MaxNorm:
    {
      gstate_ptr_->GetComm().MaxAll(&mynorm, &gnorm, 1);
      break;
    }
  }
  return gnorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::ModelEvaluator& STR::Integrator::ModelEval()
{
  check_init();
  return *modelevaluator_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::ModelEvaluator& STR::Integrator::ModelEval() const
{
  check_init();
  return *modelevaluator_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const STR::ModelEvaluator> STR::Integrator::ModelEvalPtr() const
{
  check_init();
  return modelevaluator_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Generic& STR::Integrator::Evaluator(const INPAR::STR::ModelType& mt)
{
  check_init_setup();
  return ModelEval().Evaluator(mt);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::MODELEVALUATOR::Generic& STR::Integrator::Evaluator(
    const INPAR::STR::ModelType& mt) const
{
  check_init_setup();
  return ModelEval().Evaluator(mt);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::MODELEVALUATOR::Data& STR::Integrator::EvalData() const
{
  check_init();
  return *eval_data_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Data& STR::Integrator::EvalData()
{
  check_init();
  return *eval_data_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataSDyn& STR::Integrator::s_dyn()
{
  check_init();
  return *sdyn_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::BaseDataSDyn& STR::Integrator::s_dyn() const
{
  check_init();
  return *sdyn_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::BaseDataGlobalState& STR::Integrator::global_state() const
{
  check_init();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataGlobalState& STR::Integrator::global_state()
{
  check_init();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::Dbc& STR::Integrator::dbc()
{
  check_init();
  return *dbc_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::Dbc& STR::Integrator::GetDbc() const
{
  check_init();
  return *dbc_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::Base& STR::Integrator::tim_int() const
{
  check_init();
  return *timint_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::CreateBackupState(const Epetra_Vector& dir)
{
  check_init_setup();
  ModelEval().CreateBackupState(dir);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::recover_from_backup_state()
{
  check_init_setup();
  ModelEval().recover_from_backup_state();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::Integrator::MidTimeEnergy::MidTimeEnergy(const Integrator& integrator)
    : integrator_(integrator), avg_type_(INPAR::STR::midavg_vague)
{
  /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::MidTimeEnergy::Print(std::ostream& os) const
{
  const double time_fac = integrator_.GetIntParam();

  os << "--- Contributions of time step n+1 (current)\n";
  os << "strain energy                         = " << int_energy_np_ << "\n";
  os << "external energy (dead load/potential) = " << ext_energy_np_ << "\n";
  os << "kinetic energy                        = " << kin_energy_np_ << "\n";

  if (avg_type_ == INPAR::STR::midavg_trlike)
  {
    os << "--- Contributions of time step n   (previously accepted)\n";
    os << "strain energy                         = " << int_energy_n_ << "\n";
    os << "external energy (dead load/potential) = " << ext_energy_n_ << "\n";
    os << "kinetic energy                        = " << kin_energy_n_ << "\n";
    os << std::string(40, '-') << "\n";
    os << "TimInt factor at t_n+1 (current)      = " << 1.0 - time_fac << "\n";
    os << "TimInt factor at t_n   (previous)     = " << time_fac << "\n";
  }
  os << std::string(40, '=') << "\n";
  os << "Total structural energy               = " << GetTotal() << "\n";
  os << std::string(40, '=') << "\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::Integrator::MidTimeEnergy::GetTotal() const
{
  const double fac_n = integrator_.GetIntParam();
  const double fac_np = 1.0 - fac_n;

  double total_energy = 0.0;
  total_energy = int_energy_np_ - kin_energy_np_ - ext_energy_np_;
  if (avg_type_ == INPAR::STR::midavg_trlike)
  {
    const double energy_n = int_energy_n_ - kin_energy_n_ - ext_energy_n_;
    total_energy = fac_np * total_energy + fac_n * energy_n;
  }

  return total_energy;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::Integrator::MidTimeEnergy::Average(
    const Epetra_Vector& state_np, const Epetra_Vector& state_n, const double fac_n) const
{
  const double fac_np = 1.0 - fac_n;

  Teuchos::RCP<Epetra_Vector> state_avg = Teuchos::rcp(new Epetra_Vector(state_np));
  switch (avg_type_)
  {
    case INPAR::STR::midavg_vague:
    case INPAR::STR::midavg_trlike:
      return state_avg;
    case INPAR::STR::midavg_imrlike:
    {
      state_avg->Update(fac_n, state_n, fac_np);
      return state_avg;
    }
    default:
      FOUR_C_THROW("Don't know what to do for the given MidAvg type.");
      exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::MidTimeEnergy::CopyNpToN()
{
  kin_energy_n_ = kin_energy_np_;
  int_energy_n_ = int_energy_np_;
  ext_energy_n_ = ext_energy_np_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::Integrator::MidTimeEnergy::is_correctly_configured() const
{
  FOUR_C_ASSERT(issetup_, "Call Setup() first.");

  if (avg_type_ == INPAR::STR::midavg_vague)
  {
    if (integrator_.s_dyn().GetDynamicType() != INPAR::STR::dyna_statics) return false;
  }
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::Integrator::MidTimeEnergy::StoreEnergyN() const
{
  return avg_type_ == INPAR::STR::midavg_trlike;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::MidTimeEnergy::Setup()
{
  avg_type_ = integrator_.s_dyn().get_mid_time_energy_type();
  issetup_ = true;
}

FOUR_C_NAMESPACE_CLOSE
