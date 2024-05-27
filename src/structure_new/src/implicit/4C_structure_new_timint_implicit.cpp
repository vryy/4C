/*-----------------------------------------------------------*/
/*! \file

\brief Implicit structural time integration strategy.


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_timint_implicit.hpp"

#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_pstream.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_print.hpp"
#include "4C_solver_nonlin_nox_group.hpp"
#include "4C_solver_nonlin_nox_linearsystem.hpp"
#include "4C_structure_new_impl_generic.hpp"
#include "4C_structure_new_nln_solver_factory.hpp"
#include "4C_structure_new_nln_solver_generic.hpp"
#include "4C_structure_new_predict_factory.hpp"
#include "4C_structure_new_predict_generic.hpp"
#include "4C_structure_new_timint_noxinterface.hpp"
#include "4C_structure_new_utils.hpp"

#include <NOX_Abstract_Group.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::Setup()
{
  // safety check
  check_init();

  STR::TIMINT::Base::Setup();

  // ---------------------------------------------------------------------------
  // cast the base class integrator
  // ---------------------------------------------------------------------------
  implint_ptr_ = Teuchos::rcp_dynamic_cast<STR::IMPLICIT::Generic>(integrator_ptr());

  // ---------------------------------------------------------------------------
  // build NOX interface
  // ---------------------------------------------------------------------------
  Teuchos::RCP<STR::TIMINT::NoxInterface> noxinterface_ptr =
      Teuchos::rcp(new STR::TIMINT::NoxInterface);
  noxinterface_ptr->Init(
      data_global_state_ptr(), implint_ptr_, DBCPtr(), Teuchos::rcp(this, false));
  noxinterface_ptr->Setup();

  // ---------------------------------------------------------------------------
  // build predictor
  // ---------------------------------------------------------------------------
  const enum INPAR::STR::PredEnum& predtype = DataSDyn().GetPredictorType();
  predictor_ptr_ = STR::PREDICT::BuildPredictor(predtype);
  predictor_ptr_->Init(predtype, implint_ptr_, DBCPtr(), data_global_state_ptr(), data_io_ptr(),
      DataSDyn().GetNoxParamsPtr());
  predictor_ptr_->Setup();

  // ---------------------------------------------------------------------------
  // build non-linear solver
  // ---------------------------------------------------------------------------
  const enum INPAR::STR::NonlinSolTech& nlnSolverType = DataSDyn().GetNlnSolverType();
  if (nlnSolverType == INPAR::STR::soltech_singlestep)
    std::cout << "WARNING!!! You are trying to solve implicitly using the \"singlestep\" nonlinear "
                 "solver. This is not encouraged, since it only works for linear statics analysis. "
                 "Please use NLNSOL as \"fullnewton\" for reliable result."
              << std::endl;
  nlnsolver_ptr_ = STR::NLN::SOLVER::BuildNlnSolver(nlnSolverType);
  nlnsolver_ptr_->Init(data_global_state_ptr(), data_s_dyn_ptr(), noxinterface_ptr, implint_ptr_,
      Teuchos::rcp(this, false));
  nlnsolver_ptr_->Setup();

  // set setup flag
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::set_state(const Teuchos::RCP<Epetra_Vector>& x)
{
  integrator_ptr()->set_state(*x);
  ::NOX::Epetra::Vector x_nox(x, ::NOX::Epetra::Vector::CreateView);
  nln_solver().SolutionGroup().setX(x_nox);
  set_state_in_sync_with_nox_group(true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::prepare_partition_step()
{
  check_init_setup();
  FOUR_C_THROW("Not yet implemented!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::prepare_time_step()
{
  check_init_setup();
  // update end time \f$t_{n+1}\f$ of this time step to cope with time step size adaptivity
  /* ToDo Check if this is still necessary. I moved this part to the Update(const double endtime)
   * routine, such it becomes consistent with non-adaptive update routine! See the
   * UpdateStepTime() routine for more information.                             hiermeier 12/2015
   *
  double& time_np = data_global_state().GetTimeNp();
  time_np = data_global_state().GetTimeN() + (*data_global_state().GetDeltaTime())[0]; */

  ::NOX::Abstract::Group& grp = nln_solver().SolutionGroup();
  predictor().Predict(grp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Implicit::Integrate()
{
  check_init_setup();
  FOUR_C_THROW(
      "The function is unused since the ADAPTER::StructureTimeLoop "
      "wrapper gives you all the flexibility you need.");
  return 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Implicit::IntegrateStep()
{
  check_init_setup();
  // do the predictor step
  ::NOX::Abstract::Group& grp = nln_solver().SolutionGroup();
  predictor().Predict(grp);
  return Solve();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::STR::ConvergenceStatus STR::TIMINT::Implicit::Solve()
{
  check_init_setup();
  throw_if_state_not_in_sync_with_nox_group();
  // reset the non-linear solver
  nln_solver().Reset();
  // solve the non-linear problem
  INPAR::STR::ConvergenceStatus convstatus = nln_solver().Solve();
  // return convergence status
  return PerformErrorAction(convstatus);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::update_state_incrementally(Teuchos::RCP<const Epetra_Vector> disiterinc)
{
  if (disiterinc == Teuchos::null) return;

  check_init_setup();
  throw_if_state_not_in_sync_with_nox_group();
  ::NOX::Abstract::Group& grp = nln_solver().SolutionGroup();

  auto* grp_ptr = dynamic_cast<NOX::NLN::Group*>(&grp);
  FOUR_C_ASSERT(grp_ptr != nullptr, "Dynamic cast failed!");

  // cast away const-qualifier for building the Nox Vector
  Teuchos::RCP<Epetra_Vector> mutable_disiterinc =
      Teuchos::rcp(const_cast<Epetra_Vector*>(disiterinc.get()), false);

  // wrap the displacement vector in a nox_epetra_Vector
  Teuchos::RCP<const ::NOX::Epetra::Vector> nox_disiterinc_ptr = Teuchos::rcp(
      new ::NOX::Epetra::Vector(mutable_disiterinc, ::NOX::Epetra::Vector::CreateView));

  // updated the state vector in the nox group
  grp_ptr->computeX(*grp_ptr, *nox_disiterinc_ptr, 1.0);

  // Reset the state variables
  const auto& x_eptra = dynamic_cast<const ::NOX::Epetra::Vector&>(grp_ptr->getX());
  // set the consistent state in the models (e.g. structure and contact models)
  impl_int().ResetModelStates(x_eptra.getEpetraVector());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::determine_stress_strain() { impl_int().determine_stress_strain(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc)
{
  update_state_incrementally(disiterinc);

  Evaluate();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::Evaluate()
{
  check_init_setup();
  throw_if_state_not_in_sync_with_nox_group();
  ::NOX::Abstract::Group& grp = nln_solver().SolutionGroup();

  auto* grp_ptr = dynamic_cast<NOX::NLN::Group*>(&grp);
  FOUR_C_ASSERT(grp_ptr != nullptr, "Dynamic cast failed!");

  // you definitely have to evaluate here. You might be called from a coupled
  // problem and the group might not be aware, that a different state than
  // the internally stored displacements may have changed.
  // This is a hack to get NOX to set IsValid to false.
  grp_ptr->setX(grp_ptr->getX());

  // compute the rhs vector and the stiffness matrix
  grp_ptr->computeFandJacobian();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const ::NOX::Abstract::Group& STR::TIMINT::Implicit::get_solution_group() const
{
  check_init_setup();
  return nln_solver().get_solution_group();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Abstract::Group> STR::TIMINT::Implicit::solution_group_ptr()
{
  check_init_setup();
  return Teuchos::rcpFromRef(nln_solver().SolutionGroup());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::STR::ConvergenceStatus STR::TIMINT::Implicit::PerformErrorAction(
    INPAR::STR::ConvergenceStatus nonlinsoldiv)
{
  check_init_setup();

  if (nonlinsoldiv == INPAR::STR::conv_success)
  {
    // Only relevant, if the input parameter DIVERCONT is used and set to divcontype_ == adapt_step:
    // In this case, the time step size is halved as consequence of a non-converging nonlinear
    // solver. After a prescribed number of converged time steps, the time step is doubled again.
    // The following methods checks, if the time step size can be increased again.
    check_for_time_step_increase(nonlinsoldiv);
    return INPAR::STR::conv_success;
  }
  // get ID of actual processor in parallel
  const int& myrank = data_global_state().GetMyRank();

  // what to do when nonlinear solver does not converge
  switch (GetDivergenceAction())
  {
    case INPAR::STR::divcont_stop:
    {
      // write restart output of last converged step before stopping
      Output(true);

      // we should not get here, FOUR_C_THROW for safety
      FOUR_C_THROW("Nonlinear solver did not converge! ");
      return INPAR::STR::conv_nonlin_fail;
      break;
    }
    case INPAR::STR::divcont_continue:
    {
      if (myrank == 0)
      {
        IO::cout << "\n WARNING: You are continuing your simulation although the nonlinear solver\n"
                    " did not converge in the current time step.\n"
                 << IO::endl;
      }
      return INPAR::STR::conv_success;
      break;
    }
    case INPAR::STR::divcont_repeat_step:
    {
      if (myrank == 0)
        IO::cout << "Nonlinear solver failed to converge repeat time step" << IO::endl;

      // reset step (e.g. quantities on element level or model specific stuff)
      reset_step();

      return INPAR::STR::conv_fail_repeat;
      break;
    }
    case INPAR::STR::divcont_halve_step:
    {
      if (myrank == 0)
      {
        IO::cout << "Nonlinear solver failed to converge at time t= " << GetTimeNp()
                 << ". Divide timestep in half. "
                 << "Old time step: " << GetDeltaTime() << IO::endl
                 << "New time step: " << 0.5 * GetDeltaTime() << IO::endl
                 << IO::endl;
      }

      // halve the time step size
      SetDeltaTime(GetDeltaTime() * 0.5);
      // update the number of max time steps if it does not exceed the largest possible value for
      // the type int
      if ((GetStepEnd() - GetStepNp() + 1) > std::numeric_limits<int>::max() - GetStepEnd())
        FOUR_C_THROW(" Your updated step number exceeds largest possible value for type int");
      int endstep = GetStepEnd() + (GetStepEnd() - GetStepNp()) + 1;
      SetStepEnd(endstep);
      // reset timen_ because it is set in the constructor
      SetTimeNp(GetTimeN() + GetDeltaTime());
      // reset step (e.g. quantities on element level or model specific stuff)
      reset_step();

      integrator().update_constant_state_contributions();

      return INPAR::STR::conv_fail_repeat;
      break;
    }
    case INPAR::STR::divcont_adapt_step:
    {
      if (myrank == 0)
      {
        IO::cout << "Nonlinear solver failed to converge at time t= " << GetTimeNp()
                 << ". Divide timestep in half. "
                 << "Old time step: " << GetDeltaTime() << IO::endl
                 << "New time step: " << 0.5 * GetDeltaTime() << IO::endl
                 << IO::endl;
      }

      // halve the time step size
      SetDeltaTime(GetDeltaTime() * 0.5);
      // update the number of max time steps if it does not exceed the largest possible value for
      // the type int
      if ((GetStepEnd() - GetStepNp() + 1) > std::numeric_limits<int>::max() - GetStepEnd())
        FOUR_C_THROW(" Your updated step number exceeds largest possible value for type int");
      int endstep = GetStepEnd() + (GetStepEnd() - GetStepNp()) + 1;
      SetStepEnd(endstep);
      // reset timen_ because it is set in the constructor
      SetTimeNp(GetTimeN() + GetDeltaTime());

      set_div_con_refine_level(get_div_con_refine_level() + 1);
      set_div_con_num_fine_step(0);

      if (get_div_con_refine_level() == get_max_div_con_refine_level())
        FOUR_C_THROW(
            "Maximal divercont refinement level reached. Adapt your time basic time step size!");

      // reset step (e.g. quantities on element level or model specific stuff)
      reset_step();

      integrator().update_constant_state_contributions();

      return INPAR::STR::conv_fail_repeat;
      break;
    }
    case INPAR::STR::divcont_rand_adapt_step:
    case INPAR::STR::divcont_rand_adapt_step_ele_err:
    {
      // generate random number between 0.51 and 1.99 (as mean value of random
      // numbers generated on all processors), alternating values larger
      // and smaller than 1.0
      double proc_randnum_get = ((double)rand() / (double)RAND_MAX);
      double proc_randnum = proc_randnum_get;
      double randnum = 1.0;
      const Epetra_Comm& comm = discretization()->Comm();
      comm.SumAll(&proc_randnum, &randnum, 1);
      const double numproc = comm.NumProc();
      randnum /= numproc;
      if (get_random_time_step_factor() > 1.0)
        set_random_time_step_factor(randnum * 0.49 + 0.51);
      else if (get_random_time_step_factor() < 1.0)
        set_random_time_step_factor(randnum * 0.99 + 1.0);
      else
        set_random_time_step_factor(randnum * 1.48 + 0.51);

      if (myrank == 0)
      {
        IO::cout << "Nonlinear solver failed to converge: modifying time-step size by random "
                    "number between 0.51 and 1.99 -> here: "
                 << get_random_time_step_factor() << " !" << IO::endl;
      }
      // multiply time-step size by random number
      SetDeltaTime(GetDeltaTime() * get_random_time_step_factor());
      // update maximum number of time steps
      int endstep = (1.0 / get_random_time_step_factor()) * GetStepEnd() +
                    (1.0 - (1.0 / get_random_time_step_factor())) * GetStepNp() + 1;
      if (endstep > std::numeric_limits<int>::max())
        FOUR_C_THROW(" Your updated step number exceeds largest possible value for type int");
      SetStepEnd(endstep);
      // reset timen_ because it is set in the constructor
      SetTimeNp(GetTimeN() + GetDeltaTime());
      // reset step (e.g. quantities on element level or model specific stuff)
      reset_step();

      integrator().update_constant_state_contributions();

      return INPAR::STR::conv_fail_repeat;
      break;
    }
    case INPAR::STR::divcont_adapt_penaltycontact:
    {
      // adapt penalty and search parameter
      FOUR_C_THROW("Not yet implemented for new structure time integration");
      break;
    }
    case INPAR::STR::divcont_repeat_simulation:
    {
      if (nonlinsoldiv == INPAR::STR::conv_nonlin_fail and myrank == 0)
      {
        IO::cout << "Nonlinear solver failed to converge and DIVERCONT = "
                    "repeat_simulation, hence leaving structural time integration "
                 << IO::endl;
      }
      else if (nonlinsoldiv == INPAR::STR::conv_lin_fail and myrank == 0)
      {
        IO::cout << "Linear solver failed to converge and DIVERCONT = "
                    "repeat_simulation, hence leaving structural time integration "
                 << IO::endl;
      }
      else if (nonlinsoldiv == INPAR::STR::conv_ele_fail and myrank == 0)
      {
        IO::cout << "Element failure in form of a negative Jacobian determinant and DIVERCONT = "
                    "repeat_simulation, hence leaving structural time integration "
                 << IO::endl;
      }
      return nonlinsoldiv;  // so that time loop will be aborted
      break;
    }
    default:
      FOUR_C_THROW("Unknown DIVER_CONT case");
      return INPAR::STR::conv_nonlin_fail;
      break;
  }
  return INPAR::STR::conv_success;  // make compiler happy
}  // PerformErrorAction()

/*-----------------------------------------------------------------------------*
 * check, if according to divercont flag                             meier 01/15
 * time step size can be increased
 *-----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::check_for_time_step_increase(INPAR::STR::ConvergenceStatus& status)
{
  check_init_setup();

  const int maxnumfinestep = 4;

  if (GetDivergenceAction() != INPAR::STR::divcont_adapt_step)
    return;
  else if (status == INPAR::STR::conv_success and get_div_con_refine_level() != 0)
  {
    set_div_con_num_fine_step(get_div_con_num_fine_step() + 1);

    if (get_div_con_num_fine_step() == maxnumfinestep)
    {
      // increase the step size if the remaining number of steps is a even number
      if (((GetStepEnd() - GetStepNp()) % 2) == 0 and GetStepEnd() != GetStepNp())
      {
        if (data_global_state().GetMyRank() == 0)
          IO::cout << "Nonlinear solver successful. Double timestep size!" << IO::endl;

        set_div_con_refine_level(get_div_con_refine_level() - 1);
        set_div_con_num_fine_step(0);

        SetStepEnd(GetStepEnd() - (GetStepEnd() - GetStepNp()) / 2);

        // double the time step size
        SetDeltaTime(GetDeltaTime() * 2.0);
      }
      else  // otherwise we have to wait one more time step until the step size can be increased
      {
        set_div_con_num_fine_step(get_div_con_num_fine_step() - 1);
      }
    }
    return;
  }
}  // check_for_time_step_increase()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::print_jacobian_in_matlab_format(const NOX::NLN::Group& curr_grp) const
{
  typedef CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>
      LinalgBlockSparseMatrix;

  if (not GetDataIO().is_write_jacobian_to_matlab()) return;

  // create file name
  std::stringstream filebase;

  filebase << "str_jacobian"
           << "_step-" << GetDataGlobalState().GetStepNp() << "_nlniter-"
           << nlnsolver_ptr_->GetNumNlnIterations();

  std::stringstream filename;
  filename << GetDataIO().GetOutputPtr()->Output()->FileName() << "_" << filebase.str() << ".mtl";

  if (GetDataGlobalState().GetMyRank() == 0)
    std::cout << "Writing structural jacobian to \"" << filename.str() << "\"\n";

  Teuchos::RCP<const ::NOX::Epetra::LinearSystem> linear_system = curr_grp.getLinearSystem();

  Teuchos::RCP<const NOX::NLN::LinearSystem> nln_lin_system =
      Teuchos::rcp_dynamic_cast<const NOX::NLN::LinearSystem>(linear_system, true);

  const enum NOX::NLN::LinSystem::OperatorType jac_type =
      nln_lin_system->get_jacobian_operator_type();

  Teuchos::RCP<const Epetra_Operator> jac_ptr = nln_lin_system->getJacobianOperator();

  switch (jac_type)
  {
    case NOX::NLN::LinSystem::LinalgSparseMatrix:
    {
      Teuchos::RCP<const CORE::LINALG::SparseMatrix> sparse_matrix =
          Teuchos::rcp_dynamic_cast<const CORE::LINALG::SparseMatrix>(jac_ptr, true);
      CORE::LINALG::PrintMatrixInMatlabFormat(
          filename.str().c_str(), *sparse_matrix->EpetraMatrix());

      break;
    }
    case NOX::NLN::LinSystem::LinalgBlockSparseMatrix:
    {
      Teuchos::RCP<const LinalgBlockSparseMatrix> block_matrix =
          Teuchos::rcp_dynamic_cast<const LinalgBlockSparseMatrix>(jac_ptr, true);
      CORE::LINALG::PrintBlockMatrixInMatlabFormat(filename.str(), *block_matrix);

      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported NOX::NLN::LinSystem::OperatorType: \"%s\"",
          NOX::NLN::LinSystem::OperatorType2String(jac_type).c_str());
      exit(EXIT_FAILURE);
    }
  }

  // print sparsity pattern to file
  //  CORE::LINALG::PrintSparsityToPostscript( *(SystemMatrix()->EpetraMatrix()) );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::compute_condition_number(const NOX::NLN::Group& grp) const
{
  std::string name_prefix;

  double max_rev = -1.0;
  double min_rev = -1.0;

  const INPAR::STR::ConditionNumber cond_type = GetDataIO().ConditionNumberType();
  switch (cond_type)
  {
    case INPAR::STR::ConditionNumber::gmres_estimate:
    {
      const_cast<NOX::NLN::Group&>(grp).computeJacobianConditionNumber(1000, 1.0e-3, 1000, true);

      name_prefix = "gmres_estimated_";

      break;
    }
    case INPAR::STR::ConditionNumber::max_min_ev_ratio:
    case INPAR::STR::ConditionNumber::one_norm:
    case INPAR::STR::ConditionNumber::inf_norm:
    {
      const NOX::NLN::LinSystem::ConditionNumber nox_cond_type =
          STR::NLN::Convert2NoxConditionNumberType(cond_type);
      const std::string nox_cond_type_str(
          NOX::NLN::LinSystem::ConditionNumber2String(nox_cond_type));

      const_cast<NOX::NLN::Group&>(grp).compute_serial_jacobian_condition_number(
          nox_cond_type, true);

      name_prefix = nox_cond_type_str + "_";

      if (cond_type == INPAR::STR::ConditionNumber::max_min_ev_ratio)
      {
        max_rev = grp.get_jacobian_max_real_eigenvalue();
        min_rev = grp.get_jacobian_min_real_eigenvalue();
      }

      break;
    }
    case INPAR::STR::ConditionNumber::none:
      return;
    default:
      FOUR_C_THROW("Unknown ConditionNumber type!");
      exit(EXIT_FAILURE);
  }

  const double cond_num = grp.getJacobianConditionNumber();

  // create file name
  std::stringstream filebase;

  const int iter = nlnsolver_ptr_->GetNumNlnIterations();
  filebase << "str_" << name_prefix << "condition_number_step-" << GetDataGlobalState().GetStepNp();

  std::stringstream filename;
  filename << GetDataIO().GetOutputPtr()->Output()->FileName() << "_" << filebase.str() << ".data";

  if (GetDataGlobalState().GetMyRank() == 0)
  {
    static bool first_execution = true;

    std::cout << "Writing structural condition number to \"" << filename.str() << "\"\n";

    if (first_execution)
    {
      // clear the file and write the header
      std::ofstream of(filename.str(), std::ios_base::out);
      of << std::setw(24) << "iter" << std::setw(24) << "cond" << std::setw(24) << "max_rev"
         << std::setw(24) << "min_rev\n";
      of.close();
      first_execution = false;
    }

    std::ofstream of(filename.str(), std::ios_base::out | std::ios_base::app);
    of << std::setw(24) << iter << std::setprecision(16);
    of << std::setw(24) << std::scientific << cond_num << std::setw(24) << std::scientific
       << max_rev << std::setw(24) << std::scientific << min_rev << "\n";
    of.close();
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::DynamicType STR::TIMINT::Implicit::MethodName() const
{
  return implint_ptr_->MethodName();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Implicit::MethodSteps() const { return implint_ptr_->MethodSteps(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Implicit::method_order_of_accuracy_dis() const
{
  return implint_ptr_->method_order_of_accuracy_dis();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Implicit::method_order_of_accuracy_vel() const
{
  return implint_ptr_->method_order_of_accuracy_vel();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TIMINT::Implicit::method_lin_err_coeff_dis() const
{
  return implint_ptr_->method_lin_err_coeff_dis();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TIMINT::Implicit::method_lin_err_coeff_vel() const
{
  return implint_ptr_->method_lin_err_coeff_vel();
}

FOUR_C_NAMESPACE_CLOSE
