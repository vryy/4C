/*-----------------------------------------------------------*/
/*! \file

\brief Implicit structural time integration strategy.

\maintainer Anh-Tu Vuong

\level 3

*/
/*-----------------------------------------------------------*/

#include "../drt_io/io_pstream.H"
#include "../drt_lib/drt_discret.H"

#include "str_timint_implicit.H"
#include "str_impl_generic.H"
#include "str_predict_generic.H"
#include "str_nln_solver_generic.H"
#include "str_timint_noxinterface.H"
#include "str_utils.H"

#include "../drt_io/io_gmsh.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include <NOX_Abstract_Group.H>

#include "../solver_nonlin_nox/nox_nln_linearsystem.H"
#include "../solver_nonlin_nox/nox_nln_group.H"

#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_utils_densematrix_print.H"

// factories
#include "str_predict_factory.H"
#include "str_nln_solver_factory.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::Implicit::Implicit()
    : implint_ptr_(Teuchos::null), nlnsolver_ptr_(Teuchos::null), grp_ptr_(Teuchos::null)
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::Setup()
{
  // safety check
  CheckInit();

  STR::TIMINT::Base::Setup();

  // ---------------------------------------------------------------------------
  // cast the base class integrator
  // ---------------------------------------------------------------------------
  implint_ptr_ = Teuchos::rcp_dynamic_cast<STR::IMPLICIT::Generic>(IntegratorPtr());

  // ---------------------------------------------------------------------------
  // build NOX interface
  // ---------------------------------------------------------------------------
  Teuchos::RCP<STR::TIMINT::NoxInterface> noxinterface_ptr =
      Teuchos::rcp(new STR::TIMINT::NoxInterface());
  noxinterface_ptr->Init(DataGlobalStatePtr(), implint_ptr_, DBCPtr(), Teuchos::rcp(this, false));
  noxinterface_ptr->Setup();

  // ---------------------------------------------------------------------------
  // build predictor
  // ---------------------------------------------------------------------------
  const enum INPAR::STR::PredEnum& predtype = DataSDyn().GetPredictorType();
  predictor_ptr_ = STR::PREDICT::BuildPredictor(predtype);
  predictor_ptr_->Init(predtype, implint_ptr_, DBCPtr(), DataGlobalStatePtr(), DataIOPtr(),
      DataSDyn().GetMutableNoxParamsPtr());
  predictor_ptr_->Setup();

  // ---------------------------------------------------------------------------
  // build non-linear solver
  // ---------------------------------------------------------------------------
  const enum INPAR::STR::NonlinSolTech& nlnSolverType = DataSDyn().GetNlnSolverType();
  nlnsolver_ptr_ = STR::NLN::SOLVER::BuildNlnSolver(nlnSolverType);
  nlnsolver_ptr_->Init(DataGlobalStatePtr(), DataSDynPtr(), noxinterface_ptr, implint_ptr_,
      Teuchos::rcp(this, false));
  nlnsolver_ptr_->Setup();

  // set setup flag
  issetup_ = true;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::SetState(const Teuchos::RCP<Epetra_Vector>& x)
{
  IntegratorPtr()->SetState(*x);
  NOX::Epetra::Vector x_nox(x, NOX::Epetra::Vector::CreateView);
  NlnSolver().SolutionGroup().setX(x_nox);
  SetStateInSyncWithNOXGroup(true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::PreparePartitionStep()
{
  CheckInitSetup();
  dserror("Not yet implemented!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::PrepareTimeStep()
{
  CheckInitSetup();
  // update end time \f$t_{n+1}\f$ of this time step to cope with time step size adaptivity
  /* ToDo Check if this is still necessary. I moved this part to the Update(const double endtime)
   * routine, such it becomes consistent with non-adaptive update routine! See the
   * UpdateStepTime() routine for more information.                             hiermeier 12/2015
   *
  double& time_np = DataGlobalState().GetMutableTimeNp();
  time_np = DataGlobalState().GetTimeN() + (*DataGlobalState().GetDeltaTime())[0]; */

  NOX::Abstract::Group& grp = NlnSolver().SolutionGroup();
  Predictor().Predict(grp);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Implicit::Integrate()
{
  CheckInitSetup();
  dserror(
      "The function is unused since the ADAPTER::StructureTimeLoop "
      "wrapper gives you all the flexibility you need.");
  return 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Implicit::IntegrateStep()
{
  CheckInitSetup();
  // do the predictor step
  NOX::Abstract::Group& grp = NlnSolver().SolutionGroup();
  Predictor().Predict(grp);
  return Solve();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::STR::ConvergenceStatus STR::TIMINT::Implicit::Solve()
{
  CheckInitSetup();
  ThrowIfStateNotInSyncWithNOXGroup();
  // reset the non-linear solver
  NlnSolver().Reset();
  // solve the non-linear problem
  INPAR::STR::ConvergenceStatus convstatus = NlnSolver().Solve();
  // return convergence status
  return PerformErrorAction(convstatus);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::UpdateStateIncrementally(Teuchos::RCP<const Epetra_Vector> disiterinc)
{
  if (disiterinc == Teuchos::null) return;

  CheckInitSetup();
  ThrowIfStateNotInSyncWithNOXGroup();
  NOX::Abstract::Group& grp = NlnSolver().SolutionGroup();

  NOX::NLN::Group* grp_ptr = dynamic_cast<NOX::NLN::Group*>(&grp);
  if (grp_ptr == NULL) dserror("Dynamic cast failed!");

  // cast away const-qualifier for building the Nox Vector
  Teuchos::RCP<Epetra_Vector> mutable_disiterinc =
      Teuchos::rcp(const_cast<Epetra_Vector*>(disiterinc.get()), false);

  // wrap the displacement vector in a nox_epetra_Vector
  Teuchos::RCP<const NOX::Epetra::Vector> nox_disiterinc_ptr =
      Teuchos::rcp(new NOX::Epetra::Vector(mutable_disiterinc, NOX::Epetra::Vector::CreateView));

  // updated the state vector in the nox group
  grp_ptr->computeX(*grp_ptr, *nox_disiterinc_ptr, 1.0);

  // Reset the state variables
  const NOX::Epetra::Vector& x_eptra = dynamic_cast<const NOX::Epetra::Vector&>(grp_ptr->getX());
  // set the consistent state in the models (e.g. structure and contact models)
  ImplInt().ResetModelStates(x_eptra.getEpetraVector());

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc)
{
  UpdateStateIncrementally(disiterinc);

  Evaluate();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::Evaluate()
{
  CheckInitSetup();
  ThrowIfStateNotInSyncWithNOXGroup();
  NOX::Abstract::Group& grp = NlnSolver().SolutionGroup();

  NOX::NLN::Group* grp_ptr = dynamic_cast<NOX::NLN::Group*>(&grp);
  if (grp_ptr == NULL) dserror("Dynamic cast failed!");

  // you definitely have to evaluate here. You might be called from a coupled
  // problem and the group might not be aware, that a different state than
  // the internally stored displacements may have changed.
  // This is a hack to get NOX to set IsValid to false.
  grp_ptr->setX(grp_ptr->getX());

  // compute the rhs vector and the stiffness matrix
  grp_ptr->computeFandJacobian();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const NOX::Abstract::Group& STR::TIMINT::Implicit::GetSolutionGroup() const
{
  CheckInitSetup();
  return NlnSolver().GetSolutionGroup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Abstract::Group> STR::TIMINT::Implicit::SolutionGroupPtr()
{
  CheckInitSetup();
  return Teuchos::rcpFromRef(NlnSolver().SolutionGroup());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::STR::ConvergenceStatus STR::TIMINT::Implicit::PerformErrorAction(
    INPAR::STR::ConvergenceStatus nonlinsoldiv)
{
  CheckInitSetup();

  if (nonlinsoldiv == INPAR::STR::conv_success)
  {
    // Only relevant, if the input parameter DIVERCONT is used and set to divcontype_ == adapt_step:
    // In this case, the time step size is halved as consequence of a non-converging nonlinear
    // solver. After a prescribed number of converged time steps, the time step is doubled again.
    // The following methods checks, if the time step size can be increased again.
    CheckForTimeStepIncrease(nonlinsoldiv);
    return INPAR::STR::conv_success;
  }
  // get ID of actual processor in parallel
  const int& myrank = DataGlobalState().GetMyRank();

  // what to do when nonlinear solver does not converge
  switch (GetDivergenceAction())
  {
    case INPAR::STR::divcont_stop:
    {
      // write restart output of last converged step before stopping
      Output(true);

      // we should not get here, dserror for safety
      dserror("Nonlinear solver did not converge! ");
      return INPAR::STR::conv_nonlin_fail;
      break;
    }
    case INPAR::STR::divcont_continue:
    {
      if (myrank == 0)
        IO::cout << "\n WARNING: You are continuing your simulation although the nonlinear solver\n"
                    " did not converge in the current time step.\n"
                 << IO::endl;
      return INPAR::STR::conv_success;
      break;
    }
    case INPAR::STR::divcont_repeat_step:
    {
      if (myrank == 0)
        IO::cout << "Nonlinear solver failed to converge repeat time step" << IO::endl;

      // reset step (e.g. quantities on element level or model specific stuff)
      ResetStep();

      return INPAR::STR::conv_fail_repeat;
      break;
    }
    case INPAR::STR::divcont_halve_step:
    {
      if (myrank == 0)
        IO::cout << "Nonlinear solver failed to converge at time t= " << GetTimeNp()
                 << ". Divide timestep in half. "
                 << "Old time step: " << GetDeltaTime() << IO::endl
                 << "New time step: " << 0.5 * GetDeltaTime() << IO::endl
                 << IO::endl;

      // halve the time step size
      SetDeltaTime(GetDeltaTime() * 0.5);
      // update the number of max time steps
      int endstep = GetStepEnd() + (GetStepEnd() - GetStepNp()) + 1;
      if (endstep > std::numeric_limits<int>::max())
        dserror(" Your updated step number exceeds largest possible value for type int");
      SetStepEnd(endstep);
      // reset timen_ because it is set in the constructor
      SetTimeNp(GetTimeN() + GetDeltaTime());
      // reset step (e.g. quantities on element level or model specific stuff)
      ResetStep();

      Integrator().UpdateConstantStateContributions();

      return INPAR::STR::conv_fail_repeat;
      break;
    }
    case INPAR::STR::divcont_adapt_step:
    {
      if (myrank == 0)
        IO::cout << "Nonlinear solver failed to converge at time t= " << GetTimeNp()
                 << ". Divide timestep in half. "
                 << "Old time step: " << GetDeltaTime() << IO::endl
                 << "New time step: " << 0.5 * GetDeltaTime() << IO::endl
                 << IO::endl;

      // halve the time step size
      SetDeltaTime(GetDeltaTime() * 0.5);
      // update the number of max time steps
      int endstep = GetStepEnd() + (GetStepEnd() - GetStepNp()) + 1;
      if (endstep > std::numeric_limits<int>::max())
        dserror(" Your updated step number exceeds largest possible value for type int");
      SetStepEnd(endstep);
      // reset timen_ because it is set in the constructor
      SetTimeNp(GetTimeN() + GetDeltaTime());

      SetDivConRefineLevel(GetDivConRefineLevel() + 1);
      SetDivConNumFineStep(0);

      if (GetDivConRefineLevel() == GetMaxDivConRefineLevel())
        dserror(
            "Maximal divercont refinement level reached. Adapt your time basic time step size!");

      // reset step (e.g. quantities on element level or model specific stuff)
      ResetStep();

      Integrator().UpdateConstantStateContributions();

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
      const Epetra_Comm& comm = DiscretizationInterface()->Comm();
      comm.SumAll(&proc_randnum, &randnum, 1);
      const double numproc = comm.NumProc();
      randnum /= numproc;
      if (GetRandomTimeStepFactor() > 1.0)
        SetRandomTimeStepFactor(randnum * 0.49 + 0.51);
      else if (GetRandomTimeStepFactor() < 1.0)
        SetRandomTimeStepFactor(randnum * 0.99 + 1.0);
      else
        SetRandomTimeStepFactor(randnum * 1.48 + 0.51);

      if (myrank == 0)
        IO::cout << "Nonlinear solver failed to converge: modifying time-step size by random "
                    "number between 0.51 and 1.99 -> here: "
                 << GetRandomTimeStepFactor() << " !" << IO::endl;
      // multiply time-step size by random number
      SetDeltaTime(GetDeltaTime() * GetRandomTimeStepFactor());
      // update maximum number of time steps
      int endstep = (1.0 / GetRandomTimeStepFactor()) * GetStepEnd() +
                    (1.0 - (1.0 / GetRandomTimeStepFactor())) * GetStepNp() + 1;
      if (endstep > std::numeric_limits<int>::max())
        dserror(" Your updated step number exceeds largest possible value for type int");
      SetStepEnd(endstep);
      // reset timen_ because it is set in the constructor
      SetTimeNp(GetTimeN() + GetDeltaTime());
      // reset step (e.g. quantities on element level or model specific stuff)
      ResetStep();

      Integrator().UpdateConstantStateContributions();

      return INPAR::STR::conv_fail_repeat;
      break;
    }
    case INPAR::STR::divcont_adapt_penaltycontact:
    {
      // adapt penalty and search parameter
      dserror("Not yet implemented for new structure time integration");
      break;
    }
    case INPAR::STR::divcont_repeat_simulation:
    {
      if (nonlinsoldiv == INPAR::STR::conv_nonlin_fail and myrank == 0)
        IO::cout << "Nonlinear solver failed to converge and DIVERCONT = "
                    "repeat_simulation, hence leaving structural time integration "
                 << IO::endl;
      else if (nonlinsoldiv == INPAR::STR::conv_lin_fail and myrank == 0)
        IO::cout << "Linear solver failed to converge and DIVERCONT = "
                    "repeat_simulation, hence leaving structural time integration "
                 << IO::endl;
      else if (nonlinsoldiv == INPAR::STR::conv_ele_fail and myrank == 0)
        IO::cout << "Element failure in form of a negative Jacobian determinant and DIVERCONT = "
                    "repeat_simulation, hence leaving structural time integration "
                 << IO::endl;
      return nonlinsoldiv;  // so that time loop will be aborted
      break;
    }
    default:
      dserror("Unknown DIVER_CONT case");
      return INPAR::STR::conv_nonlin_fail;
      break;
  }
  return INPAR::STR::conv_success;  // make compiler happy
}  // PerformErrorAction()

/*-----------------------------------------------------------------------------*
 * check, if according to divercont flag                             meier 01/15
 * time step size can be increased
 *-----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::CheckForTimeStepIncrease(INPAR::STR::ConvergenceStatus& status)
{
  CheckInitSetup();

  const int maxnumfinestep = 4;

  if (GetDivergenceAction() != INPAR::STR::divcont_adapt_step)
    return;
  else if (status == INPAR::STR::conv_success and GetDivConRefineLevel() != 0)
  {
    SetDivConNumFineStep(GetDivConNumFineStep() + 1);

    if (GetDivConNumFineStep() == maxnumfinestep)
    {
      // increase the step size if the remaining number of steps is a even number
      if (((GetStepEnd() - GetStepNp()) % 2) == 0 and GetStepEnd() != GetStepNp())
      {
        if (DataGlobalState().GetMyRank() == 0)
          IO::cout << "Nonlinear solver successful. Double timestep size!" << IO::endl;

        SetDivConRefineLevel(GetDivConRefineLevel() - 1);
        SetDivConNumFineStep(0);

        SetStepEnd(GetStepEnd() - (GetStepEnd() - GetStepNp()) / 2);

        // double the time step size
        SetDeltaTime(GetDeltaTime() * 2.0);
      }
      else  // otherwise we have to wait one more time step until the step size can be increased
      {
        SetDivConNumFineStep(GetDivConNumFineStep() - 1);
      }
    }
    return;
  }
}  // CheckForTimeStepIncrease()

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::PrintJacobianInMatlabFormat(const NOX::NLN::Group& curr_grp) const
{
  typedef LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> LinalgBlockSparseMatrix;

  if (not GetDataIO().IsWriteJacobianToMatlab()) return;

  // create file name
  std::stringstream filebase;

  filebase << "str_jacobian"
           << "_step-" << GetDataGlobalState().GetStepNp() << "_nlniter-"
           << nlnsolver_ptr_->GetNumNlnIterations();

  std::stringstream filename;
  filename << GetDataIO().GetOutputPtr()->Output()->FileName() << "_" << filebase.str() << ".mtl";

  if (GetDataGlobalState().GetMyRank() == 0)
    std::cout << "Writing structural jacobian to \"" << filename.str() << "\"\n";

  Teuchos::RCP<const NOX::Epetra::LinearSystem> linear_system = curr_grp.getLinearSystem();

  Teuchos::RCP<const NOX::NLN::LinearSystem> nln_lin_system =
      Teuchos::rcp_dynamic_cast<const NOX::NLN::LinearSystem>(linear_system, true);

  const enum NOX::NLN::LinSystem::OperatorType jac_type = nln_lin_system->getJacobianOperatorType();

  Teuchos::RCP<const Epetra_Operator> jac_ptr = nln_lin_system->getJacobianOperator();

  switch (jac_type)
  {
    case NOX::NLN::LinSystem::LinalgSparseMatrix:
    {
      Teuchos::RCP<const LINALG::SparseMatrix> sparse_matrix =
          Teuchos::rcp_dynamic_cast<const LINALG::SparseMatrix>(jac_ptr, true);
      LINALG::PrintMatrixInMatlabFormat(filename.str().c_str(), *sparse_matrix->EpetraMatrix());

      break;
    }
    case NOX::NLN::LinSystem::LinalgBlockSparseMatrix:
    {
      Teuchos::RCP<const LinalgBlockSparseMatrix> block_matrix =
          Teuchos::rcp_dynamic_cast<const LinalgBlockSparseMatrix>(jac_ptr, true);
      LINALG::PrintBlockMatrixInMatlabFormat(filename.str(), *block_matrix);

      break;
    }
    default:
    {
      dserror("Unsupported NOX::NLN::LinSystem::OperatorType: \"%s\"",
          NOX::NLN::LinSystem::OperatorType2String(jac_type).c_str());
      exit(EXIT_FAILURE);
    }
  }

  // print sparsity pattern to file
  //  LINALG::PrintSparsityToPostscript( *(SystemMatrix()->EpetraMatrix()) );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::ComputeConditionNumber(const NOX::NLN::Group& grp) const
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

      const_cast<NOX::NLN::Group&>(grp).computeSerialJacobianConditionNumber(nox_cond_type, true);

      name_prefix = nox_cond_type_str + "_";

      if (cond_type == INPAR::STR::ConditionNumber::max_min_ev_ratio)
      {
        max_rev = grp.getJacobianMaxRealEigenvalue();
        min_rev = grp.getJacobianMinRealEigenvalue();
      }

      break;
    }
    case INPAR::STR::ConditionNumber::none:
      return;
    default:
      dserror("Unknown ConditionNumber type!");
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
