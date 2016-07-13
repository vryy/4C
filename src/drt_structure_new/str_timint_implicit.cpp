/*-----------------------------------------------------------*/
/*!
\file str_timint_implicit.cpp

\brief Implicit structural time integration strategy.

\maintainer Michael Hiermeier

\date Aug 13, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_timint_implicit.H"
#include "str_impl_generic.H"
#include "str_predict_generic.H"
#include "str_nln_solver_generic.H"
#include "str_timint_noxinterface.H"

#include "../drt_io/io_gmsh.H"
#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

// factories
#include "str_predict_factory.H"
#include "str_nln_solver_factory.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::Implicit::Implicit()
    : implint_ptr_(Teuchos::null),
      nlnsolver_ptr_(Teuchos::null),
      grp_ptr_(Teuchos::null)
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
  noxinterface_ptr->Init(DataGlobalStatePtr(),
      implint_ptr_,DBCPtr(),Teuchos::rcp(this,false));
  noxinterface_ptr->Setup();

  // ---------------------------------------------------------------------------
  // build predictor
  // ---------------------------------------------------------------------------
  const enum INPAR::STR::PredEnum& predtype =
      DataSDyn().GetPredictorType();
  predictor_ptr_ = STR::PREDICT::BuildPredictor(predtype);
  predictor_ptr_->Init(predtype,implint_ptr_,DBCPtr(),DataGlobalStatePtr(),
      DataIOPtr(),DataSDyn().GetMutableNoxParamsPtr());
  predictor_ptr_->Setup();

  // ---------------------------------------------------------------------------
  // build non-linear solver
  // ---------------------------------------------------------------------------
  const enum INPAR::STR::NonlinSolTech& nlnSolverType =
      DataSDyn().GetNlnSolverType();
  nlnsolver_ptr_ = STR::NLN::SOLVER::BuildNlnSolver(nlnSolverType);
  nlnsolver_ptr_->Init(DataGlobalStatePtr(),
      DataSDynPtr(), noxinterface_ptr,implint_ptr_,Teuchos::rcp(this,false));
  nlnsolver_ptr_->Setup();

  // set setup flag
  issetup_ = true;

  return;
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

  // ToDo prepare contact for new time step
  // PrepareStepContact();
  NOX::Abstract::Group& grp = NlnSolver().SolutionGroup();
  Predictor().Predict(grp);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Implicit::Integrate()
{
  CheckInitSetup();
  dserror("The function is unused since the ADAPTER::StructureTimeLoop "
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
  // reset the non-linear solver
  NlnSolver().Reset();
  // solve the non-linear problem
  INPAR::STR::ConvergenceStatus convstatus = NlnSolver().Solve();
  // return convergence status
  return PerformErrorAction(convstatus);
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
Teuchos::RCP<Epetra_Vector> STR::TIMINT::Implicit::SolveRelaxationLinear()
{
  CheckInitSetup();
  dserror("FixMe: SolveRelaxationLinear() should be implemented here!");
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::STR::ConvergenceStatus STR::TIMINT::Implicit::PerformErrorAction(INPAR::STR::ConvergenceStatus nonlinsoldiv)
{
  CheckInitSetup();

  if (nonlinsoldiv==INPAR::STR::conv_success)
  {
    CheckForTimeStepIncrease(nonlinsoldiv);
    return INPAR::STR::conv_success;
  }

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
      // we should not get here, dserror for safety
      dserror("Nonlinear solver did not converge! ");
      return INPAR::STR::conv_nonlin_fail;
      break;
    }
    case INPAR::STR::divcont_repeat_step:
    {
      IO::cout << "Nonlinear solver failed to converge repeat time step"
               << IO::endl;

      // reset step (e.g. quantities on element level)
      ResetStep();

      return INPAR::STR::conv_fail_repeat;
      break;
    }
    case INPAR::STR::divcont_halve_step:
    {
      IO::cout << "Nonlinear solver failed to converge at time t= "<< GetTimeNp() << ". Divide timestep in half. "
               << "Old time step: " << GetDeltaTime() << IO::endl
               << "New time step: " << 0.5*GetDeltaTime() << IO::endl
               << IO::endl;

      // halve the time step size
      SetDeltaTime(GetDeltaTime()*0.5);
      // update the number of max time steps
      SetStepEnd(GetStepEnd() + (GetStepEnd() - GetStepNp()) +1);
      // reset timen_ because it is set in the constructor
      SetTimeNp(GetTimeN()+GetDeltaTime());
      // reset step (e.g. quantities on element level)
      ResetStep();

      return INPAR::STR::conv_fail_repeat;
      break;
    }
    case INPAR::STR::divcont_adapt_step:
    {
      // maximal possible refinementlevel
      const int maxstepmax = 1000000;
      IO::cout << "Nonlinear solver failed to converge at time t= "<< GetTimeNp() << ". Divide timestep in half. "
               << "Old time step: " << GetDeltaTime() << IO::endl
               << "New time step: " << 0.5*GetDeltaTime() << IO::endl
               << IO::endl;

      // halve the time step size
      SetDeltaTime(GetDeltaTime()*0.5);
      // update the number of max time steps
      SetStepEnd(GetStepEnd() + (GetStepEnd() - GetStepNp()) +1);
      // reset timen_ because it is set in the constructor
      SetTimeNp(GetTimeN()+GetDeltaTime());

      SetDivConRefineLevel(GetDivConRefineLevel()+1);
      SetDivConNumFineStep(0);

      if(GetDivConRefineLevel()==GetMaxDivConRefineLevel())
        dserror("Maximal divercont refinement level reached. Adapt your time basic time step size!");

      if(GetStepEnd()>maxstepmax)
        dserror("Upper level for stepmax_ reached!");

      // reset step (e.g. quantities on element level)
      ResetStep();

      return INPAR::STR::conv_fail_repeat;
      break;
    }
    case INPAR::STR::divcont_rand_adapt_step:
    case INPAR::STR::divcont_rand_adapt_step_ele_err:
    {
      // generate random number between 0.51 and 1.99 (as mean value of random
      // numbers generated on all processors), alternating values larger
      // and smaller than 1.0
      double proc_randnum_get = ((double) rand()/(double) RAND_MAX);
      double proc_randnum = proc_randnum_get;
      double randnum = 1.0;
      Discretization()->Comm().SumAll(&proc_randnum,&randnum,1);
      const double numproc = Discretization()->Comm().NumProc();
      randnum /= numproc;
      if      (GetRandomTimeStepFactor() > 1.0) SetRandomTimeStepFactor(randnum*0.49+0.51);
      else if (GetRandomTimeStepFactor() < 1.0) SetRandomTimeStepFactor(randnum*0.99+1.0);
      else                                      SetRandomTimeStepFactor(randnum*1.48+0.51);

      if (Discretization()->Comm().MyPID() == 0)
        IO::cout << "Nonlinear solver failed to converge: modifying time-step size by random number between 0.51 and 1.99 -> here: " << GetRandomTimeStepFactor() << " !" << IO::endl;
      // multiply time-step size by random number
      SetDeltaTime(GetDeltaTime()*GetRandomTimeStepFactor());
      // update maximum number of time steps
      SetStepEnd((1.0/GetRandomTimeStepFactor())*GetStepEnd() + (1.0-(1.0/GetRandomTimeStepFactor()))*GetStepNp() + 1);
      // reset timen_ because it is set in the constructor
      SetTimeNp(GetTimeN()+GetDeltaTime());
      // reset step (e.g. quantities on element level)
      ResetStep();

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
      if(nonlinsoldiv==INPAR::STR::conv_nonlin_fail)
        IO::cout << "Nonlinear solver failed to converge and DIVERCONT = "
            "repeat_simulation, hence leaving structural time integration "
            << IO::endl;
      else if (nonlinsoldiv==INPAR::STR::conv_lin_fail)
        IO::cout << "Linear solver failed to converge and DIVERCONT = "
            "repeat_simulation, hence leaving structural time integration "
            << IO::endl;
      else if (nonlinsoldiv==INPAR::STR::conv_ele_fail)
        IO::cout << "Element failure in form of a negative Jacobian determinant and DIVERCONT = "
            "repeat_simulation, hence leaving structural time integration "
            << IO::endl;
      return nonlinsoldiv; // so that time loop will be aborted
      break;
    }
    default:
      dserror("Unknown DIVER_CONT case");
    return INPAR::STR::conv_nonlin_fail;
    break;
  }
  return INPAR::STR::conv_success; // make compiler happy
} //PerformErrorAction()

/*-----------------------------------------------------------------------------*
 * check, if according to divercont flag                             meier 01/15
 * time step size can be increased
 *-----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::CheckForTimeStepIncrease(INPAR::STR::ConvergenceStatus& status)
{
  CheckInitSetup();

  const int maxnumfinestep = 4;

  if(GetDivergenceAction()!=INPAR::STR::divcont_adapt_step)
    return;
  else if(status == INPAR::STR::conv_success and GetDivConRefineLevel()!=0)
  {
    SetDivConNumFineStep(GetDivConNumFineStep()+1);

    if(GetDivConNumFineStep()==maxnumfinestep)
    {
      // increase the step size if the remaining number of steps is a even number
      if(((GetStepEnd() - GetStepNp())%2)==0 and GetStepEnd()!=GetStepNp())
      {
        IO::cout << "Nonlinear solver successful. Double timestep size!"
                 << IO::endl;

        SetDivConRefineLevel(GetDivConRefineLevel()-1);
        SetDivConNumFineStep(0);

        SetStepEnd(GetStepEnd() - (GetStepEnd() - GetStepNp())/2);

        // double the time step size
        SetDeltaTime(GetDeltaTime()*2.0);
      }
      else  //otherwise we have to wait one more time step until the step size can be increased
      {
        SetDivConNumFineStep(GetDivConNumFineStep()-1);
      }
    }
    return;
  }
} //CheckForTimeStepIncrease()
