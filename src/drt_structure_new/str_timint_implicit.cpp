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
  return NlnSolver().Solve();
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
