/*
 * str_timint_implicit.cpp
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#include "str_timint_implicit.H"
#include "str_impl_generic.H"
#include "str_impl_factory.H"
#include "str_predict_generic.H"
#include "str_nln_solver_generic.H"

// factories
#include "str_predict_factory.H"
#include "str_nln_solver_factory.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::Implicit::Implicit()
    : STR::TIMINT::Base(),
      implint_(Teuchos::null),
      nlnsolver_(Teuchos::null)
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::Setup()
{
  // safety check
  if (!IsInit())
    dserror("Init() has not been called, yet!");

  // ---------------------------------------------
  // build implicit integrator
  // ---------------------------------------------
  // get local references to the enums
  const enum INPAR::STR::DynamicType& dynType =
      DataSDyn().GetDynamicType();
  const enum INPAR::STR::PreStress& preStressType =
      DataSDyn().GetPreStressType();
  implint_ = STR::IMPLICIT::BuildImplicitIntegrator(dynType,preStressType);
  implint_->Init();
  implint_->Setup();

  // ---------------------------------------------
  // build predictor
  // ---------------------------------------------
  const enum INPAR::STR::PredEnum& predtype =
      DataSDyn().GetPredictorType();
  predictor_ = STR::PREDICT::BuildPredictor(predtype);
  predictor_->Init(predtype);
  predictor_->Setup();

  // ---------------------------------------------
  // build non-linear solver
  // ---------------------------------------------
  const enum INPAR::STR::NonlinSolTech& nlnSolverType =
      DataSDyn().GetNlnSolverType();
  nlnsolver_ = STR::NLN::SOLVER::BuildNlnSolver(nlnSolverType);
  nlnsolver_->Init(DataGlobalStatePtr(), DataSDynPtr(), Teuchos::rcp(this,false));
  nlnsolver_->Setup();

  // set isSetup flag
  issetup_ = true;

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Implicit::PrepareTimeStep()
{
  // update end time \f$t_{n+1}\f$ of this time step to cope with time step size adaptivity
  DataGlobalStatePtr()->SetTimenp(DataGlobalStatePtr()->GetUpdatedTime());

  // things that need to be done before Predict
  PrePredict();

  // TODO prepare contact for new time step
  // PrepareStepContact();
  Predictor().Predict();

  // things that need to be done after Predict
  PostPredict();

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Implicit::IntegrateStep()
{
  // do the predictor step
  Predictor().Predict();
  return Solve();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
INPAR::STR::ConvergenceStatus STR::TIMINT::Implicit::Solve()
{
  // reset the non-linear solver
  NlnSolver().Reset();
  // solve the non-linear problem
  return NlnSolver().Solve();
}
