/*
 * str_timint_implicit.cpp
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#include "str_timint_implicit.H"
#include "str_nln_solver_factory.H"
#include "str_predict_generic.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::Implicit::Implicit()
    : implint_(Teuchos::null),
      modelevaluators_(Teuchos::null),
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
  // build model evaluator
  // ---------------------------------------------
  modelevaluators_ = STR::MODELEVALUATOR::BuildModelEvaluators(DataSDyn().GetModelType());

  std::map<const enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> >::iterator me_iter;
  for (me_iter=modelevaluators_->begin();me_iter!=modelevaluators_->end();++me_iter)
  {
    me_iter->second->Init();
    me_iter->second>Setup();
  }

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
  nlnsolver_->Init(DataSDynPtr(), DataGlobalStatePtr(), Teuchos::rcp(this,false));
  nlnsolver_->Setup();

  // set isSetup flag
  issetup_ = true;

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Implicit::IntegrateStep()
{
  PreSolve();

  // ---------------------------------------------
  // Give a first educated guess
  // ---------------------------------------------
  Predictor().Predict();

  // ---------------------------------------------
  // Solve the nonlinear problem
  // ---------------------------------------------
  NlnSolver().Reset();
  int error = NlnSolver().Solve();

  PostSolve();

  return error;
}
