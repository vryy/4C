/*
 * str_timint_implicit.cpp
 *
 *  Created on: Aug 13, 2015
 *      Author: farah
 */

#include "str_timint_implicit.H"
#include "str_nln_solver_factory.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::Implicit::Implicit() : STR::TIMINT::Base(),
  nlnSolver_(Teuchos::null)
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
  implInt_ = STR::IMPLICIT::BuildImplicitIntegrator(dynType,preStressType);

  // ---------------------------------------------
  // build model evaluator
  // ---------------------------------------------
  modelEvalMap_ = STR::MODELEVALUATOR::BuildModelEvaluator();

  // ---------------------------------------------
  // build nonlinear solver
  // ---------------------------------------------
  const enum INPAR::STR::NonlinSolTech& nlnSolverType =
      DataSDyn().GetNlnSolverType();
  nlnSolver_ = STR::NLN::SOLVER::BuildNlnSolver(nlnSolverType);

  // set isSetup flag
  isSetup_ = true;

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void  STR::TIMINT::Implicit::IntegrateStep()
{
  return;
}
