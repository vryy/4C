/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_multiphysics.cpp

\brief Base class for modelevaluators in partitioned algorithms.

\maintainer Andreas Rauch

\date Nov 28, 2016

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_model_evaluator_multiphysics.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Multiphysics::Multiphysics() : active_mt_(mt_none)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Multiphysics::Init(
    const Teuchos::RCP<STR::MODELEVALUATOR::Data>& eval_data_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
    const Teuchos::RCP<STR::Integrator>& int_ptr,
    const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr, const int& dof_offset)
{
  STR::MODELEVALUATOR::Generic::Init(
      eval_data_ptr, gstate_ptr, gio_ptr, int_ptr, timint_ptr, dof_offset);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Multiphysics::Setup()
{
  CheckInit();
  issetup_ = true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Multiphysics::Reset(const Epetra_Vector& x)
{
  CheckActiveModelType();

  GetModelEvaluatorFromMap(active_mt_)->Reset(x);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Multiphysics::AssembleForce(
    Epetra_Vector& f, const double& timefac_np) const
{
  CheckActiveModelType();

  GetModelEvaluatorFromMap(active_mt_)->AssembleForce(f, timefac_np);
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Multiphysics::AssembleJacobian(
    LINALG::SparseOperator& jac, const double& timefac_np) const
{
  CheckActiveModelType();

  GetModelEvaluatorFromMap(active_mt_)->AssembleJacobian(jac, timefac_np);
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Multiphysics::EvaluateForce()
{
  CheckActiveModelType();

  GetModelEvaluatorFromMap(active_mt_)->EvaluateForce();
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Multiphysics::EvaluateStiff()
{
  CheckActiveModelType();

  GetModelEvaluatorFromMap(active_mt_)->EvaluateStiff();
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Multiphysics::EvaluateForceStiff()
{
  CheckActiveModelType();

  GetModelEvaluatorFromMap(active_mt_)->EvaluateForceStiff();
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Multiphysics::UpdateStepState(const double& timefac_n)
{
  CheckActiveModelType();

  GetModelEvaluatorFromMap(active_mt_)->UpdateStepState(timefac_n);
}
