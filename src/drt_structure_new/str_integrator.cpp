/*-----------------------------------------------------------*/
/*!
\file str_integrator.cpp

\maintainer Michael Hiermeier

\date Dec 7, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_integrator.H"
#include "../drt_lib/drt_dserror.H"

#include "str_timint_base.H"
#include "str_model_evaluator.H"

#include <Epetra_Vector.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::Integrator::Integrator()
    : isinit_(false),
      issetup_(false),
      modelevaluator_ptr_(Teuchos::null),
      gstate_ptr_(Teuchos::null),
      timint_ptr_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::Init(const Teuchos::RCP<STR::ModelEvaluator>& me_ptr)
{
  issetup_ = false;

  modelevaluator_ptr_ = me_ptr;
  gstate_ptr_ = me_ptr->GlobalStatePtr();
  timint_ptr_ = me_ptr->GetTimIntPtr();

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::CheckInit() const
{
  if (not IsInit())
    dserror("Call Init() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Integrator::CheckInitSetup() const
{
  if (not IsInit() or not IsSetup())
    dserror("Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::ModelEvaluator& STR::Integrator::ModelEval()
{
  CheckInit();
  if (modelevaluator_ptr_.is_null())
    dserror("The modelevaluator is not initialized!");

  return *modelevaluator_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::ModelEvaluator& STR::Integrator::ModelEval() const
{
  CheckInit();
  if (modelevaluator_ptr_.is_null())
    dserror("The modelevaluator is not initialized!");

  return *modelevaluator_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::BaseDataGlobalState& STR::Integrator::GlobalState() const
{
  CheckInit();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataGlobalState& STR::Integrator::GlobalState()
{
  CheckInit();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::Base& STR::Integrator::TimInt() const
{
  CheckInit();
  return *timint_ptr_;
}
