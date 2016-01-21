/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_generic.cpp

\maintainer Michael Hiermeier

\date Dec 1, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_model_evaluator_generic.H"
#include "str_timint_base.H"

#include "../drt_lib/drt_dserror.H"

#include <fenv.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Generic::Generic()
    : isinit_(false),
      issetup_(false),
      gstate_ptr_(Teuchos::null),
      discret_ptr_(Teuchos::null),
      timint_ptr_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Generic::Init(
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
    const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr)
{
  // call setup after init()
  issetup_ = false;

  gstate_ptr_ = gstate_ptr;
  discret_ptr_ = gstate_ptr->GetMutableDiscret();
  timint_ptr_ = timint_ptr;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Generic::CheckInitSetup() const
{
  if (!IsInit() or !IsSetup())
    dserror("Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Generic::CheckInit() const
{
  if (not IsInit())
    dserror("Call Init() first!");
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Generic::EvalErrorCheck(
    const Teuchos::ParameterList& p) const
{
  // Did an exception occur during the evaluation process?
  bool ok = true;
  if (fetestexcept(FE_INVALID)    or
      fetestexcept(FE_OVERFLOW)   or
      fetestexcept(FE_DIVBYZERO))
    ok = false;

  // Did the element evaluation detect an error?
  if (ok and p.isParameter("eval_error"))
   ok = p.get<bool>("eval_error");

  return ok;
}
