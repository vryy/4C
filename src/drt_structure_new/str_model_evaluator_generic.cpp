/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_generic.cpp

\brief Generic class for all model evaluators.

\maintainer Michael Hiermeier

\date Dec 1, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_model_evaluator_generic.H"
#include "str_timint_base.H"
#include "str_model_evaluator_data.H"

#include "../drt_lib/drt_dserror.H"

#include <Epetra_Comm.h>
#include <fenv.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Generic::Generic()
    : isinit_(false),
      issetup_(false),
      eval_data_ptr_(Teuchos::null),
      gstate_ptr_(Teuchos::null),
      gio_ptr_(Teuchos::null),
      discret_ptr_(Teuchos::null),
      timint_ptr_(Teuchos::null),
      dof_offset_(0)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Generic::Init(
    const Teuchos::RCP<STR::MODELEVALUATOR::Data>& eval_data_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
    const Teuchos::RCP<STR::Integrator>& int_ptr,
    const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr,
    const int& dof_offset)
{
  // call setup after init()
  issetup_ = false;

  eval_data_ptr_ = eval_data_ptr;
  gstate_ptr_    = gstate_ptr;
  gio_ptr_       = gio_ptr;
  discret_ptr_   = gstate_ptr->GetMutableDiscret();
  int_ptr_       = int_ptr;
  timint_ptr_    = timint_ptr;
  dof_offset_    = dof_offset;

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
STR::MODELEVALUATOR::Data& STR::MODELEVALUATOR::Generic::EvalData()
{
  CheckInit();
  return *eval_data_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::MODELEVALUATOR::Data>& STR::MODELEVALUATOR::Generic::EvalDataPtr()
{
  CheckInit();
  return eval_data_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataGlobalState& STR::MODELEVALUATOR::Generic::GState()
{
  CheckInit();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& STR::MODELEVALUATOR::Generic::GStatePtr()
{
  CheckInit();
  return gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::BaseDataGlobalState& STR::MODELEVALUATOR::Generic::GState()
    const
{
  CheckInit();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataIO& STR::MODELEVALUATOR::Generic::GInOutput()
{
  CheckInit();
  return *gio_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::BaseDataIO& STR::MODELEVALUATOR::Generic::GInOutput() const
{
  CheckInit();
  return *gio_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::Discretization& STR::MODELEVALUATOR::Generic::Discret()
{
  CheckInit();
  return *discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization>& STR::MODELEVALUATOR::Generic::DiscretPtr()
{
  CheckInit();
  return discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const DRT::Discretization& STR::MODELEVALUATOR::Generic::Discret() const
{
  CheckInit();
  return *discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::Integrator& STR::MODELEVALUATOR::Generic::Int()
{
  CheckInit();
  return *int_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::Base& STR::MODELEVALUATOR::Generic::TimInt() const
{
  CheckInit();
  return *timint_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const int& STR::MODELEVALUATOR::Generic::DofOffset() const
{
  CheckInit();
  return dof_offset_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Generic::EvalErrorCheck() const
{
  // --- Did an exception occur during the evaluation process? -----------------
  bool ok = true;
  if (fetestexcept(FE_INVALID)    or
      fetestexcept(FE_OVERFLOW)   or
      fetestexcept(FE_DIVBYZERO))
    ok = false;

  // --- Did the element evaluation detect an error? ---------------------------
  ok = (ok and (not eval_data_ptr_->IsEleEvalError()));
  eval_data_ptr_->SetEleEvalErrorFlag(STR::ELEMENTS::ele_error_none);

  // --- check for local errors on each proc and communicate the information ---
  int lerr = (ok ? 0 : 1);
  int gerr = 0;
  gstate_ptr_->GetComm().SumAll(&lerr,&gerr,1);
  return (gerr==0);
}
