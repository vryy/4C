/*-----------------------------------------------------------*/
/*! \file

\brief Generic class for all model evaluators.

\maintainer Matthias Mayr

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_model_evaluator_generic.H"
#include "str_timint_base.H"
#include "str_model_evaluator_data.H"
#include "str_utils.H"

#include "../drt_lib/drt_dserror.H"

#include <Epetra_Comm.h>
#include "../solver_nonlin_nox/nox_nln_floating_point_exception.H"

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
    const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr, const int& dof_offset)
{
  // call setup after init()
  issetup_ = false;

  eval_data_ptr_ = eval_data_ptr;
  gstate_ptr_ = gstate_ptr;
  gio_ptr_ = gio_ptr;
  discret_ptr_ = gstate_ptr->GetMutableDiscret();
  int_ptr_ = int_ptr;
  timint_ptr_ = timint_ptr;
  dof_offset_ = dof_offset;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Generic::CheckInitSetup() const
{
  if (!IsInit() or !IsSetup()) dserror("Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Generic::CheckInit() const
{
  if (not IsInit()) dserror("Call Init() first!");
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
const STR::MODELEVALUATOR::Data& STR::MODELEVALUATOR::Generic::EvalData() const
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
const STR::TIMINT::BaseDataGlobalState& STR::MODELEVALUATOR::Generic::GState() const
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
Teuchos::RCP<STR::TIMINT::BaseDataIO> STR::MODELEVALUATOR::Generic::GInOutputPtr()
{
  CheckInit();
  return gio_ptr_;
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
DRT::DiscretizationInterface& STR::MODELEVALUATOR::Generic::Discret()
{
  CheckInit();
  return *discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::DiscretizationInterface>& STR::MODELEVALUATOR::Generic::DiscretPtr()
{
  CheckInit();
  return discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const DRT::DiscretizationInterface& STR::MODELEVALUATOR::Generic::Discret() const
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
const STR::Integrator& STR::MODELEVALUATOR::Generic::Int() const
{
  CheckInit();
  return *int_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::Integrator>& STR::MODELEVALUATOR::Generic::IntPtr() { return int_ptr_; }

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
Teuchos::RCP<Epetra_Vector> STR::MODELEVALUATOR::Generic::GetFextIncr() const
{
  CheckInitSetup();
  const Epetra_Vector& fextn = *GState().GetFextN();
  const Epetra_Vector& fextnp = *GState().GetFextNp();

  Teuchos::RCP<Epetra_Vector> fext_incr = Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(fextnp));
  fext_incr->Update(-1.0, fextn, 1.0);

  return fext_incr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Generic::EvalErrorCheck() const
{
  // --- Did an exception occur during the evaluation process? -----------------
  bool ok = true;
  int fp_err = NOX::NLN::FloatingPointException::checkAndPrint(std::cout);
  if (fp_err)
  {
    ok = false;
    std::cout << "FLOATING POINT EXCEPTION occurred on proc "
                 "#"
              << gstate_ptr_->GetComm().MyPID() << ".\n";
  }

  // --- Did the element evaluation detect an error? ---------------------------
  ok = (ok and (not eval_data_ptr_->IsEleEvalError()));

  if (eval_data_ptr_->IsEleEvalError())
    std::cout << "ELEMENT EVALUATION failed on proc "
                 "#"
              << gstate_ptr_->GetComm().MyPID() << ".\n"
              << "(Error: "
              << STR::ELEMENTS::EvalErrorFlag2String(eval_data_ptr_->GetEleEvalErrorFlag())
              << ")\n";

  // reset the flag
  eval_data_ptr_->SetEleEvalErrorFlag(STR::ELEMENTS::ele_error_none);

  // --- check for local errors on each proc and communicate the information ---
  int lerr = (ok ? 0 : 1);
  int gerr = 0;
  gstate_ptr_->GetComm().SumAll(&lerr, &gerr, 1);
  return (gerr == 0);
}
