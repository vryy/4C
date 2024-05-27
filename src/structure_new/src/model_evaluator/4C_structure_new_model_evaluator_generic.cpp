/*-----------------------------------------------------------*/
/*! \file

\brief Generic class for all model evaluators.


\level 3

*/
/*-----------------------------------------------------------*/


#include "4C_structure_new_model_evaluator_generic.hpp"

#include "4C_solver_nonlin_nox_floating_point_exception.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_utils.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Comm.h>

FOUR_C_NAMESPACE_OPEN

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
  discret_ptr_ = gstate_ptr->GetDiscret();
  int_ptr_ = int_ptr;
  timint_ptr_ = timint_ptr;
  dof_offset_ = dof_offset;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Generic::check_init_setup() const
{
  FOUR_C_ASSERT(is_init() and is_setup(), "Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Generic::check_init() const
{
  FOUR_C_ASSERT(is_init(), "Call Init() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Data& STR::MODELEVALUATOR::Generic::eval_data()
{
  check_init();
  return *eval_data_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::MODELEVALUATOR::Data& STR::MODELEVALUATOR::Generic::EvalData() const
{
  check_init();
  return *eval_data_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::MODELEVALUATOR::Data>& STR::MODELEVALUATOR::Generic::eval_data_ptr()
{
  check_init();
  return eval_data_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataGlobalState& STR::MODELEVALUATOR::Generic::g_state()
{
  check_init();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& STR::MODELEVALUATOR::Generic::g_state_ptr()
{
  check_init();
  return gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::BaseDataGlobalState& STR::MODELEVALUATOR::Generic::GState() const
{
  check_init();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataIO& STR::MODELEVALUATOR::Generic::g_in_output()
{
  check_init();
  return *gio_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TIMINT::BaseDataIO> STR::MODELEVALUATOR::Generic::g_in_output_ptr()
{
  check_init();
  return gio_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::BaseDataIO& STR::MODELEVALUATOR::Generic::GInOutput() const
{
  check_init();
  return *gio_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::Discretization& STR::MODELEVALUATOR::Generic::discret()
{
  check_init();
  return *discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization>& STR::MODELEVALUATOR::Generic::discret_ptr()
{
  check_init();
  return discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const DRT::Discretization& STR::MODELEVALUATOR::Generic::Discret() const
{
  check_init();
  return *discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::Integrator& STR::MODELEVALUATOR::Generic::integrator()
{
  check_init();
  return *int_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::Integrator& STR::MODELEVALUATOR::Generic::Int() const
{
  check_init();
  return *int_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::Integrator>& STR::MODELEVALUATOR::Generic::int_ptr() { return int_ptr_; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::Base& STR::MODELEVALUATOR::Generic::TimInt() const
{
  check_init();
  return *timint_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const int& STR::MODELEVALUATOR::Generic::dof_offset() const
{
  check_init();
  return dof_offset_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> STR::MODELEVALUATOR::Generic::GetFextIncr() const
{
  check_init_setup();
  const Epetra_Vector& fextn = *GState().GetFextN();
  const Epetra_Vector& fextnp = *GState().GetFextNp();

  Teuchos::RCP<Epetra_Vector> fext_incr = Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(fextnp));
  fext_incr->Update(-1.0, fextn, 1.0);

  return fext_incr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Generic::eval_error_check() const
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

FOUR_C_NAMESPACE_CLOSE
