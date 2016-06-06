/*---------------------------------------------------------------------*/
/*!
\file str_model_evaluator_contactdata.cpp

\brief Concrete implementation of the contact parameter interfaces.

\maintainer Michael Hiermeier

\date Apr 18, 2016

\level 3

*/
/*---------------------------------------------------------------------*/


#include "str_model_evaluator_data.H"
#include "str_timint_implicit.H"
#include "str_nln_solver_nox.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::ContactData::ContactData()
    : isinit_(false),
      issetup_(false),
      mortar_action_(MORTAR::eval_none),
      str_data_ptr_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::ContactData::Init(
    const Teuchos::RCP<const STR::MODELEVALUATOR::Data>& str_data_ptr)
{
  issetup_ = false;
  str_data_ptr_ = str_data_ptr;
  isinit_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::ContactData::Setup()
{
  CheckInit();

  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const bool& STR::MODELEVALUATOR::ContactData::IsPredictor() const
{
  CheckInit();
  return GetGState().IsPredict();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int STR::MODELEVALUATOR::ContactData::GetNlnIter() const
{
  if (IsPredictor())
    return 0;

  bool isnox = false;
  Teuchos::RCP<const STR::NLN::SOLVER::Nox> nox_nln_ptr = Teuchos::null;
  const STR::TIMINT::Implicit* timint_impl_ptr =
      dynamic_cast<const STR::TIMINT::Implicit*>(&GetTimInt());
  if (timint_impl_ptr!=NULL)
  {
    nox_nln_ptr =
        Teuchos::rcp_dynamic_cast<const STR::NLN::SOLVER::Nox>(
            timint_impl_ptr->GetNlnSolverPtr());
    if (not nox_nln_ptr.is_null())
      isnox = true;
  }
  if (not isnox)
    dserror("The GetNlnIter() routine supports only the NOX::NLN "
        "framework at the moment.");

  return nox_nln_ptr->GetNumNlnIterations();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int STR::MODELEVALUATOR::ContactData::GetStepNp() const
{
  CheckInit();
  return GetGState().GetStepNp();
}
