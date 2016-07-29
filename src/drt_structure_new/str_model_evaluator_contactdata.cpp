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
#include "str_integrator.H"
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
