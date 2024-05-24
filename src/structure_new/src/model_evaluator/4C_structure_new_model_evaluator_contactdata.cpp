/*---------------------------------------------------------------------*/
/*! \file

\brief Concrete implementation of the contact parameter interfaces.


\level 3

*/
/*---------------------------------------------------------------------*/


#include "4C_inpar_contact.hpp"
#include "4C_structure_new_integrator.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_nln_solver_nox.hpp"
#include "4C_structure_new_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::ContactData::ContactData()
    : isinit_(false),
      issetup_(false),
      mortar_action_(MORTAR::eval_none),
      var_type_(INPAR::CONTACT::var_unknown),
      coupling_scheme_(INPAR::CONTACT::CouplingScheme::unknown),
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
  check_init();

  issetup_ = true;
}

FOUR_C_NAMESPACE_CLOSE
