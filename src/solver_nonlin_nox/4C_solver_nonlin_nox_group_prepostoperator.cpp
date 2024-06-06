/*-----------------------------------------------------------*/
/*! \file

\brief wrapper class for a user derived NOX PrePostOperator



\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_solver_nonlin_nox_group_prepostoperator.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::GROUP::PrePostOperator::PrePostOperator() : havePrePostOperator_(false)
{
  // Disallowed constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::GROUP::PrePostOperator::PrePostOperator(const PrePostOperator& ppo)
    : havePrePostOperator_(false)
{
  // Disallowed copy constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::GROUP::PrePostOperator& NOX::Nln::GROUP::PrePostOperator::operator=(
    const PrePostOperator& ppo)
{
  // disallowed assignment operator
  return *this;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::GROUP::PrePostOperator::PrePostOperator(Teuchos::ParameterList& groupOptionsSubList)
    : havePrePostOperator_(false), prePostOperatorMapPtr_(Teuchos::null)
{
  reset(groupOptionsSubList);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::GROUP::PrePostOperator::reset(Teuchos::ParameterList& groupOptionsSubList)
{
  havePrePostOperator_ = false;

  /* Check if a pre/post operator for the group is provided
   * by the user. */
  if (groupOptionsSubList.INVALID_TEMPLATE_QUALIFIER isType<Teuchos::RCP<Map>>(
          "User Defined Pre/Post Operator"))
  {
    prePostOperatorMapPtr_ = groupOptionsSubList.INVALID_TEMPLATE_QUALIFIER get<Teuchos::RCP<Map>>(
        "User Defined Pre/Post Operator");
    havePrePostOperator_ = true;
  }
}

// non-member function
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::GROUP::PrePostOperator::Map& NOX::Nln::GROUP::PrePostOp::GetMap(
    Teuchos::ParameterList& p_grp_opt)
{
  Teuchos::RCP<NOX::Nln::GROUP::PrePostOperator::Map>& mapptr =
      p_grp_opt.get<Teuchos::RCP<NOX::Nln::GROUP::PrePostOperator::Map>>(
          "User Defined Pre/Post Operator",
          Teuchos::rcp(new NOX::Nln::GROUP::PrePostOperator::Map()));

  return *mapptr;
}

FOUR_C_NAMESPACE_CLOSE
