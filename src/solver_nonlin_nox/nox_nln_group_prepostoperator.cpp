/*
 * nox_nln_group_prepostoperator.cpp
 *
 *  Created on: Oct 19, 2015
 *      Author: hiermeier
 */

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_GROUP_PREPOSTOPERATOR_CPP_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_GROUP_PREPOSTOPERATOR_CPP_


#include "nox_nln_group_prepostoperator.H"

#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::GROUP::PrePostOperator::PrePostOperator()
    : havePrePostOperator_(false)
{
  // Disallowed constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::GROUP::PrePostOperator::PrePostOperator(const PrePostOperator& ppo)
    : havePrePostOperator_(false)
{
  // Disallowed copy constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::GROUP::PrePostOperator& NOX::NLN::GROUP::PrePostOperator::
operator=(const PrePostOperator& ppo)
{
  // disallowed assignment operator
  return *this;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::GROUP::PrePostOperator::PrePostOperator(
    Teuchos::ParameterList& groupOptionsSubList)
    : havePrePostOperator_(false)
{
  reset(groupOptionsSubList);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::GROUP::PrePostOperator::reset(
    Teuchos::ParameterList& groupOptionsSubList)
{
  havePrePostOperator_ = false;

  /* Check if a pre/post processor for the linear system is provided
   * by the user.
   */
  if (groupOptionsSubList.INVALID_TEMPLATE_QUALIFIER
      isType< Teuchos::RCP<NOX::NLN::Abstract::PrePostOperator> >
      ("User Defined Pre/Post Operator"))
  {
    prePostOperatorPtr_ = groupOptionsSubList.INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RCP<NOX::NLN::Abstract::PrePostOperator> >
      ("User Defined Pre/Post Operator");
    havePrePostOperator_ = true;
  }
}


#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_GROUP_PREPOSTOPERATOR_CPP_ */
