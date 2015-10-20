/*-----------------------------------------------------------*/
/*!
\file nox_nln_linearsystem_prepostoperator.cpp

\maintainer Michael Hiermeier

\date Oct 7, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_linearsystem_prepostoperator.H"

#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::LinSystem::PrePostOperator::PrePostOperator()
    : havePrePostOperator_(false)
{
  // Disallowed constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::LinSystem::PrePostOperator::PrePostOperator(const PrePostOperator& ppo)
    : havePrePostOperator_(false)
{
  // Disallowed copy constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::LinSystem::PrePostOperator& NOX::NLN::LinSystem::PrePostOperator::
operator=(const PrePostOperator& ppo)
{
  // disallowed assignment operator
  return *this;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::LinSystem::PrePostOperator::PrePostOperator(
    Teuchos::ParameterList& linearSolverSubList)
    : havePrePostOperator_(false)
{
  reset(linearSolverSubList);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::LinSystem::PrePostOperator::~PrePostOperator()
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::LinSystem::PrePostOperator::reset(
    Teuchos::ParameterList& linearSolverSubList)
{
  havePrePostOperator_ = false;

  /* Check if a pre/post processor for the linear system is provided
   * by the user.
   */
  if (linearSolverSubList.INVALID_TEMPLATE_QUALIFIER
      isType< Teuchos::RCP<NOX::NLN::LinSystem::PrePostOp::Generic> >
      ("User Defined Pre/Post Operator"))
  {
    prePostOperatorPtr_ = linearSolverSubList.INVALID_TEMPLATE_QUALIFIER
      get< Teuchos::RCP<NOX::NLN::LinSystem::PrePostOp::Generic> >
      ("User Defined Pre/Post Operator");
    havePrePostOperator_ = true;
  }
}
