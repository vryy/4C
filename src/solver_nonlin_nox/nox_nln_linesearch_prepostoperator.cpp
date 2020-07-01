/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_linesearch_prepostoperator.H"
#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::LineSearch::PrePostOperator::PrePostOperator(Teuchos::ParameterList& linesearchSublist)
{
  reset(linesearchSublist);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::LineSearch::PrePostOperator::reset(Teuchos::ParameterList& linesearchSublist)
{
  havePrePostOperator_ = false;

  /* Check if a pre/post operator for linesearch is provided
   * by the user. */
  if (linesearchSublist.INVALID_TEMPLATE_QUALIFIER isType<Teuchos::RCP<map>>(
          "User Defined Pre/Post Operator"))
  {
    prePostOperatorMapPtr_ = linesearchSublist.INVALID_TEMPLATE_QUALIFIER get<Teuchos::RCP<map>>(
        "User Defined Pre/Post Operator");
    havePrePostOperator_ = true;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::LineSearch::PrePostOperator::map& NOX::NLN::LineSearch::PrePostOperator::GetMutableMap(
    Teuchos::ParameterList& p_ls_list)
{
  Teuchos::RCP<map>& mapptr =
      p_ls_list.get<Teuchos::RCP<map>>("User Defined Pre/Post Operator", Teuchos::rcp(new map()));

  return *mapptr;
}
