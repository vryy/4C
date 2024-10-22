// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_solver_nonlin_nox_linesearch_prepostoperator.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::LineSearch::PrePostOperator::PrePostOperator(Teuchos::ParameterList& linesearchSublist)
{
  reset(linesearchSublist);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::LineSearch::PrePostOperator::reset(Teuchos::ParameterList& linesearchSublist)
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
NOX::Nln::LineSearch::PrePostOperator::map& NOX::Nln::LineSearch::PrePostOperator::get_map(
    Teuchos::ParameterList& p_ls_list)
{
  Teuchos::RCP<map>& mapptr =
      p_ls_list.get<Teuchos::RCP<map>>("User Defined Pre/Post Operator", Teuchos::make_rcp<map>());

  return *mapptr;
}

FOUR_C_NAMESPACE_CLOSE
