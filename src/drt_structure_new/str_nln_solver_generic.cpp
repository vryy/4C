/*
 * str_nln_solver_generic.cpp
 *
 *  Created on: Sep 10, 2015
 *      Author: hiermeier
 */

#include "str_nln_solver_generic.H"

#include "str_timint_base.H"
#include "str_nln_solver_interface_required.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::NLN::SOLVER::Generic::Generic()
    : isinit_(false),
      issetup_(false),
      gstate_(Teuchos::null),
      sdyn_(Teuchos::null),
      ireq_(Teuchos::null)
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::Generic::Init(const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate,
    const Teuchos::RCP<STR::TIMINT::BaseDataSDyn> sdyn,
    const Teuchos::RCP<STR::NLN::SOLVER::INTERFACE::Required> ireq)
{
  // We have to call Setup() after Init()
  issetup_ = false;

  gstate_ = gstate;
  sdyn_ = sdyn;
  ireq_ = ireq;

  isinit_ = true;

  return;
}
