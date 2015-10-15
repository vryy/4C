/*-----------------------------------------------------------*/
/*!
\file str_nln_solver_generic.cpp

\maintainer Michael Hiermeier

\date Oct 9, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_nln_solver_generic.H"
#include "str_timint_implicit.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::NLN::SOLVER::Generic::Generic()
    : isinit_(false),
      issetup_(false),
      gstate_(Teuchos::null),
      sdyn_(Teuchos::null),
      iimpl_(Teuchos::null)
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::Generic::Init(
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate,
    const Teuchos::RCP<STR::TIMINT::BaseDataSDyn> sdyn,
    const Teuchos::RCP<STR::TIMINT::Implicit> implicit)
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // initialize internal variables
  gstate_ = gstate;
  sdyn_ = sdyn;
  iimpl_ = implicit;

  isinit_ = true;

  return;
}
