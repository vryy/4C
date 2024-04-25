/*-----------------------------------------------------------*/
/*! \file

\brief Generic class of the non-linear structural solvers.


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_nln_solver_generic.hpp"

#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_timint_implicit.hpp"
#include "4C_structure_new_timint_noxinterface.hpp"

#include <NOX_Abstract_Group.H>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::NLN::SOLVER::Generic::Generic()
    : isinit_(false),
      issetup_(false),
      gstate_ptr_(Teuchos::null),
      sdyn_ptr_(Teuchos::null),
      noxinterface_ptr_(Teuchos::null),
      group_ptr_(Teuchos::null)
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::Generic::Init(const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate,
    const Teuchos::RCP<STR::TIMINT::BaseDataSDyn>& sdyn,
    const Teuchos::RCP<STR::TIMINT::NoxInterface>& noxinterface,
    const Teuchos::RCP<STR::Integrator>& integrator,
    const Teuchos::RCP<const STR::TIMINT::Base>& timint)
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // initialize internal variables
  gstate_ptr_ = gstate;
  sdyn_ptr_ = sdyn;
  noxinterface_ptr_ = noxinterface;
  int_ptr_ = integrator;
  timint_ptr_ = timint;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Abstract::Group>& STR::NLN::SOLVER::Generic::GroupPtr()
{
  CheckInit();

  return group_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group& STR::NLN::SOLVER::Generic::Group()
{
  CheckInit();
  FOUR_C_ASSERT(!group_ptr_.is_null(), "The group pointer should be initialized beforehand!");
  return *group_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group& STR::NLN::SOLVER::Generic::SolutionGroup() { return Group(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const ::NOX::Abstract::Group& STR::NLN::SOLVER::Generic::GetSolutionGroup() const
{
  CheckInitSetup();
  FOUR_C_ASSERT(!group_ptr_.is_null(), "The group pointer should be initialized beforehand!");

  return *group_ptr_;
}

FOUR_C_NAMESPACE_CLOSE
