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
Solid::Nln::SOLVER::Generic::Generic()
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
void Solid::Nln::SOLVER::Generic::init(
    const Teuchos::RCP<Solid::TimeInt::BaseDataGlobalState>& gstate,
    const Teuchos::RCP<Solid::TimeInt::BaseDataSDyn>& sdyn,
    const Teuchos::RCP<Solid::TimeInt::NoxInterface>& noxinterface,
    const Teuchos::RCP<Solid::Integrator>& integrator,
    const Teuchos::RCP<const Solid::TimeInt::Base>& timint)
{
  // We have to call setup() after init()
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
Teuchos::RCP<::NOX::Abstract::Group>& Solid::Nln::SOLVER::Generic::group_ptr()
{
  check_init();

  return group_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group& Solid::Nln::SOLVER::Generic::group()
{
  check_init();
  FOUR_C_ASSERT(!group_ptr_.is_null(), "The group pointer should be initialized beforehand!");
  return *group_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group& Solid::Nln::SOLVER::Generic::SolutionGroup() { return group(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const ::NOX::Abstract::Group& Solid::Nln::SOLVER::Generic::get_solution_group() const
{
  check_init_setup();
  FOUR_C_ASSERT(!group_ptr_.is_null(), "The group pointer should be initialized beforehand!");

  return *group_ptr_;
}

FOUR_C_NAMESPACE_CLOSE
