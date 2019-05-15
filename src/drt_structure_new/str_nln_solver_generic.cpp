/*-----------------------------------------------------------*/
/*!

\brief Generic class of the non-linear structural solvers.

\maintainer Anh-Tu Vuong

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_nln_solver_generic.H"
#include "str_timint_implicit.H"
#include "str_timint_base.H"
#include "str_timint_noxinterface.H"

#include <NOX_Abstract_Group.H>


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

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Abstract::Group>& STR::NLN::SOLVER::Generic::GroupPtr()
{
  CheckInit();

  return group_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group& STR::NLN::SOLVER::Generic::Group()
{
  CheckInit();
  if (group_ptr_.is_null()) dserror("The group pointer should be initialized beforehand!");
  return *group_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group& STR::NLN::SOLVER::Generic::SolutionGroup() { return Group(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const NOX::Abstract::Group& STR::NLN::SOLVER::Generic::GetSolutionGroup() const
{
  CheckInitSetup();
  if (group_ptr_.is_null()) dserror("The group pointer should be initialized beforehand!");

  return *group_ptr_;
}
