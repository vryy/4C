/*----------------------------------------------------------------------------*/
/*!
\brief This file contains the adaptions of the structural time integration due
       to XFEM

\level 2

\maintainer Anh-Tu Vuong

*/
/*----------------------------------------------------------------------------*/

#include "../drt_structure_new/str_timint_base.H"
#include "../drt_xfem/xfield_state.H"

#include "../drt_structure_new/str_timint_factory.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<XFEM::XFieldState> STR::TIMINT::Base::XFieldState()
{
  return Teuchos::rcp_dynamic_cast<XFEM::XFieldState>(dataglobalstate_, true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const XFEM::XFieldState> STR::TIMINT::Base::XFieldState() const
{
  return Teuchos::rcp_dynamic_cast<const XFEM::XFieldState>(dataglobalstate_, true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::CreateNewXFieldState(
    const Teuchos::RCP<XFEM::XFieldState>& new_xstate) const
{
  // get the current global state object
  Teuchos::RCP<const XFEM::XFieldState> curr_xstate = XFieldState();

  // get the new discretization
  Teuchos::RCP<const DRT::DiscretizationInterface> full_discret = dataglobalstate_->GetDiscret();

  // build a new xstate object
  Teuchos::RCP<BaseDataGlobalState> new_gstate =
      Teuchos::rcp_dynamic_cast<BaseDataGlobalState>(new_xstate, true);
  *new_gstate = *dataglobalstate_;
  new_gstate->Setup();

  curr_xstate->TransferToNewState(*full_discret, *new_xstate);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TIMINT::Base::DestroyState()
{
  Teuchos::RCP<XFEM::XFieldState> xstate = XFieldState();

  DestroyNoxState();

  xstate->Destroy();

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::ResetXFieldNonStandardDofs()
{
  Teuchos::RCP<const DRT::DiscretizationInterface> full_discret = dataglobalstate_->GetDiscret();

  XFieldState()->ResetNonStandardDofs(*full_discret);
}
