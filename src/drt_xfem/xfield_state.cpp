/*---------------------------------------------------------------------------*/
/*! \file

\brief state management for one XFEM discretization

\level 3

\maintainer Matthias Mayr

*/
/*---------------------------------------------------------------------------*/


#include "xfield_state.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_discret_xfem.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::XFieldState::XFieldState()
    : isinit_(false),
      issetup_(false),
      wizard_(Teuchos::null),
      condition_manager_(Teuchos::null),
      xdofset_(Teuchos::null),
      xfield_discret_ptr_(Teuchos::null),
      field_discret_ptr_(Teuchos::null)
{
  // intentionally left blank
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::XFieldState::Init(const Teuchos::RCP<XFEM::ConditionManager>& condition_manager,
    const Teuchos::RCP<GEO::CutWizard>& wizard, const Teuchos::RCP<XFEM::XFEMDofSet>& xdofset,
    const Teuchos::RCP<DRT::DiscretizationInterface>& xfielddiscret,
    const Teuchos::RCP<DRT::DiscretizationInterface>& fielddiscret)
{
  // Ensure, that the Setup() routines are called afterwards.
  issetup_ = false;

  condition_manager_ = condition_manager;
  wizard_ = wizard;

  xdofset_ = xdofset;

  // store a pointer to the field discretization
  field_discret_ptr_ = fielddiscret;

  // store a pointer to the xfield discretization
  xfield_discret_ptr_ = xfielddiscret;

  isinit_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::XFieldState::SetNewState(const XFEM::XFieldState& xstate)
{
  this->isinit_ = xstate.isinit_;
  this->issetup_ = xstate.issetup_;

  this->condition_manager_ = xstate.condition_manager_;
  this->wizard_ = xstate.wizard_;

  this->xdofset_ = xstate.xdofset_;

  this->field_discret_ptr_ = xstate.field_discret_ptr_;
  this->xfield_discret_ptr_ = xstate.xfield_discret_ptr_;
}
