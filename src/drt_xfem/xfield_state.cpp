/*---------------------------------------------------------------------------*/
/*!
\file xfield_state.cpp

\brief state management for one XFEM discretization

\level 3

\maintainer Michael hiermeier

\date Jun 27, 2016
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
void XFEM::XFieldState::Init(
    const Teuchos::RCP<XFEM::ConditionManager> &     condition_manager,
    const Teuchos::RCP<GEO::CutWizard> &             wizard,
    const Teuchos::RCP<XFEM::XFEMDofSet> &           xdofset,
    const Teuchos::RCP<DRT::DiscretizationInterface> & xfielddiscret,
    const Teuchos::RCP<DRT::DiscretizationInterface> & fielddiscret)
{
  /* If the member variables have been already created, we have to destroy
   * them before we can reinitialize them. */
  if (isinit_ and issetup_)
    Destroy();

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


