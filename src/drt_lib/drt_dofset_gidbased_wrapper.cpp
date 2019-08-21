/*----------------------------------------------------------------------*/
/*! \file

 \brief Implementation of a dofset using a GID based mapping

\level 2

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#include "drt_dofset_gidbased_wrapper.H"

#include "drt_dofset_proxy.H"
#include "drt_discret.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::DofSetGIDBasedWrapper::DofSetGIDBasedWrapper(
    Teuchos::RCP<DRT::Discretization> sourcedis, Teuchos::RCP<DRT::DofSetInterface> sourcedofset)
    : DofSetBase(),
      sourcedis_(sourcedis),
      sourcedofset_(sourcedofset),
      isassigned_(sourcedofset->Filled())
{
  if (sourcedofset_ == Teuchos::null) dserror("Source dof set is null pointer.");
  if (sourcedis_ == Teuchos::null) dserror("Source discretization is null pointer.");

  sourcedofset_->Register(this);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::DofSetGIDBasedWrapper::~DofSetGIDBasedWrapper()
{
  if (sourcedofset_ != Teuchos::null) sourcedofset_->Unregister(this);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetGIDBasedWrapper::Reset()
{
  isassigned_ = false;
  NotifyReset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSetGIDBasedWrapper::AssignDegreesOfFreedom(
    const Discretization& dis, const unsigned dspos, const int start)
{
  NotifyAssigned();
  return start;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetGIDBasedWrapper::NotifyAssigned()
{
  if (sourcedis_->NodeColMap() == NULL) dserror("No NodeColMap on sourcedis");
  if (sourcedis_->ElementColMap() == NULL) dserror("No ElementColMap on sourcedis");

  isassigned_ = sourcedofset_->Filled();

  // call base class
  DRT::DofSetBase::NotifyAssigned();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetGIDBasedWrapper::Disconnect(DofSetInterface* dofset)
{
  if (dofset == sourcedofset_.get())
  {
    sourcedofset_ = Teuchos::null;
    sourcedis_ = Teuchos::null;
  }
  else
    dserror("cannot disconnect from non-connected DofSet");

  // clear my Teuchos::rcps.
  Reset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetGIDBasedWrapper::CheckIsAssigned() const
{
  // checks in debug mode only
  dsassert(isassigned_,
      "AssignDegreesOfFreedom was not called on parent dofset of this proxy,\n"
      "and/or this proxy was not notified.");
  dsassert(sourcedofset_ != Teuchos::null, "dofset_ pointer is NULL");

  return;
}
