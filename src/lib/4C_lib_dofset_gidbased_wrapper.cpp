/*----------------------------------------------------------------------*/
/*! \file

 \brief Implementation of a dofset using a GID based mapping

\level 2


*/
/*----------------------------------------------------------------------*/

#include "4C_lib_dofset_gidbased_wrapper.hpp"

#include "4C_lib_discret.hpp"
#include "4C_lib_dofset_proxy.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::DofSetGIDBasedWrapper::DofSetGIDBasedWrapper(
    Teuchos::RCP<DRT::Discretization> sourcedis, Teuchos::RCP<DRT::DofSetInterface> sourcedofset)
    : DofSetBase(),
      sourcedis_(sourcedis),
      sourcedofset_(sourcedofset),
      isassigned_(sourcedofset->Filled())
{
  if (sourcedofset_ == Teuchos::null) FOUR_C_THROW("Source dof set is null pointer.");
  if (sourcedis_ == Teuchos::null) FOUR_C_THROW("Source discretization is null pointer.");

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
  if (sourcedis_->NodeColMap() == nullptr) FOUR_C_THROW("No NodeColMap on sourcedis");
  if (sourcedis_->ElementColMap() == nullptr) FOUR_C_THROW("No ElementColMap on sourcedis");

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
    FOUR_C_THROW("cannot disconnect from non-connected DofSet");

  // clear my Teuchos::rcps.
  Reset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetGIDBasedWrapper::CheckIsAssigned() const
{
  // checks in debug mode only
  FOUR_C_ASSERT(isassigned_,
      "AssignDegreesOfFreedom was not called on parent dofset of this proxy,\n"
      "and/or this proxy was not notified.");
  FOUR_C_ASSERT(sourcedofset_ != Teuchos::null, "dofset_ pointer is nullptr");

  return;
}

FOUR_C_NAMESPACE_CLOSE
