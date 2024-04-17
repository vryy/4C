/*----------------------------------------------------------------------*/
/*! \file

\brief Proxy to a set of degrees of freedom

\level 1


*/
/*----------------------------------------------------------------------*/


#include "baci_lib_dofset_proxy.hpp"

#include "baci_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::DofSetProxy::DofSetProxy(DofSetInterface* dofset)
    : dofset_(dofset), isassigned_(dofset->Filled())
{
  dofset->Register(this);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::DofSetProxy::~DofSetProxy()
{
  if (dofset_ != nullptr) dofset_->Unregister(this);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetProxy::AddDofSettoList()
{
  // We do nothing here as a proxy does not show up in the dof set list.
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::DofSetProxy::AssignDegreesOfFreedom(
    const Discretization& dis, const unsigned dspos, const int start)
{
  // This method does nothing, because the DofSetProxy is not supposed to assign dofs itself.
  // Instead, the original dofset assigns dofs when FillComplete() is called on its discretization.
  // This invokes the call to AssignDegreesOfFreedom on the original dofset. In
  // AssignDegreesOfFreedom NotifyAssigned() is called. This calls NotifyAssigned() on all
  // registered proxies.
  NotifyAssigned();
  return start;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetProxy::NotifyAssigned()
{
  if (dofset_ == nullptr)
    dserror("dofset_ pointer is nullptr");
  else
    isassigned_ = dofset_->Filled();

  DRT::DofSetBase::NotifyAssigned();
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetProxy::Reset()
{
  isassigned_ = false;
  NotifyReset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetProxy::Disconnect(DofSetInterface* dofset)
{
  if (dofset == dofset_)
    dofset_ = nullptr;
  else
    dserror("cannot disconnect from non-connected DofSet");

  // clear my Teuchos::rcps.
  Reset();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::DofSetProxy::Filled() const
{
  if (dofset_) return dofset_->Filled();

  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DofSetProxy::CheckIsAssigned() const
{
  // checks in debug mode only
  dsassert(isassigned_,
      "AssignDegreesOfFreedom was not called on parent dofset of this proxy,\n"
      "and/or this proxy was not notified.");
  dsassert(dofset_ != nullptr, "dofset_ pointer is nullptr");

  return;
}

FOUR_C_NAMESPACE_CLOSE
