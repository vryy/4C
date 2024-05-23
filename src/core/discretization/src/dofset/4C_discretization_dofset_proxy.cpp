/*----------------------------------------------------------------------*/
/*! \file

\brief Proxy to a set of degrees of freedom

\level 1


*/
/*----------------------------------------------------------------------*/


#include "4C_discretization_dofset_proxy.hpp"

#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::Dofsets::DofSetProxy::DofSetProxy(DofSetInterface* dofset)
    : dofset_(dofset), isassigned_(dofset->Filled())
{
  dofset->Register(this);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::Dofsets::DofSetProxy::~DofSetProxy()
{
  if (dofset_ != nullptr) dofset_->Unregister(this);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Dofsets::DofSetProxy::AddDofSettoList()
{
  // We do nothing here as a proxy does not show up in the dof set list.
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int CORE::Dofsets::DofSetProxy::assign_degrees_of_freedom(
    const DRT::Discretization& dis, const unsigned dspos, const int start)
{
  // This method does nothing, because the DofSetProxy is not supposed to assign dofs itself.
  // Instead, the original dofset assigns dofs when FillComplete() is called on its discretization.
  // This invokes the call to assign_degrees_of_freedom on the original dofset. In
  // assign_degrees_of_freedom NotifyAssigned() is called. This calls NotifyAssigned() on all
  // registered proxies.
  NotifyAssigned();
  return start;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Dofsets::DofSetProxy::NotifyAssigned()
{
  if (dofset_ == nullptr)
    FOUR_C_THROW("dofset_ pointer is nullptr");
  else
    isassigned_ = dofset_->Filled();

  DofSetBase::NotifyAssigned();
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Dofsets::DofSetProxy::Reset()
{
  isassigned_ = false;
  NotifyReset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Dofsets::DofSetProxy::Disconnect(DofSetInterface* dofset)
{
  if (dofset == dofset_)
    dofset_ = nullptr;
  else
    FOUR_C_THROW("cannot disconnect from non-connected DofSet");

  // clear my Teuchos::rcps.
  Reset();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CORE::Dofsets::DofSetProxy::Filled() const
{
  if (dofset_) return dofset_->Filled();

  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Dofsets::DofSetProxy::CheckIsAssigned() const
{
  // checks in debug mode only
  FOUR_C_ASSERT(isassigned_,
      "assign_degrees_of_freedom was not called on parent dofset of this proxy,\n"
      "and/or this proxy was not notified.");
  FOUR_C_ASSERT(dofset_ != nullptr, "dofset_ pointer is nullptr");

  return;
}

FOUR_C_NAMESPACE_CLOSE
