/*----------------------------------------------------------------------*/
/*! \file

\brief Proxy to a set of degrees of freedom

\level 1


*/
/*----------------------------------------------------------------------*/


#include "4C_fem_dofset_proxy.hpp"

#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::DOFSets::DofSetProxy::DofSetProxy(DofSetInterface* dofset)
    : dofset_(dofset), isassigned_(dofset->filled())
{
  dofset->Register(this);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::DOFSets::DofSetProxy::~DofSetProxy()
{
  if (dofset_ != nullptr) dofset_->unregister(this);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetProxy::add_dof_setto_list()
{
  // We do nothing here as a proxy does not show up in the dof set list.
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Core::DOFSets::DofSetProxy::assign_degrees_of_freedom(
    const Core::FE::Discretization& dis, const unsigned dspos, const int start)
{
  // This method does nothing, because the DofSetProxy is not supposed to assign dofs itself.
  // Instead, the original dofset assigns dofs when fill_complete() is called on its discretization.
  // This invokes the call to assign_degrees_of_freedom on the original dofset. In
  // assign_degrees_of_freedom NotifyAssigned() is called. This calls NotifyAssigned() on all
  // registered proxies.
  notify_assigned();
  return start;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetProxy::notify_assigned()
{
  if (dofset_ == nullptr)
    FOUR_C_THROW("dofset_ pointer is nullptr");
  else
    isassigned_ = dofset_->filled();

  DofSetBase::notify_assigned();
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetProxy::reset()
{
  isassigned_ = false;
  notify_reset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetProxy::disconnect(DofSetInterface* dofset)
{
  if (dofset == dofset_)
    dofset_ = nullptr;
  else
    FOUR_C_THROW("cannot disconnect from non-connected DofSet");

  // clear my Teuchos::rcps.
  reset();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Core::DOFSets::DofSetProxy::filled() const
{
  if (dofset_) return dofset_->filled();

  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetProxy::check_is_assigned() const
{
  // checks in debug mode only
  FOUR_C_ASSERT(isassigned_,
      "assign_degrees_of_freedom was not called on parent dofset of this proxy,\n"
      "and/or this proxy was not notified.");
  FOUR_C_ASSERT(dofset_ != nullptr, "dofset_ pointer is nullptr");

  return;
}

FOUR_C_NAMESPACE_CLOSE
