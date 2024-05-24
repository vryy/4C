/*----------------------------------------------------------------------*/
/*! \file

 \brief Implementation of a dofset using a GID based mapping

\level 2


*/
/*----------------------------------------------------------------------*/

#include "4C_discretization_dofset_gidbased_wrapper.hpp"

#include "4C_discretization_dofset_proxy.hpp"
#include "4C_lib_discret.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::Dofsets::DofSetGIDBasedWrapper::DofSetGIDBasedWrapper(
    Teuchos::RCP<DRT::Discretization> sourcedis,
    Teuchos::RCP<CORE::Dofsets::DofSetInterface> sourcedofset)
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
CORE::Dofsets::DofSetGIDBasedWrapper::~DofSetGIDBasedWrapper()
{
  if (sourcedofset_ != Teuchos::null) sourcedofset_->Unregister(this);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Dofsets::DofSetGIDBasedWrapper::Reset()
{
  isassigned_ = false;
  NotifyReset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int CORE::Dofsets::DofSetGIDBasedWrapper::assign_degrees_of_freedom(
    const DRT::Discretization& dis, const unsigned dspos, const int start)
{
  NotifyAssigned();
  return start;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Dofsets::DofSetGIDBasedWrapper::NotifyAssigned()
{
  if (sourcedis_->NodeColMap() == nullptr) FOUR_C_THROW("No NodeColMap on sourcedis");
  if (sourcedis_->ElementColMap() == nullptr) FOUR_C_THROW("No ElementColMap on sourcedis");

  isassigned_ = sourcedofset_->Filled();

  // call base class
  DofSetBase::NotifyAssigned();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::Dofsets::DofSetGIDBasedWrapper::Disconnect(DofSetInterface* dofset)
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
void CORE::Dofsets::DofSetGIDBasedWrapper::check_is_assigned() const
{
  // checks in debug mode only
  FOUR_C_ASSERT(isassigned_,
      "assign_degrees_of_freedom was not called on parent dofset of this proxy,\n"
      "and/or this proxy was not notified.");
  FOUR_C_ASSERT(sourcedofset_ != Teuchos::null, "dofset_ pointer is nullptr");

  return;
}

FOUR_C_NAMESPACE_CLOSE
