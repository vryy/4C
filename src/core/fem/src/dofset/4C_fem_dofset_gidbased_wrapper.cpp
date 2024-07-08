/*----------------------------------------------------------------------*/
/*! \file

 \brief Implementation of a dofset using a GID based mapping

\level 2


*/
/*----------------------------------------------------------------------*/

#include "4C_fem_dofset_gidbased_wrapper.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_dofset_proxy.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::DOFSets::DofSetGIDBasedWrapper::DofSetGIDBasedWrapper(
    Teuchos::RCP<Core::FE::Discretization> sourcedis,
    Teuchos::RCP<Core::DOFSets::DofSetInterface> sourcedofset)
    : DofSetBase(),
      sourcedis_(sourcedis),
      sourcedofset_(sourcedofset),
      isassigned_(sourcedofset->filled())
{
  if (sourcedofset_ == Teuchos::null) FOUR_C_THROW("Source dof set is null pointer.");
  if (sourcedis_ == Teuchos::null) FOUR_C_THROW("Source discretization is null pointer.");

  sourcedofset_->register_proxy(this);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::DOFSets::DofSetGIDBasedWrapper::~DofSetGIDBasedWrapper()
{
  if (sourcedofset_ != Teuchos::null) sourcedofset_->unregister(this);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetGIDBasedWrapper::reset()
{
  isassigned_ = false;
  notify_reset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Core::DOFSets::DofSetGIDBasedWrapper::assign_degrees_of_freedom(
    const Core::FE::Discretization& dis, const unsigned dspos, const int start)
{
  notify_assigned();
  return start;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetGIDBasedWrapper::notify_assigned()
{
  if (sourcedis_->node_col_map() == nullptr) FOUR_C_THROW("No NodeColMap on sourcedis");
  if (sourcedis_->element_col_map() == nullptr) FOUR_C_THROW("No ElementColMap on sourcedis");

  isassigned_ = sourcedofset_->filled();

  // call base class
  DofSetBase::notify_assigned();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetGIDBasedWrapper::disconnect(DofSetInterface* dofset)
{
  if (dofset == sourcedofset_.get())
  {
    sourcedofset_ = Teuchos::null;
    sourcedis_ = Teuchos::null;
  }
  else
    FOUR_C_THROW("cannot disconnect from non-connected DofSet");

  // clear my Teuchos::rcps.
  reset();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::DOFSets::DofSetGIDBasedWrapper::check_is_assigned() const
{
  // checks in debug mode only
  FOUR_C_ASSERT(isassigned_,
      "assign_degrees_of_freedom was not called on parent dofset of this proxy,\n"
      "and/or this proxy was not notified.");
  FOUR_C_ASSERT(sourcedofset_ != Teuchos::null, "dofset_ pointer is nullptr");

  return;
}

FOUR_C_NAMESPACE_CLOSE
