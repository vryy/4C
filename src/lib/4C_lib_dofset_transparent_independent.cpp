/*----------------------------------------------------------------------*/
/*! \file

\brief transparent independent dofset

\level 2


*/
/*----------------------------------------------------------------------*/


#include "4C_lib_dofset_transparent_independent.hpp"

#include "4C_lib_dofset.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN



DRT::TransparentIndependentDofSet::TransparentIndependentDofSet(
    Teuchos::RCP<DRT::Discretization> sourcedis, bool parallel)
    : TransparentDofSet(sourcedis, parallel)
{
  return;
}

int DRT::TransparentIndependentDofSet::AssignDegreesOfFreedom(
    const DRT::Discretization& dis, const unsigned dspos, const int start)
{
  // first, we call the standard AssignDegreesOfFreedom from the base class
  int count = DRT::IndependentDofSet::AssignDegreesOfFreedom(dis, dspos, start);

  if (!parallel_)
  {
    TransferDegreesOfFreedom(*sourcedis_, dis, start);
  }
  else
  {
    ParallelTransferDegreesOfFreedom(*sourcedis_, dis, start);
  }

  // tell all proxies (again!)
  NotifyAssigned();

  return count;
}

int DRT::TransparentIndependentDofSet::NumDofPerNode(const DRT::Node& node) const
{
  return DRT::DofSet::NumDofPerNode(node);
}

FOUR_C_NAMESPACE_CLOSE
