/*----------------------------------------------------------------------*/
/*! \file

\brief transparent independent dofset

\level 2


*/
/*----------------------------------------------------------------------*/


#include "4C_discretization_dofset_transparent_independent.hpp"

#include "4C_discretization_dofset.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN



CORE::Dofsets::TransparentIndependentDofSet::TransparentIndependentDofSet(
    Teuchos::RCP<DRT::Discretization> sourcedis, bool parallel)
    : TransparentDofSet(sourcedis, parallel)
{
  return;
}

int CORE::Dofsets::TransparentIndependentDofSet::AssignDegreesOfFreedom(
    const DRT::Discretization& dis, const unsigned dspos, const int start)
{
  // first, we call the standard AssignDegreesOfFreedom from the base class
  int count = IndependentDofSet::AssignDegreesOfFreedom(dis, dspos, start);

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

int CORE::Dofsets::TransparentIndependentDofSet::NumDofPerNode(const DRT::Node& node) const
{
  return DofSet::NumDofPerNode(node);
}

FOUR_C_NAMESPACE_CLOSE
