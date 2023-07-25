/*----------------------------------------------------------------------*/
/*! \file

\brief transparent independent dofset

\level 2


*/
/*----------------------------------------------------------------------*/


#include "baci_lib_dofset_transparent_independent.H"
#include "baci_lib_dofset.H"

#include "baci_linalg_utils_sparse_algebra_math.H"



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
