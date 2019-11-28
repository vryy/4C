/*----------------------------------------------------------------------*/
/*! \file

\brief This file implements a class that holds different nodal fibers

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#include "nodal_fiber_holder.H"

void DRT::FIBER::NodalFiberHolder::SetFiber(
    DRT::FIBER::FiberType type, const std::vector<LINALG::Matrix<3, 1>>& fiber)
{
  map_.insert(
      std::pair<DRT::FIBER::FiberType, const std::vector<LINALG::Matrix<3, 1>>>(type, fiber));
}

const std::vector<LINALG::Matrix<3, 1>>& DRT::FIBER::NodalFiberHolder::GetFiber(
    DRT::FIBER::FiberType type) const
{
  return map_.at(type);
}

unsigned int DRT::FIBER::NodalFiberHolder::Size() const { return map_.size(); }
