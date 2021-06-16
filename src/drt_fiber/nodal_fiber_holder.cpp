/*----------------------------------------------------------------------*/
/*! \file

\brief This file implements a class that holds different nodal fibers

\level 3


*/
/*----------------------------------------------------------------------*/

#include "nodal_fiber_holder.H"
#include "../drt_lib/drt_parobject.H"
#include "drt_fiber_node.H"

void DRT::FIBER::NodalFiberHolder::SetFiber(
    DRT::FIBER::FiberType type, const std::vector<LINALG::Matrix<3, 1>>& fiber)
{
  fibers_.insert(
      std::pair<DRT::FIBER::FiberType, const std::vector<LINALG::Matrix<3, 1>>>(type, fiber));
}

const std::vector<LINALG::Matrix<3, 1>>& DRT::FIBER::NodalFiberHolder::GetFiber(
    DRT::FIBER::FiberType type) const
{
  return fibers_.at(type);
}

std::vector<LINALG::Matrix<3, 1>>& DRT::FIBER::NodalFiberHolder::GetFiberMutual(
    DRT::FIBER::FiberType type)
{
  return fibers_.at(type);
}

void DRT::FIBER::NodalFiberHolder::SetAngle(AngleType type, const std::vector<double>& angle)
{
  angles_.insert(std::pair<DRT::FIBER::AngleType, const std::vector<double>>(type, angle));
}

const std::vector<double>& DRT::FIBER::NodalFiberHolder::GetAngle(DRT::FIBER::AngleType type) const
{
  return angles_.at(type);
}

unsigned int DRT::FIBER::NodalFiberHolder::FibersSize() const { return fibers_.size(); }

unsigned int DRT::FIBER::NodalFiberHolder::AnglesSize() const { return angles_.size(); }

bool DRT::FIBER::NodalFiberHolder::ContainsFiber(DRT::FIBER::FiberType type) const
{
  return fibers_.find(type) != fibers_.end();
}
