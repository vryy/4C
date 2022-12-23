/*----------------------------------------------------------------------*/
/*! \file

\brief This file implements a class that holds different nodal fibers

\level 3


*/
/*----------------------------------------------------------------------*/

#include "nodal_fiber_holder.H"
#include "parobject.H"
#include "fiber_node.H"

void DRT::FIBER::NodalFiberHolder::SetCoordinateSystemDirection(
    DRT::FIBER::CoordinateSystemDirection type, const std::vector<LINALG::Matrix<3, 1>>& fiber)
{
  coordinateSystemDirections_.insert(
      std::pair<DRT::FIBER::CoordinateSystemDirection, const std::vector<LINALG::Matrix<3, 1>>>(
          type, fiber));
}

const std::vector<LINALG::Matrix<3, 1>>& DRT::FIBER::NodalFiberHolder::GetCoordinateSystemDirection(
    DRT::FIBER::CoordinateSystemDirection type) const
{
  return coordinateSystemDirections_.at(type);
}

std::vector<LINALG::Matrix<3, 1>>& DRT::FIBER::NodalFiberHolder::GetCoordinateSystemDirectionMutual(
    DRT::FIBER::CoordinateSystemDirection type)
{
  return coordinateSystemDirections_.at(type);
}

void DRT::FIBER::NodalFiberHolder::AddFiber(const std::vector<LINALG::Matrix<3, 1>>& fiber)
{
  fibers_.emplace_back(fiber);
}

const std::vector<LINALG::Matrix<3, 1>>& DRT::FIBER::NodalFiberHolder::GetFiber(
    std::size_t fiberid) const
{
  return fibers_.at(fiberid);
}

std::vector<LINALG::Matrix<3, 1>>& DRT::FIBER::NodalFiberHolder::GetFiberMutual(std::size_t fiberid)
{
  return fibers_.at(fiberid);
}

void DRT::FIBER::NodalFiberHolder::SetAngle(AngleType type, const std::vector<double>& angle)
{
  angles_.insert(std::pair<DRT::FIBER::AngleType, const std::vector<double>>(type, angle));
}

const std::vector<double>& DRT::FIBER::NodalFiberHolder::GetAngle(DRT::FIBER::AngleType type) const
{
  return angles_.at(type);
}

std::size_t DRT::FIBER::NodalFiberHolder::FibersSize() const { return fibers_.size(); }

std::size_t DRT::FIBER::NodalFiberHolder::CoordinateSystemSize() const
{
  return coordinateSystemDirections_.size();
}

std::size_t DRT::FIBER::NodalFiberHolder::AnglesSize() const { return angles_.size(); }

bool DRT::FIBER::NodalFiberHolder::ContainsCoordinateSystemDirection(
    DRT::FIBER::CoordinateSystemDirection type) const
{
  return coordinateSystemDirections_.find(type) != coordinateSystemDirections_.end();
}
