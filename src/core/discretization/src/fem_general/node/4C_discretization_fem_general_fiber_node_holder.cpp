/*----------------------------------------------------------------------*/
/*! \file

\brief This file implements a class that holds different nodal fibers

\level 3


*/
/*----------------------------------------------------------------------*/

#include "4C_discretization_fem_general_fiber_node_holder.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_discretization_fem_general_fiber_node.hpp"

FOUR_C_NAMESPACE_OPEN

void CORE::Nodes::NodalFiberHolder::set_coordinate_system_direction(
    CORE::Nodes::CoordinateSystemDirection type,
    const std::vector<CORE::LINALG::Matrix<3, 1>>& fiber)
{
  coordinate_system_directions_.insert(std::pair<CORE::Nodes::CoordinateSystemDirection,
      const std::vector<CORE::LINALG::Matrix<3, 1>>>(type, fiber));
}

const std::vector<CORE::LINALG::Matrix<3, 1>>&
CORE::Nodes::NodalFiberHolder::get_coordinate_system_direction(
    CORE::Nodes::CoordinateSystemDirection type) const
{
  return coordinate_system_directions_.at(type);
}

std::vector<CORE::LINALG::Matrix<3, 1>>&
CORE::Nodes::NodalFiberHolder::get_coordinate_system_direction_mutual(
    CORE::Nodes::CoordinateSystemDirection type)
{
  return coordinate_system_directions_.at(type);
}

void CORE::Nodes::NodalFiberHolder::AddFiber(const std::vector<CORE::LINALG::Matrix<3, 1>>& fiber)
{
  fibers_.emplace_back(fiber);
}

const std::vector<CORE::LINALG::Matrix<3, 1>>& CORE::Nodes::NodalFiberHolder::GetFiber(
    std::size_t fiberid) const
{
  return fibers_.at(fiberid);
}

std::vector<CORE::LINALG::Matrix<3, 1>>& CORE::Nodes::NodalFiberHolder::GetFiberMutual(
    std::size_t fiberid)
{
  return fibers_.at(fiberid);
}

void CORE::Nodes::NodalFiberHolder::SetAngle(AngleType type, const std::vector<double>& angle)
{
  angles_.insert(std::pair<CORE::Nodes::AngleType, const std::vector<double>>(type, angle));
}

const std::vector<double>& CORE::Nodes::NodalFiberHolder::GetAngle(
    CORE::Nodes::AngleType type) const
{
  return angles_.at(type);
}

std::size_t CORE::Nodes::NodalFiberHolder::FibersSize() const { return fibers_.size(); }

std::size_t CORE::Nodes::NodalFiberHolder::coordinate_system_size() const
{
  return coordinate_system_directions_.size();
}

std::size_t CORE::Nodes::NodalFiberHolder::AnglesSize() const { return angles_.size(); }

bool CORE::Nodes::NodalFiberHolder::contains_coordinate_system_direction(
    CORE::Nodes::CoordinateSystemDirection type) const
{
  return coordinate_system_directions_.find(type) != coordinate_system_directions_.end();
}

FOUR_C_NAMESPACE_CLOSE
