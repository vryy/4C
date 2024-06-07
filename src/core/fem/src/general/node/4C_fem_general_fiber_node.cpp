/*----------------------------------------------------------------------*/
/*! \file

\brief   This is basically a (3d-) node with an additional fiber direction.

\level 2
*----------------------------------------------------------------------*/

#include "4C_fem_general_fiber_node.hpp"

#include "4C_fem_general_fiber_node_holder.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN

Core::Nodes::FiberNodeType Core::Nodes::FiberNodeType::instance_;


Core::Communication::ParObject* Core::Nodes::FiberNodeType::Create(const std::vector<char>& data)
{
  std::vector<double> dummy_coords(3, 999.0);
  std::map<CoordinateSystemDirection, std::array<double, 3>> coordinateSystemDirections;
  std::vector<std::array<double, 3>> fibers;
  std::map<AngleType, double> angles;
  auto* object = new FiberNode(-1, dummy_coords, coordinateSystemDirections, fibers, angles, -1);
  object->Unpack(data);
  return object;
}

Core::Nodes::FiberNode::FiberNode(int id, const std::vector<double>& coords,
    std::map<CoordinateSystemDirection, std::array<double, 3>> coordinateSystemDirections,
    std::vector<std::array<double, 3>> fibers, std::map<AngleType, double> angles, const int owner)
    : Core::Nodes::Node(id, coords, owner),
      coordinateSystemDirections_(std::move(coordinateSystemDirections)),
      fibers_(std::move(fibers)),
      angles_(std::move(angles))
{
}

/*
  Deep copy the derived class and return pointer to it
*/
Core::Nodes::FiberNode* Core::Nodes::FiberNode::Clone() const
{
  auto* newfn = new Core::Nodes::FiberNode(*this);

  return newfn;
}

/*
  Pack this class so it can be communicated

  Pack and Unpack are used to communicate this fiber node

*/
void Core::Nodes::FiberNode::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  Core::Nodes::Node::AddtoPack(data, type);
  // add base class of fiber node
  Core::Nodes::Node::Pack(data);

  // Add fiber data
  Core::Communication::ParObject::AddtoPack(data, fibers_);
  Core::Communication::ParObject::AddtoPack(data, coordinateSystemDirections_);
  Core::Communication::ParObject::AddtoPack(data, angles_);
}

/*
  Unpack data from a char vector into this class

  Pack and Unpack are used to communicate this fiber node
*/
void Core::Nodes::FiberNode::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Node
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Core::Nodes::Node::Unpack(basedata);

  // extract fiber data
  Core::Communication::ParObject::ExtractfromPack(position, data, fibers_);
  Core::Communication::ParObject::ExtractfromPack(position, data, coordinateSystemDirections_);
  Core::Communication::ParObject::ExtractfromPack(position, data, angles_);
}

/*
  Print this fiber node
*/
void Core::Nodes::FiberNode::Print(std::ostream& os) const
{
  os << "Fiber Node :";
  Core::Nodes::Node::Print(os);
  os << "(" << fibers_.size() << " fibers, " << angles_.size() << " angles)";
}

FOUR_C_NAMESPACE_CLOSE
