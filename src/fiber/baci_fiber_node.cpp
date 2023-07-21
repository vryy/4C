/*----------------------------------------------------------------------*/
/*! \file

\brief   This is basically a (3d-) node with an additional fiber direction.

\level 2
*----------------------------------------------------------------------*/

#include <utility>

#include "baci_fiber_node.H"
#include "baci_fiber_nodal_fiber_holder.H"

DRT::FIBER::FiberNodeType DRT::FIBER::FiberNodeType::instance_;


DRT::ParObject* DRT::FIBER::FiberNodeType::Create(const std::vector<char>& data)
{
  std::array<double, 3> dummy_coords = {999., 999., 999.};
  std::map<FIBER::CoordinateSystemDirection, std::array<double, 3>> coordinateSystemDirections;
  std::vector<std::array<double, 3>> fibers;
  std::map<FIBER::AngleType, double> angles;
  auto* object =
      new DRT::FIBER::FiberNode(-1, dummy_coords, coordinateSystemDirections, fibers, angles, -1);
  object->Unpack(data);
  return object;
}

DRT::FIBER::FiberNode::FiberNode(int id, std::array<double, 3> coords,
    std::map<FIBER::CoordinateSystemDirection, std::array<double, 3>> coordinateSystemDirections,
    std::vector<std::array<double, 3>> fibers, std::map<FIBER::AngleType, double> angles,
    const int owner)
    : DRT::Node(id, coords.data(), owner),
      coordinateSystemDirections_(std::move(coordinateSystemDirections)),
      fibers_(std::move(fibers)),
      angles_(std::move(angles))
{
}

/*
  Deep copy the derived class and return pointer to it
*/
DRT::FIBER::FiberNode* DRT::FIBER::FiberNode::Clone() const
{
  auto* newfn = new DRT::FIBER::FiberNode(*this);

  return newfn;
}

/*
  Pack this class so it can be communicated

  Pack and Unpack are used to communicate this fiber node

*/
void DRT::FIBER::FiberNode::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  DRT::Node::AddtoPack(data, type);
  // add base class of fiber node
  DRT::Node::Pack(data);

  // Add fiber data
  DRT::ParObject::AddtoPack(data, fibers_);
  DRT::ParObject::AddtoPack(data, coordinateSystemDirections_);
  DRT::ParObject::AddtoPack(data, angles_);
}

/*
  Unpack data from a char vector into this class

  Pack and Unpack are used to communicate this fiber node
*/
void DRT::FIBER::FiberNode::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class Node
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::Node::Unpack(basedata);

  // extract fiber data
  DRT::ParObject::ExtractfromPack(position, data, fibers_);
  DRT::ParObject::ExtractfromPack(position, data, coordinateSystemDirections_);
  DRT::ParObject::ExtractfromPack(position, data, angles_);
}

/*
  Print this fiber node
*/
void DRT::FIBER::FiberNode::Print(std::ostream& os) const
{
  os << "Fiber Node :";
  DRT::Node::Print(os);
  os << "(" << fibers_.size() << " fibers, " << angles_.size() << " angles)";
}
