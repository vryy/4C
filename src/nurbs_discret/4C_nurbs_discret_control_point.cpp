/*----------------------------------------------------------------------*/
/*! \file

\brief   This is basically a (3d-) node with an additional weight.
   The weight is required for the evaluation of the nurbs
   basis functions.

   note that X() is not the coordinate of some grid point
   anymore, it's just the control point position


\level 2

*----------------------------------------------------------------------*/

#include "4C_nurbs_discret_control_point.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::Nurbs::ControlPointType Discret::Nurbs::ControlPointType::instance_;


Core::Communication::ParObject* Discret::Nurbs::ControlPointType::Create(
    const std::vector<char>& data)
{
  std::vector<double> dummycoord(3, 999.0);
  double dummyweight = 999.;
  Discret::Nurbs::ControlPoint* object =
      new Discret::Nurbs::ControlPoint(-1, dummycoord, dummyweight, -1);
  object->Unpack(data);
  return object;
}

/*
  Standard ctor
 */
Discret::Nurbs::ControlPoint::ControlPoint(
    int id, const std::vector<double>& coords, const double weight, const int owner)
    : Core::Nodes::Node(id, coords, owner), w_(weight)
{
  return;
}

/*
  Copy Constructor

  Makes a deep copy of a control point

*/
Discret::Nurbs::ControlPoint::ControlPoint(const Discret::Nurbs::ControlPoint& old)
    : Core::Nodes::Node(old), w_(old.W())
{
  return;
}

/*
  Deep copy the derived class and return pointer to it

*/
Discret::Nurbs::ControlPoint* Discret::Nurbs::ControlPoint::Clone() const
{
  Discret::Nurbs::ControlPoint* newcp = new Discret::Nurbs::ControlPoint(*this);

  return newcp;
}


/*
  Pack this class so it can be communicated

  Pack and Unpack are used to communicate this control point

*/
void Discret::Nurbs::ControlPoint::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  Core::Nodes::Node::add_to_pack(data, type);
  // add base class of control point
  Core::Nodes::Node::Pack(data);
  // add weight
  Core::Nodes::Node::add_to_pack(data, &w_, sizeof(double));

  return;
}

/*
  Unpack data from a char vector into this class

  Pack and Unpack are used to communicate this control point
*/
void Discret::Nurbs::ControlPoint::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Node
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Core::Nodes::Node::Unpack(basedata);
  // extract weight
  Core::Nodes::Node::extract_from_pack(position, data, w_);

  return;
}

/*
  Print this control point
*/
void Discret::Nurbs::ControlPoint::Print(std::ostream& os) const
{
  os << "Control Point :";
  Core::Nodes::Node::Print(os);
  os << "\n+ additional weight ";
  os << w_ << "\n";
  return;
}

FOUR_C_NAMESPACE_CLOSE
