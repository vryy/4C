/*----------------------------------------------------------------------*/
/*! \file

\brief specialized node for immersed problems.

\level 2

*----------------------------------------------------------------------*/

#include "4C_discretization_fem_general_immersed_node.hpp"

FOUR_C_NAMESPACE_OPEN


CORE::Nodes::ImmersedNodeType CORE::Nodes::ImmersedNodeType::instance_;


/*----------------------------------------------------------------------*
 |  kind of ctor (public)                                   rauch 11/14 |
 *----------------------------------------------------------------------*/
CORE::COMM::ParObject* CORE::Nodes::ImmersedNodeType::Create(const std::vector<char>& data)
{
  std::vector<double> dummycoord(3, 999.0);
  Node* object = new CORE::Nodes::ImmersedNode(-1, dummycoord, -1);
  object->Unpack(data);
  return object;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           rauch 11/14 |
 *----------------------------------------------------------------------*/
CORE::Nodes::ImmersedNode::ImmersedNode(
    const int id, const std::vector<double>& coords, const int owner)
    : Node(id, coords, owner),
      ismatched_(false),
      IsBoundaryImmersed_(false),
      IsPseudoBoundary_(false)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      rauch 11/14 |
 *----------------------------------------------------------------------*/
CORE::Nodes::ImmersedNode::ImmersedNode(const CORE::Nodes::ImmersedNode& old)
    : Node(old),
      ismatched_(old.ismatched_),
      IsBoundaryImmersed_(old.IsBoundaryImmersed_),
      IsPseudoBoundary_(old.IsPseudoBoundary_)
{
  FOUR_C_THROW("ERROR: ImmersedNode copy-ctor has not been used before. Be careful when using it.");
  return;
}



/*----------------------------------------------------------------------*
 |  Deep copy this instance of Node and return pointer to it (public)   |
 |                                                          rauch 11/14 |
 *----------------------------------------------------------------------*/
CORE::Nodes::ImmersedNode* CORE::Nodes::ImmersedNode::Clone() const
{
  CORE::Nodes::ImmersedNode* newnode = new CORE::Nodes::ImmersedNode(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 |  << operator                                             rauch 11/14 |
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const CORE::Nodes::ImmersedNode& immersednode)
{
  immersednode.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                             rauch 11/14 |
 *----------------------------------------------------------------------*/
void CORE::Nodes::ImmersedNode::Print(std::ostream& os) const
{
  os << "Immersed ";
  Node::Print(os);

  if (IsBoundaryImmersed())
    os << " Immersed Boundary  ";
  else
    os << " NOT Immersed Boundary ";

  if (IsMatched())
    os << " Matched ";
  else
    os << " NOT Matched ";

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          rauch 11/14 |
 *----------------------------------------------------------------------*/
void CORE::Nodes::ImmersedNode::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // add base class CORE::Nodes::Node
  Node::Pack(data);

  // add IsBoundaryImmersed_
  AddtoPack(data, IsBoundaryImmersed_);
  // add ismatched_
  AddtoPack(data, ismatched_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          rauch 11/14 |
 *----------------------------------------------------------------------*/
void CORE::Nodes::ImmersedNode::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class CORE::Nodes::Node
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Node::Unpack(basedata);

  // isimersedboundary_
  IsBoundaryImmersed_ = ExtractInt(position, data);
  // ismatched_
  ismatched_ = ExtractInt(position, data);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}


/*----------------------------------------------------------------------*
 |  Visualization Data                                         (public) |
 |                                                          rauch 03/17 |
 *----------------------------------------------------------------------*/
void CORE::Nodes::ImmersedNode::VisNames(std::map<std::string, int>& names)
{
  names.insert(std::pair<std::string, int>("IsBoundaryImmersedNode", 1));
  return;
}


/*----------------------------------------------------------------------*
 |  Query data to be visualized by BINIO                       (public) |
 |                                                          rauch 03/17 |
 *----------------------------------------------------------------------*/
bool CORE::Nodes::ImmersedNode::VisData(const std::string& name, std::vector<double>& data)
{
  if (name == "IsBoundaryImmersedNode")
  {
    if ((int)data.size() < 1) FOUR_C_THROW("Size mismatch");
    data[0] = IsBoundaryImmersed();
    return true;
  }
  return false;
}

FOUR_C_NAMESPACE_CLOSE
