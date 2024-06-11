/*----------------------------------------------------------------------*/
/*! \file

\brief specialized node for immersed problems.

\level 2

*----------------------------------------------------------------------*/

#include "4C_fem_general_immersed_node.hpp"

FOUR_C_NAMESPACE_OPEN


Core::Nodes::ImmersedNodeType Core::Nodes::ImmersedNodeType::instance_;


/*----------------------------------------------------------------------*
 |  kind of ctor (public)                                   rauch 11/14 |
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Core::Nodes::ImmersedNodeType::Create(const std::vector<char>& data)
{
  std::vector<double> dummycoord(3, 999.0);
  Node* object = new Core::Nodes::ImmersedNode(-1, dummycoord, -1);
  object->Unpack(data);
  return object;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           rauch 11/14 |
 *----------------------------------------------------------------------*/
Core::Nodes::ImmersedNode::ImmersedNode(
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
Core::Nodes::ImmersedNode::ImmersedNode(const Core::Nodes::ImmersedNode& old)
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
Core::Nodes::ImmersedNode* Core::Nodes::ImmersedNode::Clone() const
{
  Core::Nodes::ImmersedNode* newnode = new Core::Nodes::ImmersedNode(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 |  << operator                                             rauch 11/14 |
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const Core::Nodes::ImmersedNode& immersednode)
{
  immersednode.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                             rauch 11/14 |
 *----------------------------------------------------------------------*/
void Core::Nodes::ImmersedNode::Print(std::ostream& os) const
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
void Core::Nodes::ImmersedNode::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // add base class Core::Nodes::Node
  Node::Pack(data);

  // add IsBoundaryImmersed_
  add_to_pack(data, IsBoundaryImmersed_);
  // add ismatched_
  add_to_pack(data, ismatched_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          rauch 11/14 |
 *----------------------------------------------------------------------*/
void Core::Nodes::ImmersedNode::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Core::Nodes::Node
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
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
void Core::Nodes::ImmersedNode::VisNames(std::map<std::string, int>& names)
{
  names.insert(std::pair<std::string, int>("IsBoundaryImmersedNode", 1));
  return;
}


/*----------------------------------------------------------------------*
 |  Query data to be visualized by BINIO                       (public) |
 |                                                          rauch 03/17 |
 *----------------------------------------------------------------------*/
bool Core::Nodes::ImmersedNode::VisData(const std::string& name, std::vector<double>& data)
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
