// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_general_immersed_node.hpp"

#include "4C_comm_pack_helpers.hpp"

FOUR_C_NAMESPACE_OPEN


Core::Nodes::ImmersedNodeType Core::Nodes::ImmersedNodeType::instance_;


/*----------------------------------------------------------------------*
 |  kind of ctor (public)                                   rauch 11/14 |
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Core::Nodes::ImmersedNodeType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  std::vector<double> dummycoord(3, 999.0);
  Node* object = new Core::Nodes::ImmersedNode(-1, dummycoord, -1);
  object->unpack(buffer);
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
Core::Nodes::ImmersedNode* Core::Nodes::ImmersedNode::clone() const
{
  Core::Nodes::ImmersedNode* newnode = new Core::Nodes::ImmersedNode(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 |  << operator                                             rauch 11/14 |
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const Core::Nodes::ImmersedNode& immersednode)
{
  immersednode.print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                             rauch 11/14 |
 *----------------------------------------------------------------------*/
void Core::Nodes::ImmersedNode::print(std::ostream& os) const
{
  os << "Immersed ";
  Node::print(os);

  if (is_boundary_immersed())
    os << " Immersed Boundary  ";
  else
    os << " NOT Immersed Boundary ";

  if (is_matched())
    os << " Matched ";
  else
    os << " NOT Matched ";

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          rauch 11/14 |
 *----------------------------------------------------------------------*/
void Core::Nodes::ImmersedNode::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // add base class Core::Nodes::Node
  Node::pack(data);

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
void Core::Nodes::ImmersedNode::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class Core::Nodes::Node
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Communication::UnpackBuffer basedata_buffer(basedata);
  Node::unpack(basedata_buffer);

  // isimersedboundary_
  IsBoundaryImmersed_ = extract_int(buffer);
  // ismatched_
  ismatched_ = extract_int(buffer);

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
  return;
}


/*----------------------------------------------------------------------*
 |  Visualization Data                                         (public) |
 |                                                          rauch 03/17 |
 *----------------------------------------------------------------------*/
void Core::Nodes::ImmersedNode::vis_names(std::map<std::string, int>& names)
{
  names.insert(std::pair<std::string, int>("IsBoundaryImmersedNode", 1));
  return;
}


/*----------------------------------------------------------------------*
 |  Query data to be visualized by BINIO                       (public) |
 |                                                          rauch 03/17 |
 *----------------------------------------------------------------------*/
bool Core::Nodes::ImmersedNode::vis_data(const std::string& name, std::vector<double>& data)
{
  if (name == "IsBoundaryImmersedNode")
  {
    if ((int)data.size() < 1) FOUR_C_THROW("Size mismatch");
    data[0] = is_boundary_immersed();
    return true;
  }
  return false;
}

FOUR_C_NAMESPACE_CLOSE
