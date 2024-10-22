// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_crosslinker_node.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_mat_crosslinkermat.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

CrossLinking::CrosslinkerNodeType CrossLinking::CrosslinkerNodeType::instance_;



/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Communication::ParObject* CrossLinking::CrosslinkerNodeType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  std::vector<double> dummycoord(3, 999.0);
  auto* crosslinker = new CrossLinking::CrosslinkerNode(-1, dummycoord, -1);
  crosslinker->unpack(buffer);
  return crosslinker;
}

/*----------------------------------------------------------------------------*
 |  ctor (public)                                              eichinger 10/16|
 *----------------------------------------------------------------------------*/
CrossLinking::CrosslinkerNodeDataContainer::CrosslinkerNodeDataContainer() : numbond_(0)
{
  clbspots_.clear();
  std::pair<int, int> pair;
  pair.first = -1;
  pair.second = -1;
  clbspots_.push_back(pair);  // first binding spot of crosslinker
  clbspots_.push_back(pair);  // second binding spot of crosslinker

  return;
}

/*----------------------------------------------------------------------------*
 |  Pack data                                                        (public) |
 |                                                             eichinger 10/16|
 *----------------------------------------------------------------------------*/
void CrossLinking::CrosslinkerNodeDataContainer::pack(Core::Communication::PackBuffer& data) const
{
  // add numbond
  add_to_pack(data, numbond_);
  // add clbspots_
  add_to_pack(data, clbspots_);

  return;
}

/*----------------------------------------------------------------------------*
 |  Unpack data                                                      (public) |
 |                                                             eichinger 10/16|
 *----------------------------------------------------------------------------*/
void CrossLinking::CrosslinkerNodeDataContainer::unpack(Core::Communication::UnpackBuffer& buffer)
{
  // numbond
  extract_from_pack(buffer, numbond_);
  // clbspots_
  extract_from_pack(buffer, clbspots_);

  return;
}

/*----------------------------------------------------------------------------*
 *  ctor (public)                                              eichinger 10/16|
 *----------------------------------------------------------------------------*/
CrossLinking::CrosslinkerNode::CrosslinkerNode(
    int id, const std::vector<double>& coords, const int owner)
    : Core::Nodes::Node(id, coords, owner), mat_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------------*
 |  copy-ctor (public)                                         eichinger 10/16|
 *----------------------------------------------------------------------------*/
CrossLinking::CrosslinkerNode::CrosslinkerNode(const CrossLinking::CrosslinkerNode& old)
    : Core::Nodes::Node(old)
{
  FOUR_C_THROW(
      "Copy constructor of CrosslinkerNodeDataContainer needs to "
      "implemented first");
  return;
}

/*----------------------------------------------------------------------------*
 |  Deep copy this instance and return pointer to it (public)                 |
 |                                                             eichinger 10/16|
 *----------------------------------------------------------------------------*/
CrossLinking::CrosslinkerNode* CrossLinking::CrosslinkerNode::clone() const
{
  CrossLinking::CrosslinkerNode* newnode = new CrossLinking::CrosslinkerNode(*this);
  return newnode;
}

/*----------------------------------------------------------------------------*
 |  << operator                                                eichinger 10/16|
 *----------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const CrossLinking::CrosslinkerNode& crosslinker_node)
{
  crosslinker_node.print(os);
  return os;
}

/*----------------------------------------------------------------------------*
 |  print this CrossslinkerNode (public)                       eichinger 10/16|
 *----------------------------------------------------------------------------*/
void CrossLinking::CrosslinkerNode::print(std::ostream& os) const
{
  // Print id and coordinates
  os << "Crosslinker ";
  Core::Nodes::Node::print(os);


  return;
}

/*----------------------------------------------------------------------------*
 |  Pack data                                                        (public) |
 |                                                             eichinger 10/16|
 *----------------------------------------------------------------------------*/
void CrossLinking::CrosslinkerNode::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class Core::Nodes::Node
  Core::Nodes::Node::pack(data);

  // add material
  bool hasmat = (mat_ != Teuchos::null);
  add_to_pack(data, hasmat);
  if (hasmat) mat_->pack(data);

  return;
}

/*----------------------------------------------------------------------------*
 |  Unpack data                                                      (public) |
 |                                                             eichinger 10/16|
 *----------------------------------------------------------------------------*/
void CrossLinking::CrosslinkerNode::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class Core::Nodes::Node
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer basedata_buffer(basedata);
  Core::Nodes::Node::unpack(basedata_buffer);

  // mat
  bool hasmat = extract_int(buffer);
  if (hasmat)
  {
    std::vector<char> tmp;
    extract_from_pack(buffer, tmp);
    Core::Communication::UnpackBuffer buffer_tmp(tmp);
    Core::Communication::ParObject* o = Core::Communication::factory(buffer_tmp);
    Mat::CrosslinkerMat* mat = dynamic_cast<Mat::CrosslinkerMat*>(o);
    if (mat == nullptr) FOUR_C_THROW("failed to unpack material");
    // unpack material
    mat_ = Teuchos::RCP(mat);
  }
  else
  {
    mat_ = Teuchos::null;
  }

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
  return;
}

///*----------------------------------------------------------------------------*
// |  Initialize data container                                  eichinger 10/16|
// *----------------------------------------------------------------------------*/
// void CrossLinking::CrosslinkerNode::initialize_data_container()
//{
//  // only initialize if not yet done
//  if (cldata_ ==Teuchos::null)
//    cldata_ = Teuchos::rcp(new CrossLinking::CrosslinkerNodeDataContainer());
//
//  return;
//}

/*----------------------------------------------------------------------------*
 |  create material class (public)                             eichinger 10/16|
 *---------------------------------------------------------------------------*/
void CrossLinking::CrosslinkerNode::set_material(int const matnum)
{
  Teuchos::RCP<Mat::CrosslinkerMat> mat =
      Teuchos::rcp_dynamic_cast<Mat::CrosslinkerMat>(Mat::factory(matnum));
  if (mat == Teuchos::null) FOUR_C_THROW("Invalid material given to crosslinker node. \n");
  mat_ = mat;
}

/*----------------------------------------------------------------------------*
 |  create material class (public)                             eichinger 10/16|
 *---------------------------------------------------------------------------*/
void CrossLinking::CrosslinkerNode::set_material(Teuchos::RCP<Core::Mat::Material> material)
{
  Teuchos::RCP<Mat::CrosslinkerMat> mat = Teuchos::rcp_dynamic_cast<Mat::CrosslinkerMat>(material);
  if (mat == Teuchos::null) FOUR_C_THROW("Invalid material given to crosslinker node. \n");
  mat_ = mat;
}

///*----------------------------------------------------------------------------*
// |  Reset data container                                       eichinger 10/16|
// *----------------------------------------------------------------------------*/
// void CrossLinking::CrosslinkerNode::ResetDataContainer()
//{
//  // reset to Teuchos::null
//  cldata_  = Teuchos::null;
//
//  return;
//}

FOUR_C_NAMESPACE_CLOSE
