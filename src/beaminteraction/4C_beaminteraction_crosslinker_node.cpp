/*-----------------------------------------------------------*/
/*! \file

\brief A class for a crosslinker node


\date Oct, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_beaminteraction_crosslinker_node.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_mat_crosslinkermat.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

CrossLinking::CrosslinkerNodeType CrossLinking::CrosslinkerNodeType::instance_;



/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Communication::ParObject* CrossLinking::CrosslinkerNodeType::Create(
    const std::vector<char>& data)
{
  std::vector<double> dummycoord(3, 999.0);
  auto* crosslinker = new CrossLinking::CrosslinkerNode(-1, dummycoord, -1);
  crosslinker->unpack(data);
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
  Core::Communication::ParObject::add_to_pack(data, numbond_);
  // add clbspots_
  Core::Communication::ParObject::add_to_pack(data, clbspots_);

  return;
}

/*----------------------------------------------------------------------------*
 |  Unpack data                                                      (public) |
 |                                                             eichinger 10/16|
 *----------------------------------------------------------------------------*/
void CrossLinking::CrosslinkerNodeDataContainer::unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  // numbond
  Core::Communication::ParObject::extract_from_pack(position, data, numbond_);
  // clbspots_
  Core::Communication::ParObject::extract_from_pack(position, data, clbspots_);

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
CrossLinking::CrosslinkerNode* CrossLinking::CrosslinkerNode::Clone() const
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
  int type = UniqueParObjectId();
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
void CrossLinking::CrosslinkerNode::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Core::Nodes::Node
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Core::Nodes::Node::unpack(basedata);

  // mat
  bool hasmat = extract_int(position, data);
  if (hasmat)
  {
    std::vector<char> tmp;
    extract_from_pack(position, data, tmp);
    Core::Communication::ParObject* o = Core::Communication::Factory(tmp);
    Mat::CrosslinkerMat* mat = dynamic_cast<Mat::CrosslinkerMat*>(o);
    if (mat == nullptr) FOUR_C_THROW("failed to unpack material");
    // unpack material
    mat_ = Teuchos::rcp(mat);
  }
  else
  {
    mat_ = Teuchos::null;
  }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
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
void CrossLinking::CrosslinkerNode::SetMaterial(int const matnum)
{
  Teuchos::RCP<Mat::CrosslinkerMat> mat =
      Teuchos::rcp_dynamic_cast<Mat::CrosslinkerMat>(Mat::Factory(matnum));
  if (mat == Teuchos::null) FOUR_C_THROW("Invalid material given to crosslinker node. \n");
  mat_ = mat;
}

/*----------------------------------------------------------------------------*
 |  create material class (public)                             eichinger 10/16|
 *---------------------------------------------------------------------------*/
void CrossLinking::CrosslinkerNode::SetMaterial(Teuchos::RCP<Core::Mat::Material> material)
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
