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

CROSSLINKING::CrosslinkerNodeType CROSSLINKING::CrosslinkerNodeType::instance_;



/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CORE::COMM::ParObject* CROSSLINKING::CrosslinkerNodeType::Create(const std::vector<char>& data)
{
  std::vector<double> dummycoord(3, 999.0);
  auto* crosslinker = new CROSSLINKING::CrosslinkerNode(-1, dummycoord, -1);
  crosslinker->Unpack(data);
  return crosslinker;
}

/*----------------------------------------------------------------------------*
 |  ctor (public)                                              eichinger 10/16|
 *----------------------------------------------------------------------------*/
CROSSLINKING::CrosslinkerNodeDataContainer::CrosslinkerNodeDataContainer() : numbond_(0)
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
void CROSSLINKING::CrosslinkerNodeDataContainer::Pack(CORE::COMM::PackBuffer& data) const
{
  // add numbond
  CORE::COMM::ParObject::AddtoPack(data, numbond_);
  // add clbspots_
  CORE::COMM::ParObject::AddtoPack(data, clbspots_);

  return;
}

/*----------------------------------------------------------------------------*
 |  Unpack data                                                      (public) |
 |                                                             eichinger 10/16|
 *----------------------------------------------------------------------------*/
void CROSSLINKING::CrosslinkerNodeDataContainer::Unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  // numbond
  CORE::COMM::ParObject::ExtractfromPack(position, data, numbond_);
  // clbspots_
  CORE::COMM::ParObject::ExtractfromPack(position, data, clbspots_);

  return;
}

/*----------------------------------------------------------------------------*
 *  ctor (public)                                              eichinger 10/16|
 *----------------------------------------------------------------------------*/
CROSSLINKING::CrosslinkerNode::CrosslinkerNode(
    int id, const std::vector<double>& coords, const int owner)
    : CORE::Nodes::Node(id, coords, owner), mat_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------------*
 |  copy-ctor (public)                                         eichinger 10/16|
 *----------------------------------------------------------------------------*/
CROSSLINKING::CrosslinkerNode::CrosslinkerNode(const CROSSLINKING::CrosslinkerNode& old)
    : CORE::Nodes::Node(old)
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
CROSSLINKING::CrosslinkerNode* CROSSLINKING::CrosslinkerNode::Clone() const
{
  CROSSLINKING::CrosslinkerNode* newnode = new CROSSLINKING::CrosslinkerNode(*this);
  return newnode;
}

/*----------------------------------------------------------------------------*
 |  << operator                                                eichinger 10/16|
 *----------------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const CROSSLINKING::CrosslinkerNode& crosslinker_node)
{
  crosslinker_node.Print(os);
  return os;
}

/*----------------------------------------------------------------------------*
 |  print this CrossslinkerNode (public)                       eichinger 10/16|
 *----------------------------------------------------------------------------*/
void CROSSLINKING::CrosslinkerNode::Print(std::ostream& os) const
{
  // Print id and coordinates
  os << "Crosslinker ";
  CORE::Nodes::Node::Print(os);


  return;
}

/*----------------------------------------------------------------------------*
 |  Pack data                                                        (public) |
 |                                                             eichinger 10/16|
 *----------------------------------------------------------------------------*/
void CROSSLINKING::CrosslinkerNode::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class CORE::Nodes::Node
  CORE::Nodes::Node::Pack(data);

  // add material
  bool hasmat = (mat_ != Teuchos::null);
  AddtoPack(data, hasmat);
  if (hasmat) mat_->Pack(data);

  return;
}

/*----------------------------------------------------------------------------*
 |  Unpack data                                                      (public) |
 |                                                             eichinger 10/16|
 *----------------------------------------------------------------------------*/
void CROSSLINKING::CrosslinkerNode::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class CORE::Nodes::Node
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  CORE::Nodes::Node::Unpack(basedata);

  // mat
  bool hasmat = ExtractInt(position, data);
  if (hasmat)
  {
    std::vector<char> tmp;
    ExtractfromPack(position, data, tmp);
    CORE::COMM::ParObject* o = CORE::COMM::Factory(tmp);
    MAT::CrosslinkerMat* mat = dynamic_cast<MAT::CrosslinkerMat*>(o);
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
// void CROSSLINKING::CrosslinkerNode::initialize_data_container()
//{
//  // only initialize if not yet done
//  if (cldata_ ==Teuchos::null)
//    cldata_ = Teuchos::rcp(new CROSSLINKING::CrosslinkerNodeDataContainer());
//
//  return;
//}

/*----------------------------------------------------------------------------*
 |  create material class (public)                             eichinger 10/16|
 *---------------------------------------------------------------------------*/
void CROSSLINKING::CrosslinkerNode::SetMaterial(int const matnum)
{
  Teuchos::RCP<MAT::CrosslinkerMat> mat =
      Teuchos::rcp_dynamic_cast<MAT::CrosslinkerMat>(MAT::Factory(matnum));
  if (mat == Teuchos::null) FOUR_C_THROW("Invalid material given to crosslinker node. \n");
  mat_ = mat;
}

/*----------------------------------------------------------------------------*
 |  create material class (public)                             eichinger 10/16|
 *---------------------------------------------------------------------------*/
void CROSSLINKING::CrosslinkerNode::SetMaterial(Teuchos::RCP<CORE::MAT::Material> material)
{
  Teuchos::RCP<MAT::CrosslinkerMat> mat = Teuchos::rcp_dynamic_cast<MAT::CrosslinkerMat>(material);
  if (mat == Teuchos::null) FOUR_C_THROW("Invalid material given to crosslinker node. \n");
  mat_ = mat;
}

///*----------------------------------------------------------------------------*
// |  Reset data container                                       eichinger 10/16|
// *----------------------------------------------------------------------------*/
// void CROSSLINKING::CrosslinkerNode::ResetDataContainer()
//{
//  // reset to Teuchos::null
//  cldata_  = Teuchos::null;
//
//  return;
//}

FOUR_C_NAMESPACE_CLOSE
