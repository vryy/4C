/*-----------------------------------------------------------*/
/*! \file

\brief A class for a crosslinker node


\date Oct, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#include "baci_beaminteraction_crosslinker_node.H"

#include "baci_comm_utils_factory.H"
#include "baci_mat_crosslinkermat.H"
#include "baci_utils_exceptions.H"

BACI_NAMESPACE_OPEN

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
    : DRT::Node(id, coords, owner), mat_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------------*
 |  copy-ctor (public)                                         eichinger 10/16|
 *----------------------------------------------------------------------------*/
CROSSLINKING::CrosslinkerNode::CrosslinkerNode(const CROSSLINKING::CrosslinkerNode& old)
    : DRT::Node(old)
{
  dserror(
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
  DRT::Node::Print(os);


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
  // add base class DRT::Node
  DRT::Node::Pack(data);

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

  // extract base class DRT::Node
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::Node::Unpack(basedata);

  // mat
  bool hasmat = ExtractInt(position, data);
  if (hasmat)
  {
    std::vector<char> tmp;
    ExtractfromPack(position, data, tmp);
    CORE::COMM::ParObject* o = CORE::COMM::Factory(tmp);
    MAT::CrosslinkerMat* mat = dynamic_cast<MAT::CrosslinkerMat*>(o);
    if (mat == nullptr) dserror("failed to unpack material");
    // unpack material
    mat_ = Teuchos::rcp(mat);
  }
  else
  {
    mat_ = Teuchos::null;
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

///*----------------------------------------------------------------------------*
// |  Initialize data container                                  eichinger 10/16|
// *----------------------------------------------------------------------------*/
// void CROSSLINKING::CrosslinkerNode::InitializeDataContainer()
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
      Teuchos::rcp_dynamic_cast<MAT::CrosslinkerMat>(MAT::Material::Factory(matnum));
  if (mat == Teuchos::null) dserror("Invalid material given to crosslinker node. \n");
  mat_ = mat;
}

/*----------------------------------------------------------------------------*
 |  create material class (public)                             eichinger 10/16|
 *---------------------------------------------------------------------------*/
void CROSSLINKING::CrosslinkerNode::SetMaterial(Teuchos::RCP<MAT::Material> material)
{
  Teuchos::RCP<MAT::CrosslinkerMat> mat = Teuchos::rcp_dynamic_cast<MAT::CrosslinkerMat>(material);
  if (mat == Teuchos::null) dserror("Invalid material given to crosslinker node. \n");
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

BACI_NAMESPACE_CLOSE
