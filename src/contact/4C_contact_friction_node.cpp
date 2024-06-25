/*-----------------------------------------------------------------------*/
/*! \file
\brief A class for a frictional contact node
\level 2
*/
/*-----------------------------------------------------------------------*/

#include "4C_contact_friction_node.hpp"

#include "4C_contact_defines.hpp"
#include "4C_contact_element.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

CONTACT::FriNodeType CONTACT::FriNodeType::instance_;

Core::Communication::ParObject* CONTACT::FriNodeType::Create(const std::vector<char>& data)
{
  std::vector<double> x(3, 0.0);
  std::vector<int> dofs(0);

  // TODO: friplus = true for all nodes!!! change this with pack/unpack
  CONTACT::FriNode* node = new CONTACT::FriNode(0, x, 0, dofs, false, false, true);
  node->unpack(data);

  return node;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             mgit 01/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNodeDataContainer::FriNodeDataContainer()
    : slip_(false), slipold_(false), drowsold_(0), drowsoldLTL_(0)
{
  for (int i = 0; i < 3; ++i)
  {
    jump()[i] = 0.0;
    traction()[i] = 0.0;
    tractionold()[i] = 0.0;
    tractionLTL()[i] = 0.0;
    tractionoldLTL()[i] = 0.0;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNodeDataContainer::pack(Core::Communication::PackBuffer& data) const
{
  // add jump_
  Core::Communication::ParObject::add_to_pack(data, jump_, 3 * sizeof(double));
  // add slip_
  Core::Communication::ParObject::add_to_pack(data, slip_);
  // add slip_
  Core::Communication::ParObject::add_to_pack(data, slipold_);
  // add traction_
  Core::Communication::ParObject::add_to_pack(data, traction_, 3 * sizeof(double));
  // add tractionold_
  Core::Communication::ParObject::add_to_pack(data, tractionold_, 3 * sizeof(double));
  // add traction_
  Core::Communication::ParObject::add_to_pack(data, tractionLTL_, 3 * sizeof(double));
  // add tractionold_
  Core::Communication::ParObject::add_to_pack(data, tractionoldLTL_, 3 * sizeof(double));

  // add drowsold_,mrowsold_,mnodesold_
  int hasdata = drowsold_.size();

  Core::Communication::ParObject::add_to_pack(data, hasdata);

  if (hasdata != 0)
  {
    int dentries = (int)drowsold_.size();
    Core::Communication::ParObject::add_to_pack(data, dentries);
    Core::Communication::ParObject::add_to_pack(data, drowsold_);
    Core::Communication::ParObject::add_to_pack(data, mrowsold_);
    Core::Communication::ParObject::add_to_pack(data, mnodesold_);
  }

  int hasdata2 = drowsoldLTL_.size();
  Core::Communication::ParObject::add_to_pack(data, hasdata2);
  if (hasdata2 != 0)
  {
    int dentries = (int)drowsoldLTL_.size();
    Core::Communication::ParObject::add_to_pack(data, dentries);
    Core::Communication::ParObject::add_to_pack(data, drowsoldLTL_);
    Core::Communication::ParObject::add_to_pack(data, mrowsoldLTL_);
  }
  // add derivjump
  int hasdataderivjump = derivjump_.size();
  Core::Communication::ParObject::add_to_pack(data, hasdataderivjump);

  if (hasdataderivjump != 0)
  {
    for (int i = 0; i < hasdataderivjump; i++)
      Core::Communication::ParObject::add_to_pack(data, (derivjump_[i]));
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNodeDataContainer::unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  // jump_
  Core::Communication::ParObject::extract_from_pack(position, data, jump_, 3 * sizeof(double));
  // slip_
  slip_ = Core::Communication::ParObject::extract_int(position, data);
  // slipold_
  slipold_ = Core::Communication::ParObject::extract_int(position, data);
  // traction_
  Core::Communication::ParObject::extract_from_pack(position, data, traction_, 3 * sizeof(double));
  // tractionold_
  Core::Communication::ParObject::extract_from_pack(
      position, data, tractionold_, 3 * sizeof(double));
  // traction_
  Core::Communication::ParObject::extract_from_pack(
      position, data, tractionLTL_, 3 * sizeof(double));
  // tractionold_
  Core::Communication::ParObject::extract_from_pack(
      position, data, tractionoldLTL_, 3 * sizeof(double));

  // drowsold_,mrowsold_,mnodesold_
  int hasdata;
  Core::Communication::ParObject::extract_from_pack(position, data, hasdata);

  if (hasdata != 0)
  {
    int dentries = Core::Communication::ParObject::extract_int(position, data);

    drowsold_.resize(dentries);
    Core::Communication::ParObject::extract_from_pack(position, data, drowsold_);
    Core::Communication::ParObject::extract_from_pack(position, data, mrowsold_);
    Core::Communication::ParObject::extract_from_pack(position, data, mnodesold_);
  }

  // drowsold_,mrowsold_,mnodesold_
  int hasdata2;
  Core::Communication::ParObject::extract_from_pack(position, data, hasdata2);

  if (hasdata2 != 0)
  {
    int dentries = Core::Communication::ParObject::extract_int(position, data);

    drowsoldLTL_.resize(dentries);
    Core::Communication::ParObject::extract_from_pack(position, data, drowsoldLTL_);
    Core::Communication::ParObject::extract_from_pack(position, data, mrowsoldLTL_);
  }

  // and derivjump_
  int hasdataderivjump;
  Core::Communication::ParObject::extract_from_pack(position, data, hasdataderivjump);

  if (hasdataderivjump != 0)
  {
    derivjump_.resize(hasdataderivjump);
    for (int i = 0; i < hasdataderivjump; i++)
    {
      Core::Communication::ParObject::extract_from_pack(position, data, derivjump_[i]);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             mgit 07/11|
 *----------------------------------------------------------------------*/
CONTACT::FriNodeWearDataContainer::FriNodeWearDataContainer()
    : weightedwear_(0.0), deltaweightedwear_(0.0)
{
  wcurr_[0] = 0.0;
  wold_[0] = 0.0;
  waccu_[0] = 0.0;
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           farah 10/13|
 *----------------------------------------------------------------------*/
void CONTACT::FriNodeWearDataContainer::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::ParObject::add_to_pack(data, weightedwear_);
  Core::Communication::ParObject::add_to_pack(data, deltaweightedwear_);

  // add d2row
  int hasdata = d2rows_.size();
  Core::Communication::ParObject::add_to_pack(data, hasdata);

  if (hasdata != 0)
  {
    for (int i = 0; i < hasdata; i++)
    {
      Core::Communication::ParObject::add_to_pack(data, (d2rows_[i]));
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           farah 10/13|
 *----------------------------------------------------------------------*/
void CONTACT::FriNodeWearDataContainer::unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  Core::Communication::ParObject::extract_from_pack(position, data, weightedwear_);
  Core::Communication::ParObject::extract_from_pack(position, data, deltaweightedwear_);

  // d2rows_
  int hasdata;
  Core::Communication::ParObject::extract_from_pack(position, data, hasdata);

  if (hasdata != 0)
  {
    d2rows_.resize(hasdata);
    for (int i = 0; i < hasdata; i++)
    {
      Core::Communication::ParObject::extract_from_pack(position, data, d2rows_[i]);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             mgit 02/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNode::FriNode(int id, const std::vector<double>& coords, const int owner,
    const std::vector<int>& dofs, const bool isslave, const bool initactive, const bool friplus)
    : CONTACT::Node(id, coords, owner, dofs, isslave, initactive), wear_(friplus)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        mgit 02/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNode::FriNode(const CONTACT::FriNode& old) : CONTACT::Node(old), wear_(old.wear_)
{
  // not yet used and thus not necessarily consistent
  FOUR_C_THROW("FriNode copy-ctor not yet implemented");

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of FriNode and return pointer to it (public)|
 |                                                            mgit 02/10|
 *----------------------------------------------------------------------*/
CONTACT::FriNode* CONTACT::FriNode::Clone() const
{
  CONTACT::FriNode* newnode = new CONTACT::FriNode(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 |  << operator                                               mgit 02/10|
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const CONTACT::FriNode& frinode)
{
  frinode.print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                               mgit 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::print(std::ostream& os) const
{
  // Print id and coordinates
  os << "Contact ";
  CONTACT::Node::print(os);
  if (IsSlave())
    if (IsInitActive()) os << " InitActive ";

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // add base class Mortar::Node
  CONTACT::Node::pack(data);

  // add data_
  bool hasdata = (fridata_ != Teuchos::null);
  add_to_pack(data, hasdata);
  if (hasdata) fridata_->pack(data);

  bool hasweardata = (weardata_ != Teuchos::null);
  add_to_pack(data, hasweardata);
  if (hasweardata) weardata_->pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class CONTACT::Node
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  CONTACT::Node::unpack(basedata);

  // **************************
  // FriData
  bool hasdata = extract_int(position, data);
  if (hasdata)
  {
    fridata_ = Teuchos::rcp(new CONTACT::FriNodeDataContainer());
    fridata_->unpack(position, data);
  }
  else
    fridata_ = Teuchos::null;

  // **************************
  // FriDataPlus
  bool hasdataplus = extract_int(position, data);
  if (hasdataplus)
  {
    weardata_ = Teuchos::rcp(new CONTACT::FriNodeWearDataContainer());
    weardata_->unpack(position, data);
  }
  else
    weardata_ = Teuchos::null;

  // Check
  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 | calculate the apparent coefficient of friction            seitz 11/15|
 *----------------------------------------------------------------------*/
double CONTACT::FriNode::FrCoeff(const double& frcoeff_in)
{
  // return the friction coefficient, if we do not have a TSI problem
  if (cTSIdata_ == Teuchos::null) return frcoeff_in;

  // in TSI case, the friction coefficient is temperature dependent
  else
  {
    double maxT = std::max(TSIData().Temp(), TSIData().TempMaster());
    return frcoeff_in * (maxT - cTSIdata_->Temp_Dam()) * (maxT - cTSIdata_->Temp_Dam()) /
           ((cTSIdata_->Temp_Dam() - cTSIdata_->Temp_Ref()) *
               (cTSIdata_->Temp_Dam() - cTSIdata_->Temp_Ref()));
  }
}

/*----------------------------------------------------------------------*
 | calculate derivative of apparent coefficient of friction  seitz 11/15|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::derivFrCoeffTemp(
    const double& frcoeff_in, std::map<int, double>& derivT, std::map<int, double>& derivDisp)
{
  derivT.clear();
  derivDisp.clear();

  // if we do not have a TSI problem, the friction coefficient is constant
  if (cTSIdata_ == Teuchos::null) return;

  double T_dam = cTSIdata_->Temp_Dam();
  double T_ref = cTSIdata_->Temp_Ref();
  if (cTSIdata_->Temp() > cTSIdata_->TempMaster())
  {
    double maxT = TSIData().Temp();
    derivT[Dofs()[0]] += 2. * frcoeff_in * (maxT - T_dam) / ((T_dam - T_ref) * (T_dam - T_ref));
    derivDisp.clear();
  }
  else
  {
    double maxT = TSIData().TempMaster();
    for (std::map<int, double>::const_iterator i = TSIData().DerivTempMasterTemp().begin();
         i != TSIData().DerivTempMasterTemp().end(); ++i)
      derivT[i->first] +=
          2. * frcoeff_in * (maxT - T_dam) / ((T_dam - T_ref) * (T_dam - T_ref)) * i->second;
    for (std::map<int, double>::const_iterator i = TSIData().DerivTempMasterDisp().begin();
         i != TSIData().DerivTempMasterDisp().end(); ++i)
      derivDisp[i->first] +=
          2. * frcoeff_in * (maxT - T_dam) / ((T_dam - T_ref) * (T_dam - T_ref)) * i->second;
  }
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'SNodes' set                        gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddSNode(int node)
{
  FriData().GetSNodes().insert(node);

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'MNodes' set                        gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddMNode(int node)
{
  FriData().GetMNodes().insert(node);

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'D2' map                              farah 06/13|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddD2Value(int& row, int& col, double& val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave() == true) FOUR_C_THROW("AddD2Value: function called for slave node %i", Id());

  // check if this has been called before
  if ((int)WearData().GetD2().size() == 0) WearData().GetD2().resize(NumDof());

  // check row index input
  if ((int)WearData().GetD2().size() <= row)
    FOUR_C_THROW("AddD2Value: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  std::map<int, double>& d2map = WearData().GetD2()[row];
  d2map[col] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'DerivJump' map                     gitterle 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddDerivJumpValue(int& row, const int& col, double val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave() == false) FOUR_C_THROW("AddJumpValue: function called for master node %i", Id());
  if (IsOnBound() == true) FOUR_C_THROW("AddJumpValue: function called for boundary node %i", Id());

  // check if this has been called before
  if ((int)FriData().GetDerivJump().size() == 0) FriData().GetDerivJump().resize(NumDof());

  // check row index input
  if ((int)FriData().GetDerivJump().size() <= row)
    FOUR_C_THROW("AddDerivJumpValue: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  std::map<int, double>& zmap = FriData().GetDerivJump()[row];
  zmap[col] += val;

  return;
}

/*----------------------------------------------------------------------*
|  Add jump value from gp-wise integration of slip          farah 08/13|
*----------------------------------------------------------------------*/
void CONTACT::FriNode::AddJumpValue(double val, int k)
{
  // check if this is a master node or slave boundary node
  if (IsSlave() == false) FOUR_C_THROW("AddJumpValue: function called for master node %i", Id());
  if (IsOnBound() == true) FOUR_C_THROW("AddJumpValue: function called for boundary node %i", Id());

  FriData().jump_var()[k] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'T' map                               farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddTValue(int& row, int& col, double& val)
{
  // check if this is a master node or slave boundary node
  //  if (IsSlave()==false)
  //    FOUR_C_THROW("AddTValue: function called for master node %i", Id());

  // check if this has been called before
  if ((int)WearData().GetT().size() == 0) WearData().GetT().resize(NumDof());

  // check row index input
  if ((int)WearData().GetT().size() <= row)
    FOUR_C_THROW("AddTValue: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  std::map<int, double>& tmap = WearData().GetT()[row];
  tmap[col] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'E' map                               farah 09/13|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::AddEValue(int& row, int& col, double& val)
{
  // check if this is a master node or slave boundary node
  //  if (IsSlave()==false)
  //    FOUR_C_THROW("AddEValue: function called for master node %i", Id());

  // check if this has been called before
  if ((int)WearData().GetE().size() == 0) WearData().GetE().resize(NumDof());

  // check row index input
  if ((int)WearData().GetE().size() <= row)
    FOUR_C_THROW("AddEValue: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  std::map<int, double>& emap = WearData().GetE()[row];
  emap[col] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Store nodal entries of D and M to old ones             gitterle 12/08|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::StoreDMOld()
{
  // clear and zero nodal vectors
  FriData().GetDOld().clear();
  FriData().GetMOld().clear();
  FriData().GetDOldLTL().clear();
  FriData().GetMOldLTL().clear();

  // write drows_ to drowsold_
  FriData().GetDOld() = MoData().GetD();
  FriData().GetMOld() = MoData().GetM();

  // write drows_ to drowsold_
  FriData().GetDOldLTL() = MoData().GetDltl();
  FriData().GetMOldLTL() = MoData().GetMltl();

  // also vectors containing the according master nodes
  FriData().GetMNodesOld().clear();
  FriData().GetMNodesOld() = FriData().GetMNodes();

  return;
}

/*----------------------------------------------------------------------*
 |  Store nodal entries penalty tractions to old ones      gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::StoreTracOld()
{
  // write entries to old ones
  for (int j = 0; j < 3; ++j) FriData().tractionold()[j] = FriData().traction()[j];

  for (int j = 0; j < 3; ++j) FriData().tractionoldLTL()[j] = FriData().tractionLTL()[j];

  return;
}

/*-----------------------------------------------------------------------*
 |  Set the value of deltawear                             gitterle 12/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::add_delta_weighted_wear_value(double& val)
{
  // add given value to deltawear_
  WearData().DeltaWeightedWear() += val;
}

/*----------------------------------------------------------------------*
 |  Initialize data container                             gitterle 10/09|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::initialize_data_container()
{
  // get maximum size of lin vectors
  linsize_ = 0;
  for (int i = 0; i < NumElement(); ++i)
    for (int j = 0; j < Elements()[i]->num_node(); ++j)
      linsize_ += Elements()[i]->NumDofPerNode(*(Elements()[i]->Nodes()[j]));

  // get maximum size of lin vectors
  dentries_ = 0;
  for (int i = 0; i < NumElement(); ++i)
    for (int j = 0; j < Elements()[i]->num_node(); ++j)
      dentries_ += Elements()[i]->NumDofPerNode(*(Elements()[i]->Nodes()[j]));

  // only initialize if not yet done
  if (modata_ == Teuchos::null && codata_ == Teuchos::null && fridata_ == Teuchos::null)
  {
    modata_ = Teuchos::rcp(new Mortar::NodeDataContainer());
    codata_ = Teuchos::rcp(new CONTACT::NodeDataContainer());
    fridata_ = Teuchos::rcp(new CONTACT::FriNodeDataContainer());
  }

  // initialize data container for wear and tsi problems
  if (wear_ == true)
  {
    if (weardata_ == Teuchos::null)
      weardata_ = Teuchos::rcp(new CONTACT::FriNodeWearDataContainer());
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Reset data container                                      popp 09/10|
 *----------------------------------------------------------------------*/
void CONTACT::FriNode::ResetDataContainer()
{
  // reset to Teuchos::null
  fridata_ = Teuchos::null;
  weardata_ = Teuchos::null;
  codata_ = Teuchos::null;
  modata_ = Teuchos::null;

  return;
}

FOUR_C_NAMESPACE_CLOSE
