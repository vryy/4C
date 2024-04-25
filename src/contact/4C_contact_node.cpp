/*----------------------------------------------------------------------*/
/*! \file
\brief A class for a contact node

\level 2

*/
/*----------------------------------------------------------------------*/

#include "4C_contact_node.hpp"

#include "4C_contact_aug_contact_integrator_utils.hpp"
#include "4C_contact_defines.hpp"
#include "4C_contact_element.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

CONTACT::NodeType CONTACT::NodeType::instance_;


CORE::COMM::ParObject* CONTACT::NodeType::Create(const std::vector<char>& data)
{
  std::vector<double> x(3, 0.0);
  std::vector<int> dofs(0);
  auto* node = new CONTACT::Node(0, x, 0, dofs, false, false);
  node->Unpack(data);
  return node;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mgit 02/10 |
 *----------------------------------------------------------------------*/
CONTACT::NodeDataContainer::NodeDataContainer()
    : grow_(1.0e12),
      gnts_(1.0e12),
      glts_(1.0e12),
      activeold_(false),
      derivn_(0, 0),     // init deriv normal to length 0 with 0 entries per direction
      derivtxi_(0, 0),   // init deriv txi    to length 0 with 0 entries per direction
      derivteta_(0, 0),  // init deriv teta   to length 0 with 0 entries per direction
      alpha_(0),
      kappa_(1.0)
{
  // set all tangent entries to 0.0
  for (int i = 0; i < 3; ++i)
  {
    txi()[i] = 0.0;
    teta()[i] = 0.0;
    Getgltl()[i] = 1.0e12;
    Getjumpltl()[i] = 1.0e12;
  }
  GetDerivGltl().resize(3);
  GetDerivJumpltl().resize(3);
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::NodeDataContainer::Pack(CORE::COMM::PackBuffer& data) const
{
  // add txi_
  CORE::COMM::ParObject::AddtoPack(data, txi_, 3 * sizeof(double));
  // add teta_
  CORE::COMM::ParObject::AddtoPack(data, teta_, 3 * sizeof(double));
  // add grow_
  CORE::COMM::ParObject::AddtoPack(data, grow_);
  // add kappa_
  CORE::COMM::ParObject::AddtoPack(data, kappa_);
  // add activeold_
  CORE::COMM::ParObject::AddtoPack(data, activeold_);
  // add n_old_
  CORE::COMM::ParObject::AddtoPack(data, n_old_, 3 * sizeof(double));

  // no need to pack derivs_
  // (these will evaluated anew anyway)

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::NodeDataContainer::Unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  // txi_
  CORE::COMM::ParObject::ExtractfromPack(position, data, txi_, 3 * sizeof(double));
  // teta_
  CORE::COMM::ParObject::ExtractfromPack(position, data, teta_, 3 * sizeof(double));
  // grow_
  CORE::COMM::ParObject::ExtractfromPack(position, data, grow_);
  // kappa_
  CORE::COMM::ParObject::ExtractfromPack(position, data, kappa_);
  // activeold_
  activeold_ = CORE::COMM::ParObject::ExtractInt(position, data);
  // n_old_
  CORE::COMM::ParObject::ExtractfromPack(position, data, n_old_, 3 * sizeof(double));

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::NodeDataContainer::NodeDataContainer(Node& parentNode)
    : mentries_(-1),
      kappa_(1.0e12),
      w_gap_(1.0e12),
      aug_a_(1.0e12),
      d_aug_a_(0),
      parent_node_(parentNode)
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::AUG::NodeDataContainer::NodeDataContainer(
    Node& parentNode, const int slMaElementAreaRatio, const bool isTriangleOnMaster)
    : kappa_(1.0e12), w_gap_(1.0e12), aug_a_(1.0e12), parent_node_(parentNode)
{
  mentries_ = ApproximateMEntries(slMaElementAreaRatio, isTriangleOnMaster);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int CONTACT::AUG::NodeDataContainer::ApproximateMEntries(
    const int slMaElementAreaRatio, const bool isTriangleOnMaster) const
{
  // number of adjacent slave elements
  // The max function catches the rare case of a pure ghost node, i.e. the node
  // belongs to a different proc than all surrounding elements/nodes on this proc
  const int numSlEle = std::max(parent_node_.NumElement(), 1);

  // numSlEle slave elements project approximately into numMaEle if there is no
  // off-set
  int numMaEle = numSlEle * slMaElementAreaRatio;

  /* These numMaEle elements are approximately surrounded by numSurround
   * elements + numCorner elements */
  int numMaSurrounding = -1;
  int numCorner = -1;
  if (isTriangleOnMaster)
  {
    numMaSurrounding = 4 * std::ceil(std::sqrt(static_cast<double>(numMaEle) / 2.0));
    numCorner = 2;
  }
  else
  {
    numMaSurrounding = 2 * std::ceil(std::sqrt(static_cast<double>(numMaEle)));
    numCorner = 1;
  }

  // approximate number of relevant master elements
  numMaEle += numMaSurrounding + numCorner;

  // approximate number of dentries per element
  const int dentries_per_element = parent_node_.GetNumDentries() / numSlEle;

  // approximate number of mentries
  double mentries = dentries_per_element * numMaEle;

  double exp2 = std::ceil(std::log(mentries) / std::log(2.0));
  exp2 = std::max(exp2, 1.0);

  return static_cast<int>(std::pow(2.0, exp2));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::NodeDataContainer::Setup()
{
  if (mentries_ == -1) FOUR_C_THROW("mentries_ must be initialized!");

  //  const int linsize = parentNode_.GetLinsize();
  const int dentries = parent_node_.GetNumDentries();

  d_aug_a_.resize(dentries);
  d_kappa_.resize(dentries);

  CORE::GEN::reset(dentries, dd_aug_a_);
  CORE::GEN::reset(dentries, dd_kappa_);

  d_wgap_sl_.resize(dentries);
  d_wgap_ma_.resize(mentries_);

  d_wgap_sl_complete_ = Teuchos::rcp(new Deriv1stMap(dentries));
  d_wgap_ma_complete_ = Teuchos::rcp(new Deriv1stMap(mentries_));

  CORE::GEN::reset(dentries, dd_wgap_sl_);
  CORE::GEN::reset(mentries_, dd_wgap_ma_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::NodeDataContainer::Pack(CORE::COMM::PackBuffer& data) const
{
  // add maxNumMasterElements
  CORE::COMM::ParObject::AddtoPack(data, mentries_);
  // add kappa_
  CORE::COMM::ParObject::AddtoPack(data, kappa_);
  // add grow_
  CORE::COMM::ParObject::AddtoPack(data, w_gap_);
  // add augA_
  CORE::COMM::ParObject::AddtoPack(data, aug_a_);

  // no need to pack derivs_
  // (these will evaluated new anyway)

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::NodeDataContainer::Unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  // maxNumMasterElements
  CORE::COMM::ParObject::ExtractfromPack(position, data, mentries_);
  // kappa_
  CORE::COMM::ParObject::ExtractfromPack(position, data, kappa_);
  // grow_
  CORE::COMM::ParObject::ExtractfromPack(position, data, w_gap_);
  // augA_
  CORE::COMM::ParObject::ExtractfromPack(position, data, aug_a_);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::NodeDataContainer::Complete()
{
  GetDeriv1st_WGapSl().complete();
  GetDeriv1st_WGapMa().complete();

  GetDeriv2nd_WGapSl().complete();
  GetDeriv2nd_WGapMa().complete();

  GetDeriv1st_WGapSl_Complete().complete();
  GetDeriv1st_WGapMa_Complete().complete();

  GetDeriv1st_Kappa().complete();
  GetDeriv2nd_Kappa().complete();
  GetDeriv1st_A().complete();

  debug_.Complete();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::AUG::NodeDataContainer::Debug::Complete()
{
  d_.complete();
  dd_.complete();

  CORE::GEN::complete(d_vec_);
  CORE::GEN::complete(dd_vec_);
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             ager 08/14|
 *----------------------------------------------------------------------*/
CONTACT::NodePoroDataContainer::NodePoroDataContainer()
{
  ncouprow_ = 0.0;
  for (int i = 0; i < 3; ++i)
  {
    fvel()[i] = 0.0;
    svel()[i] = 0.0;
    poroLM()[i] = 0.0;
  }
  *fpres() = 0.0;
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            ager 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::NodePoroDataContainer::Pack(CORE::COMM::PackBuffer& data) const
{
  // add fvel
  CORE::COMM::ParObject::AddtoPack(data, fvel_, 3 * sizeof(double));
  // add fpres
  CORE::COMM::ParObject::AddtoPack(data, fpres_);
  // add svel
  CORE::COMM::ParObject::AddtoPack(data, svel_, 3 * sizeof(double));
  // add poroLM
  CORE::COMM::ParObject::AddtoPack(data, porolm_, 3 * sizeof(double));
  // add ncoup
  CORE::COMM::ParObject::AddtoPack(data, ncouprow_);
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            ager 08/14|
 *----------------------------------------------------------------------*/
void CONTACT::NodePoroDataContainer::Unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  // fvel
  CORE::COMM::ParObject::ExtractfromPack(position, data, fvel_, 3 * sizeof(double));
  // fpres
  CORE::COMM::ParObject::ExtractfromPack(position, data, fpres_);
  // svel
  CORE::COMM::ParObject::ExtractfromPack(position, data, svel_, 3 * sizeof(double));
  // poroLM
  CORE::COMM::ParObject::ExtractfromPack(position, data, porolm_, 3 * sizeof(double));
  // ncoup
  CORE::COMM::ParObject::ExtractfromPack(position, data, ncouprow_);
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            seitz 08/15|
 *----------------------------------------------------------------------*/
CONTACT::NodeTSIDataContainer::NodeTSIDataContainer(double t_ref, double t_dam)
    : temp_(-1.e12), t_ref_(t_ref), t_dam_(t_dam), thermo_lm_(0.), temp_master_(-1.e12)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           seitz 08/15|
 *----------------------------------------------------------------------*/
void CONTACT::NodeTSIDataContainer::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::ParObject::AddtoPack(data, temp_master_);
  CORE::COMM::ParObject::AddtoPack(data, t_ref_);
  CORE::COMM::ParObject::AddtoPack(data, t_dam_);
  CORE::COMM::ParObject::AddtoPack(data, derivTempMasterDisp_);
  CORE::COMM::ParObject::AddtoPack(data, derivTempMasterTemp_);
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           seitz 08/15|
 *----------------------------------------------------------------------*/
void CONTACT::NodeTSIDataContainer::Unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  CORE::COMM::ParObject::ExtractfromPack(position, data, temp_master_);
  CORE::COMM::ParObject::ExtractfromPack(position, data, t_ref_);
  CORE::COMM::ParObject::ExtractfromPack(position, data, t_dam_);
  CORE::COMM::ParObject::ExtractfromPack(position, data, derivTempMasterDisp_);
  CORE::COMM::ParObject::ExtractfromPack(position, data, derivTempMasterTemp_);
  return;
}

/*----------------------------------------------------------------------*
 |  clear data                                                 (public) |
 |                                                           seitz 08/15|
 *----------------------------------------------------------------------*/
void CONTACT::NodeTSIDataContainer::Clear()
{
  temp_master_ = -1.e12;
  derivTempMasterDisp_.clear();
  derivTempMasterTemp_.clear();
  return;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::Node::Node(int id, const std::vector<double>& coords, const int owner,
    const std::vector<int>& dofs, const bool isslave, const bool initactive)
    : MORTAR::Node(id, coords, owner, dofs, isslave),
      active_(false),
      initactive_(initactive),
      involvedm_(false),
      linsize_(0),  // length of linearization
      codata_(Teuchos::null),
      augdata_(Teuchos::null),
      coporodata_(Teuchos::null),
      cTSIdata_(Teuchos::null)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::Node::Node(const CONTACT::Node& old)
    : MORTAR::Node(old),
      active_(old.active_),
      initactive_(old.initactive_),
      involvedm_(false),
      linsize_(0),
      codata_(Teuchos::null),
      augdata_(Teuchos::null),
      coporodata_(Teuchos::null),
      cTSIdata_(Teuchos::null)
{
  // not yet used and thus not necessarily consistent
  FOUR_C_THROW("Node copy-ctor not yet implemented");

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Node and return pointer to it (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::Node* CONTACT::Node::Clone() const
{
  CONTACT::Node* newnode = new CONTACT::Node(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const CONTACT::Node& cnode)
{
  cnode.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::Node::Print(std::ostream& os) const
{
  // Print id and coordinates
  os << "Contact ";
  MORTAR::Node::Print(os);
  if (IsSlave())
    if (IsInitActive()) os << " InitActive ";
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::Node::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // add base class MORTAR::Node
  MORTAR::Node::Pack(data);

  // add active_
  AddtoPack(data, active_);
  // add initactive_
  AddtoPack(data, initactive_);
  // add involved
  AddtoPack(data, involvedm_);
  // add linsize_
  AddtoPack(data, linsize_);
  // add data_
  bool hasdata = (codata_ != Teuchos::null);
  AddtoPack(data, hasdata);
  if (hasdata) codata_->Pack(data);

  // add augdata_
  bool hasdataaug = (augdata_ != Teuchos::null);
  AddtoPack(data, hasdataaug);
  if (hasdataaug) augdata_->Pack(data);

  // add porodata_
  bool hasdataporo = (coporodata_ != Teuchos::null);
  AddtoPack(data, hasdataporo);
  if (hasdataporo) coporodata_->Pack(data);

  // add tsidata
  bool hasTSIdata = (cTSIdata_ != Teuchos::null);
  AddtoPack(data, (int)hasTSIdata);
  if (hasTSIdata) cTSIdata_->Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::Node::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class MORTAR::Node
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  MORTAR::Node::Unpack(basedata);

  // active_
  active_ = ExtractInt(position, data);
  // isslave_
  initactive_ = ExtractInt(position, data);
  // isslave_
  involvedm_ = ExtractInt(position, data);
  // isslave_
  linsize_ = ExtractInt(position, data);

  // data_
  bool hasdata = ExtractInt(position, data);
  if (hasdata)
  {
    codata_ = Teuchos::rcp(new CONTACT::NodeDataContainer());
    codata_->Unpack(position, data);
  }
  else
  {
    codata_ = Teuchos::null;
  }

  // augdata_
  bool hasdataaug = ExtractInt(position, data);
  if (hasdataaug)
  {
    augdata_ = Teuchos::rcp(new CONTACT::AUG::NodeDataContainer(*this));
    augdata_->Unpack(position, data);
    augdata_->Setup();
  }
  else
    augdata_ = Teuchos::null;

  // porodata_
  bool hasdataporo = ExtractInt(position, data);
  if (hasdataporo)
  {
    coporodata_ = Teuchos::rcp(new CONTACT::NodePoroDataContainer());
    coporodata_->Unpack(position, data);
  }
  else
  {
    coporodata_ = Teuchos::null;
  }

  // TSI data
  bool hasTSIdata = (bool)ExtractInt(position, data);
  if (hasTSIdata)
  {
    cTSIdata_ = Teuchos::rcp(new CONTACT::NodeTSIDataContainer());
    cTSIdata_->Unpack(position, data);
  }
  else
    cTSIdata_ = Teuchos::null;

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the weighted gap                           popp 01/08|
 *----------------------------------------------------------------------*/
void CONTACT::Node::AddgValue(double& val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave() == false) FOUR_C_THROW("AddgValue: function called for master node %i", Id());
  if (IsOnBound() == true) FOUR_C_THROW("AddgValue: function called for boundary node %i", Id());

  // initialize if called for the first time
  if (Data().Getg() == 1.0e12) Data().Getg() = 0.0;

  // add given value to grow_
  Data().Getg() += val;

  return;
}


/*----------------------------------------------------------------------*
 |  Add a value to the weighted gap                      hiermeier 04/14|
 *----------------------------------------------------------------------*/
void CONTACT::Node::AddWGapValue(const double val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave() == false) FOUR_C_THROW("AddWGapValue: function called for master node %i", Id());
  if (IsOnBound() == true) FOUR_C_THROW("AddWGapValue: function called for boundary node %i", Id());

  // initialize if called for the first time
  if (AugData().GetWGap() == 1.0e12) AugData().GetWGap() = 0;

  // add given value to wGap_
  AugData().GetWGap() += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the nts gap                               farah 01/16|
 *----------------------------------------------------------------------*/
void CONTACT::Node::AddntsGapValue(double& val)
{
  // check if this is a master node or slave boundary node
  // initialize if called for the first time
  if (Data().Getgnts() == 1.0e12) Data().Getgnts() = 0;

  // add given value to wGap_
  Data().Getgnts() += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the lts gap                               farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::Node::AddltsGapValue(double& val)
{
  // initialize if called for the first time
  if (Data().Getglts() == 1.0e12) Data().Getglts() = 0;

  // add given value to wGap_
  Data().Getglts() += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the ltl gap                               farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::Node::AddltlGapValue(double* val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave() == false) FOUR_C_THROW("function called for master node %i", Id());
  if (!IsOnEdge()) FOUR_C_THROW("function call for non edge node! %i", Id());

  // initialize if called for the first time
  if (Data().Getgltl()[0] == 1.0e12 or Data().Getgltl()[1] == 1.0e12 or
      Data().Getgltl()[2] == 1.0e12)
  {
    Data().Getgltl()[0] = 0.0;
    Data().Getgltl()[1] = 0.0;
    Data().Getgltl()[2] = 0.0;
  }

  // add given value to wGap_
  Data().Getgltl()[0] += val[0];
  Data().Getgltl()[1] += val[1];
  Data().Getgltl()[2] += val[2];

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the ltl jump                               farah 07/16|
 *----------------------------------------------------------------------*/
void CONTACT::Node::AddltlJumpValue(double* val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave() == false) FOUR_C_THROW("function called for master node %i", Id());
  if (!IsOnEdge()) FOUR_C_THROW("function call for non edge node! %i", Id());

  // initialize if called for the first time
  if (Data().Getjumpltl()[0] == 1.0e12 or Data().Getjumpltl()[1] == 1.0e12 or
      Data().Getjumpltl()[2] == 1.0e12)
  {
    Data().Getjumpltl()[0] = 0.0;
    Data().Getjumpltl()[1] = 0.0;
    Data().Getjumpltl()[2] = 0.0;
  }

  // add given value to wGap_
  Data().Getjumpltl()[0] += val[0];
  Data().Getjumpltl()[1] += val[1];
  Data().Getjumpltl()[2] += val[2];

  return;
}


/*----------------------------------------------------------------------*
 |  Add a value to scaling factor kappa                  hiermeier 04/14|
 *----------------------------------------------------------------------*/
void CONTACT::Node::AddKappaValue(double& val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave() == false) FOUR_C_THROW("AddKappaValue: function called for master node %i", Id());
  if (IsOnBound() == true)
    FOUR_C_THROW("AddKappaValue: function called for boundary node %i", Id());

  // initialize if called for the first time
  if (AugData().GetKappa() == 1.0e12) AugData().GetKappa() = 0.0;

  // add given value to kappa_
  AugData().GetKappa() += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the 'DerivZ' map                           popp 06/09|
 *----------------------------------------------------------------------*/
void CONTACT::Node::AddDerivZValue(int& row, const int& col, double val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave() == false) FOUR_C_THROW("AddZValue: function called for master node %i", Id());
  if (IsOnBound() == true) FOUR_C_THROW("AddZValue: function called for boundary node %i", Id());

  // check if this has been called before
  if ((int)Data().GetDerivZ().size() == 0) Data().GetDerivZ().resize(NumDof());

  // check row index input
  if ((int)Data().GetDerivZ().size() <= row)
    FOUR_C_THROW("AddDerivZValue: tried to access invalid row index!");

  // add the pair (col,val) to the given row
  std::map<int, double>& zmap = Data().GetDerivZ()[row];
  zmap[col] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  Initialize data container                             gitterle 02/10|
 *----------------------------------------------------------------------*/
void CONTACT::Node::InitializeDataContainer()
{
  // get maximum size of lin vectors
  linsize_ = 0;
  // get maximum size of nodal D-entries
  dentries_ = 0;
  std::set<int> sIdCheck;
  std::pair<std::set<int>::iterator, bool> check;
  for (int i = 0; i < NumElement(); ++i)
  {
    const int* snodeIds = Elements()[i]->NodeIds();
    const DRT::Node* const* nodes = Elements()[i]->Nodes();
    for (int j = 0; j < Elements()[i]->NumNode(); ++j)
    {
      const int numdof = Elements()[i]->NumDofPerNode(*(nodes[j]));

      check = sIdCheck.insert(snodeIds[j]);
      if (check.second)
      {
        dentries_ += numdof;
      }
      linsize_ += numdof;
    }
  }
  // set a minimal number for pure ghost nodes
  if (NumElement() == 0)
  {
    dentries_ = 3;
    linsize_ = 3;
  }

  // only initialize if not yet done
  if (modata_ == Teuchos::null && codata_ == Teuchos::null)
  {
    codata_ = Teuchos::rcp(new CONTACT::NodeDataContainer());
    modata_ = Teuchos::rcp(new MORTAR::NodeDataContainer());
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::Node::InitializeAugDataContainer(
    const int slMaElementAreaRatio, const bool isTriangleOnMaster)
{
  // Do it only, if the container has not been initialized, yet.
  if (augdata_.is_null())
  {
    augdata_ = Teuchos::rcp(
        new CONTACT::AUG::NodeDataContainer(*this, slMaElementAreaRatio, isTriangleOnMaster));

    augdata_->Setup();
  }
}

/*-----------------------------------------------------------------------*
 |  Initialize poro data container                             ager 07/14|
 *----------------------------------------------------------------------*/
void CONTACT::Node::InitializePoroDataContainer()
{
  // only initialize if not yet done

  if (coporodata_ == Teuchos::null)
  {
    coporodata_ = Teuchos::rcp(new CONTACT::NodePoroDataContainer());
  }

  return;
}

/*-----------------------------------------------------------------------*
 |  Initialize ehl data container                             seitz 11/17|
 *----------------------------------------------------------------------*/
void CONTACT::Node::InitializeEhlDataContainer()
{
  // only initialize if not yet done

  if (cEHLdata_ == Teuchos::null)
  {
    cEHLdata_ = Teuchos::rcp(new CONTACT::NodeEhlDataContainer());
  }

  return;
}

/*-----------------------------------------------------------------------*
 |  Initialize TSI data container                             seitz 08/15|
 *----------------------------------------------------------------------*/
void CONTACT::Node::InitializeTSIDataContainer(double t_ref, double t_dam)
{
  // only initialize if not yet done

  if (cTSIdata_ == Teuchos::null)
    cTSIdata_ = Teuchos::rcp(new CONTACT::NodeTSIDataContainer(t_ref, t_dam));

  return;
}

/*----------------------------------------------------------------------*
 |  Reset data container                                      popp 09/10|
 *----------------------------------------------------------------------*/
void CONTACT::Node::ResetDataContainer()
{
  // reset to Teuchos::null
  codata_ = Teuchos::null;
  modata_ = Teuchos::null;
  coporodata_ = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------*
 |  Build averaged nodal edge tangents                       farah 11/16|
 *----------------------------------------------------------------------*/
void CONTACT::Node::BuildAveragedEdgeTangent()
{
  for (int j = 0; j < 3; ++j)
  {
    MoData().EdgeTangent()[j] = 0.0;
  }
  int nseg = NumElement();
  DRT::Element** adjeles = Elements();

  //**************************************************
  //              CALCULATE EDGES
  //**************************************************
  // empty vector of slave element pointers
  std::vector<Teuchos::RCP<MORTAR::Element>> lineElementsS;
  std::set<std::pair<int, int>> donebefore;

  // loop over all surface elements
  for (int surfele = 0; surfele < nseg; ++surfele)
  {
    Element* cele = dynamic_cast<Element*>(adjeles[surfele]);

    if (cele->Shape() == CORE::FE::CellType::quad4)
    {
      for (int j = 0; j < 4; ++j)
      {
        int nodeIds[2] = {0, 0};
        int nodeLIds[2] = {0, 0};

        if (j == 0)
        {
          nodeIds[0] = cele->NodeIds()[0];
          nodeIds[1] = cele->NodeIds()[1];

          nodeLIds[0] = 0;
          nodeLIds[1] = 1;
        }
        else if (j == 1)
        {
          nodeIds[0] = cele->NodeIds()[1];
          nodeIds[1] = cele->NodeIds()[2];

          nodeLIds[0] = 1;
          nodeLIds[1] = 2;
        }
        else if (j == 2)
        {
          nodeIds[0] = cele->NodeIds()[2];
          nodeIds[1] = cele->NodeIds()[3];

          nodeLIds[0] = 2;
          nodeLIds[1] = 3;
        }
        else if (j == 3)
        {
          nodeIds[0] = cele->NodeIds()[3];
          nodeIds[1] = cele->NodeIds()[0];

          nodeLIds[0] = 3;
          nodeLIds[1] = 0;
        }

        // check if both nodes on edge geometry
        bool node0Edge = dynamic_cast<MORTAR::Node*>(cele->Nodes()[nodeLIds[0]])->IsOnEdge();
        bool node1Edge = dynamic_cast<MORTAR::Node*>(cele->Nodes()[nodeLIds[1]])->IsOnEdge();

        if (!node0Edge or !node1Edge) continue;

        // create pair
        std::pair<int, int> actIDs = std::pair<int, int>(nodeIds[0], nodeIds[1]);
        std::pair<int, int> actIDstw = std::pair<int, int>(nodeIds[1], nodeIds[0]);

        // check if processed before
        std::set<std::pair<int, int>>::iterator iter = donebefore.find(actIDs);
        std::set<std::pair<int, int>>::iterator itertw = donebefore.find(actIDstw);

        // if not then create ele
        if (iter == donebefore.end() and itertw == donebefore.end())
        {
          // add to set of processed nodes
          donebefore.insert(actIDs);
          donebefore.insert(actIDstw);

          // create line ele:
          Teuchos::RCP<MORTAR::Element> lineEle = Teuchos::rcp(
              new MORTAR::Element(j, cele->Owner(), CORE::FE::CellType::line2, 2, nodeIds, false));

          // get nodes
          std::array<DRT ::Node*, 2> nodes = {
              cele->Nodes()[nodeLIds[0]], cele->Nodes()[nodeLIds[1]]};
          lineEle->BuildNodalPointers(nodes.data());

          // init data container for dual shapes
          lineEle->InitializeDataContainer();

          // push back into vector
          lineElementsS.push_back(lineEle);
        }
      }  // end edge loop
    }
    else
      FOUR_C_THROW("only quad4!");
  }  // end surfele loop


  //**************************************************
  //      GET EDGES WHICH ARE CONNECTED TO NODE
  //**************************************************
  std::vector<int> dummy;
  for (int edges = 0; edges < (int)lineElementsS.size(); ++edges)
  {
    for (int count = 0; count < lineElementsS[edges]->NumNode(); ++count)
    {
      if (lineElementsS[edges]->Nodes()[count]->Id() == Id()) dummy.push_back(edges);
    }
  }

  if (dummy.size() > 2) std::cout << "WARNING: multiple edge definitions possible!" << std::endl;

  if (dummy.size() < 1) FOUR_C_THROW("ERROR!");

  //**************************************************
  //      CALC ADJACENT TANGENTS
  //**************************************************
  Node* n1 = nullptr;
  Node* n2 = nullptr;

  if (lineElementsS[dummy[0]]->Nodes()[0]->Id() != Id())
    n1 = dynamic_cast<Node*>(lineElementsS[dummy[0]]->Nodes()[0]);
  else if (lineElementsS[dummy[0]]->Nodes()[1]->Id() != Id())
    n1 = dynamic_cast<Node*>(lineElementsS[dummy[0]]->Nodes()[1]);
  else
    FOUR_C_THROW("ERROR");

  if (dummy.size() < 2)
  {
    // get node itself if no adjacent edge elements could be found
    n2 = this;
  }
  else
  {
    if (lineElementsS[dummy[1]]->Nodes()[0]->Id() != Id())
      n2 = dynamic_cast<Node*>(lineElementsS[dummy[0]]->Nodes()[0]);
    else if (lineElementsS[dummy[1]]->Nodes()[1]->Id() != Id())
      n2 = dynamic_cast<Node*>(lineElementsS[dummy[0]]->Nodes()[1]);
    else
      FOUR_C_THROW("ERROR");
  }


  std::array<double, 3> tmp1 = {0.0, 0.0, 0.0};
  // difference
  tmp1[0] = n1->xspatial()[0] - n2->xspatial()[0];
  tmp1[1] = n1->xspatial()[1] - n2->xspatial()[1];
  tmp1[2] = n1->xspatial()[2] - n2->xspatial()[2];
  const double length = sqrt(tmp1[0] * tmp1[0] + tmp1[1] * tmp1[1] + tmp1[2] * tmp1[2]);
  if (length < 1e-12) FOUR_C_THROW("ERROR");
  MoData().EdgeTangent()[0] = tmp1[0] / length;
  MoData().EdgeTangent()[1] = tmp1[1] / length;
  MoData().EdgeTangent()[2] = tmp1[2] / length;

  //**************************************************
  //      LINEARIZATION
  //**************************************************
  typedef CORE::GEN::Pairedvector<int, double>::const_iterator _CI;

  for (int j = 0; j < (int)((Data().GetDerivTangent()).size()); ++j)
    (Data().GetDerivTangent())[j].clear();
  (Data().GetDerivTangent()).resize(0, 0);
  if ((int)Data().GetDerivTangent().size() == 0) Data().GetDerivTangent().resize(3, 2 * 100);

  std::vector<CORE::GEN::Pairedvector<int, double>> lint(3, 100);  // added all sizes
  if (n1 != nullptr)
  {
    lint[0][n1->Dofs()[0]] += 1;
    lint[1][n1->Dofs()[1]] += 1;
    lint[2][n1->Dofs()[2]] += 1;
  }
  if (n2 != nullptr)
  {
    lint[0][n2->Dofs()[0]] -= 1;
    lint[1][n2->Dofs()[1]] -= 1;
    lint[2][n2->Dofs()[2]] -= 1;
  }
  // first part
  std::vector<CORE::GEN::Pairedvector<int, double>> Lin1(3, 100);  // added all sizes
  for (_CI p = lint[0].begin(); p != lint[0].end(); ++p) Lin1[0][p->first] += p->second / length;
  for (_CI p = lint[1].begin(); p != lint[1].end(); ++p) Lin1[1][p->first] += p->second / length;
  for (_CI p = lint[2].begin(); p != lint[2].end(); ++p) Lin1[2][p->first] += p->second / length;

  CORE::GEN::Pairedvector<int, double> Lin2(100);  // added all sizes
  for (_CI p = lint[0].begin(); p != lint[0].end(); ++p)
    Lin2[p->first] += p->second * MoData().EdgeTangent()[0];
  for (_CI p = lint[1].begin(); p != lint[1].end(); ++p)
    Lin2[p->first] += p->second * MoData().EdgeTangent()[1];
  for (_CI p = lint[2].begin(); p != lint[2].end(); ++p)
    Lin2[p->first] += p->second * MoData().EdgeTangent()[2];

  std::vector<CORE::GEN::Pairedvector<int, double>> Lin3(3, 100);  // added all sizes
  for (_CI p = Lin2.begin(); p != Lin2.end(); ++p)
    Lin3[0][p->first] += p->second * MoData().EdgeTangent()[0] / (length * length * length);
  for (_CI p = Lin2.begin(); p != Lin2.end(); ++p)
    Lin3[1][p->first] += p->second * MoData().EdgeTangent()[1] / (length * length * length);
  for (_CI p = Lin2.begin(); p != Lin2.end(); ++p)
    Lin3[2][p->first] += p->second * MoData().EdgeTangent()[2] / (length * length * length);

  for (_CI p = Lin1[0].begin(); p != Lin1[0].end(); ++p)
    Data().GetDerivTangent()[0][p->first] += p->second;
  for (_CI p = Lin1[1].begin(); p != Lin1[1].end(); ++p)
    Data().GetDerivTangent()[1][p->first] += p->second;
  for (_CI p = Lin1[2].begin(); p != Lin1[2].end(); ++p)
    Data().GetDerivTangent()[2][p->first] += p->second;

  for (_CI p = Lin3[0].begin(); p != Lin3[0].end(); ++p)
    Data().GetDerivTangent()[0][p->first] -= p->second;
  for (_CI p = Lin3[1].begin(); p != Lin3[1].end(); ++p)
    Data().GetDerivTangent()[1][p->first] -= p->second;
  for (_CI p = Lin3[2].begin(); p != Lin3[2].end(); ++p)
    Data().GetDerivTangent()[2][p->first] -= p->second;

  // std::cout << "tangent = " << MoData().EdgeTangent()[0] << "  " << MoData().EdgeTangent()[1] <<
  // "  " << MoData().EdgeTangent()[2] << std::endl;

  // bye bye
  return;
}

/*----------------------------------------------------------------------*
 |  Build averaged nodal normal + tangents                    popp 12/07|
 *----------------------------------------------------------------------*/
void CONTACT::Node::BuildAveragedNormal()
{
  // reset normal and tangents when this method is called
  for (int j = 0; j < 3; ++j)
  {
    MoData().n()[j] = 0.0;
    Data().txi()[j] = 0.0;
    Data().teta()[j] = 0.0;
  }

  int nseg = NumElement();
  DRT::Element** adjeles = Elements();

  // temporary vector to store nodal normal
  std::array<double, 3> n_tmp = {0., 0., 0.};
  CORE::LINALG::SerialDenseMatrix elens(6, nseg);

  // we need to store some stuff here
  //**********************************************************************
  // elens(0,i): x-coord of element normal
  // elens(1,i): y-coord of element normal
  // elens(2,i): z-coord of element normal
  // elens(3,i): id of adjacent element i
  // elens(4,i): length of element normal
  // elens(5,i): length/area of element itself
  //**********************************************************************

  // loop over all adjacent elements
  for (int i = 0; i < nseg; ++i)
  {
    Element* adjcele = dynamic_cast<Element*>(adjeles[i]);

    // build element normal at current node
    // (we have to pass in the index i to be able to store the
    // normal and other information at the right place in elens)
    adjcele->BuildNormalAtNode(Id(), i, elens);

    // add (weighted) element normal to nodal normal n
    for (int j = 0; j < 3; ++j) n_tmp[j] += elens(j, i) / elens(4, i);
  }

  // modify normal in case of symmetry condition
  for (int i = 0; i < 3; i++)
    if (DbcDofs()[i]) n_tmp[i] = 0.;

  // create unit normal vector
  double length = sqrt(n_tmp[0] * n_tmp[0] + n_tmp[1] * n_tmp[1] + n_tmp[2] * n_tmp[2]);
  if (length < 1e-12)
  {
    std::cout << "normal zero: node slave= " << IsSlave() << "  length= " << length << std::endl;
    FOUR_C_THROW("Nodal normal length 0, node ID %i", Id());
  }
  else
  {
    for (int j = 0; j < 3; ++j)
    {
      MoData().n()[j] = n_tmp[j] / length;
    }
  }

  // create unit tangent vectors
  // (note that this definition is not unique in 3D!)
  double ltxi = 1.0;

  if (NumDof() == 2)
  {
    // simple definition for txi
    Data().txi()[0] = -MoData().n()[1];
    Data().txi()[1] = MoData().n()[0];
    Data().txi()[2] = 0.0;

    // teta is z-axis
    Data().teta()[0] = 0.0;
    Data().teta()[1] = 0.0;
    Data().teta()[2] = 1.0;
  }
  else if (NumDof() == 3)
  {
#ifdef CONTACTPSEUDO2D
    // we want to treat a 3D mesh as pseudo 2D contact problem
    // with all nodes fixed in z-direction
    // thus, the second tangent is fixed to (0,0,1)
    Data().teta()[0] = 0.0;
    Data().teta()[1] = 0.0;
    Data().teta()[2] = 1.0;

    // txi follows from corkscrew rule (txi = teta x n)
    Data().txi()[0] = Data().teta()[1] * MoData().n()[2] - Data().teta()[2] * MoData().n()[1];
    Data().txi()[1] = Data().teta()[2] * MoData().n()[0] - Data().teta()[0] * MoData().n()[2];
    Data().txi()[2] = Data().teta()[0] * MoData().n()[1] - Data().teta()[1] * MoData().n()[0];
#else

    if (abs(MoData().n()[0]) > 1.0e-4 || abs(MoData().n()[1]) > 1.0e-4)
    {
      Data().txi()[0] = -MoData().n()[1];
      Data().txi()[1] = MoData().n()[0];
      Data().txi()[2] = 0.0;
    }
    else
    {
      Data().txi()[0] = 0.0;
      Data().txi()[1] = -MoData().n()[2];
      Data().txi()[2] = MoData().n()[1];
    }

    ltxi = sqrt(Data().txi()[0] * Data().txi()[0] + Data().txi()[1] * Data().txi()[1] +
                Data().txi()[2] * Data().txi()[2]);
    if (ltxi < 1e-12)
    {
      std::cout << "tangent 1 zero: node slave= " << IsSlave() << "  length= " << ltxi << std::endl;
      FOUR_C_THROW("Nodal tangent length 0, node ID %i", Id());
    }
    else
    {
      for (int j = 0; j < 3; ++j) Data().txi()[j] /= ltxi;
    }



    // teta follows from corkscrew rule (teta = n x txi)
    Data().teta()[0] = MoData().n()[1] * Data().txi()[2] - MoData().n()[2] * Data().txi()[1];
    Data().teta()[1] = MoData().n()[2] * Data().txi()[0] - MoData().n()[0] * Data().txi()[2];
    Data().teta()[2] = MoData().n()[0] * Data().txi()[1] - MoData().n()[1] * Data().txi()[0];

#endif
  }
  else
    FOUR_C_THROW("Contact problems must be either 2D or 3D");

  // build linearization of averaged nodal normal and tangents
  DerivAveragedNormal(elens, length, ltxi);

  return;
}

/*----------------------------------------------------------------------*
 |  Build directional deriv. of nodal normal + tangents       popp 09/08|
 *----------------------------------------------------------------------*/
void CONTACT::Node::DerivAveragedNormal(
    CORE::LINALG::SerialDenseMatrix& elens, double length, double ltxi)
{
  int nseg = NumElement();
  DRT::Element** adjeles = Elements();

  // prepare nodal storage maps for derivative
  if ((int)Data().GetDerivN().size() == 0) Data().GetDerivN().resize(3, linsize_);
  if ((int)Data().GetDerivTxi().size() == 0) Data().GetDerivTxi().resize(3, linsize_);
  if ((int)Data().GetDerivTeta().size() == 0) Data().GetDerivTeta().resize(3, linsize_);

  // loop over all adjacent elements
  for (int i = 0; i < nseg; ++i)
  {
    Element* adjcele = dynamic_cast<Element*>(adjeles[i]);

    // build element normal derivative at current node
    adjcele->DerivNormalAtNode(Id(), i, elens, Data().GetDerivN());
  }

  // modify normal in case of symmetry condition
  for (int i = 0; i < 3; i++)
    if (DbcDofs()[i]) Data().GetDerivN()[i].clear();

  // normalize directional derivative
  // (length differs for weighted/unweighted case but not the procedure!)
  // (be careful with reference / copy of derivative maps!)
  typedef CORE::GEN::Pairedvector<int, double>::const_iterator CI;
  CORE::GEN::Pairedvector<int, double>& derivnx = Data().GetDerivN()[0];
  CORE::GEN::Pairedvector<int, double>& derivny = Data().GetDerivN()[1];
  CORE::GEN::Pairedvector<int, double>& derivnz = Data().GetDerivN()[2];
  CORE::GEN::Pairedvector<int, double> cderivnx = Data().GetDerivN()[0];
  CORE::GEN::Pairedvector<int, double> cderivny = Data().GetDerivN()[1];
  CORE::GEN::Pairedvector<int, double> cderivnz = Data().GetDerivN()[2];
  const double nxnx = MoData().n()[0] * MoData().n()[0];
  const double nxny = MoData().n()[0] * MoData().n()[1];
  const double nxnz = MoData().n()[0] * MoData().n()[2];
  const double nyny = MoData().n()[1] * MoData().n()[1];
  const double nynz = MoData().n()[1] * MoData().n()[2];
  const double nznz = MoData().n()[2] * MoData().n()[2];

  // build a vector with all keys from x,y,z maps
  // (we need this in order not to miss any entry!)
  std::vector<int> allkeysn;
  for (CI p = derivnx.begin(); p != derivnx.end(); ++p)
  {
    bool found = false;
    for (int j = 0; j < (int)allkeysn.size(); ++j)
      if ((p->first) == allkeysn[j]) found = true;
    if (!found) allkeysn.push_back(p->first);
  }
  for (CI p = derivny.begin(); p != derivny.end(); ++p)
  {
    bool found = false;
    for (int j = 0; j < (int)allkeysn.size(); ++j)
      if ((p->first) == allkeysn[j]) found = true;
    if (!found) allkeysn.push_back(p->first);
  }
  for (CI p = derivnz.begin(); p != derivnz.end(); ++p)
  {
    bool found = false;
    for (int j = 0; j < (int)allkeysn.size(); ++j)
      if ((p->first) == allkeysn[j]) found = true;
    if (!found) allkeysn.push_back(p->first);
  }

  // normalize x-components
  for (int j = 0; j < (int)allkeysn.size(); ++j)
  {
    double val = cderivnx[allkeysn[j]];
    derivnx[allkeysn[j]] =
        (val - nxnx * val - nxny * cderivny[allkeysn[j]] - nxnz * cderivnz[allkeysn[j]]) / length;
  }

  // normalize y-components
  for (int j = 0; j < (int)allkeysn.size(); ++j)
  {
    double val = cderivny[allkeysn[j]];
    derivny[allkeysn[j]] =
        (val - nxny * cderivnx[allkeysn[j]] - nyny * val - nynz * cderivnz[allkeysn[j]]) / length;
  }

  // normalize z-components
  for (int j = 0; j < (int)allkeysn.size(); ++j)
  {
    double val = cderivnz[allkeysn[j]];
    derivnz[allkeysn[j]] =
        (val - nxnz * cderivnx[allkeysn[j]] - nynz * cderivny[allkeysn[j]] - nznz * val) / length;
  }

  //**********************************************************************
  // tangent derivatives 2D
  //**********************************************************************
  if (NumDof() == 2)
  {
    // get directional derivative of nodal tangent txi "for free"
    // (we just have to use the orthogonality of n and t)
    // the directional derivative of nodal tangent teta is 0
    CORE::GEN::Pairedvector<int, double>& derivtxix = Data().GetDerivTxi()[0];
    CORE::GEN::Pairedvector<int, double>& derivtxiy = Data().GetDerivTxi()[1];

    for (CI p = derivny.begin(); p != derivny.end(); ++p) derivtxix[p->first] = -(p->second);
    for (CI p = derivnx.begin(); p != derivnx.end(); ++p) derivtxiy[p->first] = (p->second);
  }

  //**********************************************************************
  // tangent derivatives 3D
  //**********************************************************************
  else
  {
#ifdef CONTACTPSEUDO2D
    // trivial tangent derivative teta
    // this is 0 as teta is fixed to (0,0,1)

    // get normalized tangent derivative txi
    // use corkscrew rule from BuildAveragedNormal()
    CORE::GEN::Pairedvector<int, double>& derivtxix = Data().GetDerivTxi()[0];
    CORE::GEN::Pairedvector<int, double>& derivtxiy = Data().GetDerivTxi()[1];
    CORE::GEN::Pairedvector<int, double>& derivtxiz = Data().GetDerivTxi()[2];

    for (CI p = derivnx.begin(); p != derivnx.end(); ++p)
    {
      derivtxiy[p->first] += Data().teta()[2] * (p->second);
      derivtxiz[p->first] -= Data().teta()[1] * (p->second);
    }
    for (CI p = derivny.begin(); p != derivny.end(); ++p)
    {
      derivtxix[p->first] -= Data().teta()[2] * (p->second);
      derivtxiz[p->first] += Data().teta()[0] * (p->second);
    }
    for (CI p = derivnz.begin(); p != derivnz.end(); ++p)
    {
      derivtxix[p->first] += Data().teta()[1] * (p->second);
      derivtxiy[p->first] -= Data().teta()[0] * (p->second);
    }
  }
#else
    // unnormalized tangent derivative txi
    // use definitions for txi from BuildAveragedNormal()
    if (abs(MoData().n()[0]) > 1.0e-4 || abs(MoData().n()[1]) > 1.0e-4)
    {
      CORE::GEN::Pairedvector<int, double>& derivtxix = Data().GetDerivTxi()[0];
      CORE::GEN::Pairedvector<int, double>& derivtxiy = Data().GetDerivTxi()[1];

      for (CI p = derivny.begin(); p != derivny.end(); ++p) derivtxix[p->first] -= (p->second);

      for (CI p = derivnx.begin(); p != derivnx.end(); ++p) derivtxiy[p->first] += (p->second);
    }
    else
    {
      CORE::GEN::Pairedvector<int, double>& derivtxiy = Data().GetDerivTxi()[1];
      CORE::GEN::Pairedvector<int, double>& derivtxiz = Data().GetDerivTxi()[2];

      for (CI p = derivnz.begin(); p != derivnz.end(); ++p) derivtxiy[p->first] -= (p->second);

      for (CI p = derivny.begin(); p != derivny.end(); ++p) derivtxiz[p->first] += (p->second);
    }

    // normalize txi directional derivative
    // (identical to normalization of normal derivative)
    typedef CORE::GEN::Pairedvector<int, double>::const_iterator CI;
    CORE::GEN::Pairedvector<int, double>& derivtxix = Data().GetDerivTxi()[0];
    CORE::GEN::Pairedvector<int, double>& derivtxiy = Data().GetDerivTxi()[1];
    CORE::GEN::Pairedvector<int, double>& derivtxiz = Data().GetDerivTxi()[2];
    CORE::GEN::Pairedvector<int, double> cderivtxix = Data().GetDerivTxi()[0];
    CORE::GEN::Pairedvector<int, double> cderivtxiy = Data().GetDerivTxi()[1];
    CORE::GEN::Pairedvector<int, double> cderivtxiz = Data().GetDerivTxi()[2];
    const double txtx = Data().txi()[0] * Data().txi()[0];
    const double txty = Data().txi()[0] * Data().txi()[1];
    const double txtz = Data().txi()[0] * Data().txi()[2];
    const double tyty = Data().txi()[1] * Data().txi()[1];
    const double tytz = Data().txi()[1] * Data().txi()[2];
    const double tztz = Data().txi()[2] * Data().txi()[2];

    // build a vector with all keys from x,y,z maps
    // (we need this in order not to miss any entry!)
    std::vector<int> allkeyst;
    for (CI p = derivtxix.begin(); p != derivtxix.end(); ++p)
    {
      bool found = false;
      for (int j = 0; j < (int)allkeyst.size(); ++j)
        if ((p->first) == allkeyst[j]) found = true;
      if (!found) allkeyst.push_back(p->first);
    }
    for (CI p = derivtxiy.begin(); p != derivtxiy.end(); ++p)
    {
      bool found = false;
      for (int j = 0; j < (int)allkeyst.size(); ++j)
        if ((p->first) == allkeyst[j]) found = true;
      if (!found) allkeyst.push_back(p->first);
    }
    for (CI p = derivtxiz.begin(); p != derivtxiz.end(); ++p)
    {
      bool found = false;
      for (int j = 0; j < (int)allkeyst.size(); ++j)
        if ((p->first) == allkeyst[j]) found = true;
      if (!found) allkeyst.push_back(p->first);
    }

    // normalize x-components
    for (int j = 0; j < (int)allkeyst.size(); ++j)
    {
      double val = cderivtxix[allkeyst[j]];
      derivtxix[allkeyst[j]] =
          (val - txtx * val - txty * cderivtxiy[allkeyst[j]] - txtz * cderivtxiz[allkeyst[j]]) /
          ltxi;
    }

    // normalize y-components
    for (int j = 0; j < (int)allkeyst.size(); ++j)
    {
      double val = cderivtxiy[allkeyst[j]];
      derivtxiy[allkeyst[j]] =
          (val - txty * cderivtxix[allkeyst[j]] - tyty * val - tytz * cderivtxiz[allkeyst[j]]) /
          ltxi;
    }

    // normalize z-components
    for (int j = 0; j < (int)allkeyst.size(); ++j)
    {
      double val = cderivtxiz[allkeyst[j]];
      derivtxiz[allkeyst[j]] =
          (val - txtz * cderivtxix[allkeyst[j]] - tytz * cderivtxiy[allkeyst[j]] - tztz * val) /
          ltxi;
    }

    // get normalized tangent derivative teta
    // use corkscrew rule from BuildAveragedNormal()
    CORE::GEN::Pairedvector<int, double>& derivtetax = Data().GetDerivTeta()[0];
    CORE::GEN::Pairedvector<int, double>& derivtetay = Data().GetDerivTeta()[1];
    CORE::GEN::Pairedvector<int, double>& derivtetaz = Data().GetDerivTeta()[2];

    for (CI p = derivnx.begin(); p != derivnx.end(); ++p)
    {
      derivtetay[p->first] -= Data().txi()[2] * (p->second);
      derivtetaz[p->first] += Data().txi()[1] * (p->second);
    }
    for (CI p = derivny.begin(); p != derivny.end(); ++p)
    {
      derivtetax[p->first] += Data().txi()[2] * (p->second);
      derivtetaz[p->first] -= Data().txi()[0] * (p->second);
    }
    for (CI p = derivnz.begin(); p != derivnz.end(); ++p)
    {
      derivtetax[p->first] -= Data().txi()[1] * (p->second);
      derivtetay[p->first] += Data().txi()[0] * (p->second);
    }
    for (CI p = derivtxix.begin(); p != derivtxix.end(); ++p)
    {
      derivtetay[p->first] += MoData().n()[2] * (p->second);
      derivtetaz[p->first] -= MoData().n()[1] * (p->second);
    }
    for (CI p = derivtxiy.begin(); p != derivtxiy.end(); ++p)
    {
      derivtetax[p->first] -= MoData().n()[2] * (p->second);
      derivtetaz[p->first] += MoData().n()[0] * (p->second);
    }
    for (CI p = derivtxiz.begin(); p != derivtxiz.end(); ++p)
    {
      derivtetax[p->first] += MoData().n()[1] * (p->second);
      derivtetay[p->first] -= MoData().n()[0] * (p->second);
    }
  }
#endif  // #ifdef CONTACTPSEUDO2D

  return;
}

/*----------------------------------------------------------------------*
 |  Add a value to the NCoup of this node                      ager 06/14|
 *----------------------------------------------------------------------*/
void CONTACT::Node::AddNcoupValue(double& val)
{
  // check if this is a master node or slave boundary node
  if (IsSlave() == false) FOUR_C_THROW("AddNcoupValue: function called for master node %i", Id());
  if (IsOnBound() == true)
    FOUR_C_THROW("AddNcoupValue: function called for boundary node %i", Id());

  // add given value to ncoup
  PoroData().GetnCoup() += val;
  return;
}

/*----------------------------------------------------------------------*
 |  Store nodal normals to old ones                         seitz 05/17 |
 *----------------------------------------------------------------------*/
void CONTACT::Node::StoreOldNormal()
{
  // write entries to old ones
  for (int j = 0; j < 3; ++j) Data().Normal_old()[j] = MoData().n()[j];

  return;
}

FOUR_C_NAMESPACE_CLOSE
