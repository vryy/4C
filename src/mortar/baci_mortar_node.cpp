/*-----------------------------------------------------------------------*/
/*! \file
\brief A class for a mortar coupling node

\level 2

*/
/*-----------------------------------------------------------------------*/

#include "baci_mortar_node.H"

#include "baci_lib_discret.H"
#include "baci_mortar_element.H"
#include "baci_utils_exceptions.H"


MORTAR::MortarNodeType MORTAR::MortarNodeType::instance_;


CORE::COMM::ParObject* MORTAR::MortarNodeType::Create(const std::vector<char>& data)
{
  std::vector<double> x(3, 0.0);
  std::vector<int> dofs(0);
  auto* node = new MORTAR::MortarNode(0, x, 0, dofs, false);
  node->Unpack(data);
  return node;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MORTAR::MortarNodeDataContainer::MortarNodeDataContainer()
    : drows_(0), drows_nts_(0), drows_lts_(0), drows_ltl_(0)
{
  for (int i = 0; i < 3; ++i)
  {
    n()[i] = 0.0;
    EdgeTangent()[i] = 0.0;
    lm()[i] = 0.0;
    lmold()[i] = 0.0;
    lmuzawa()[i] = 0.0;
    GetDscale() = 0.0;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarNodeDataContainer::Pack(CORE::COMM::PackBuffer& data) const
{
  // add n_
  CORE::COMM::ParObject::AddtoPack(data, n_, 3 * sizeof(double));
  // add edgetangent_
  CORE::COMM::ParObject::AddtoPack(data, edgeTangent_, 3 * sizeof(double));
  // add lm_
  CORE::COMM::ParObject::AddtoPack(data, lm_, 3 * sizeof(double));
  // add lmold_
  CORE::COMM::ParObject::AddtoPack(data, lmold_, 3 * sizeof(double));
  // add lmuzawa_
  CORE::COMM::ParObject::AddtoPack(data, lmuzawa_, 3 * sizeof(double));

  // no need to pack drows_, mrows_ and mmodrows_
  // (these will evaluated anew anyway)
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarNodeDataContainer::Unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  // n_
  CORE::COMM::ParObject::ExtractfromPack(position, data, n_, 3 * sizeof(double));
  // edgetangent_
  CORE::COMM::ParObject::ExtractfromPack(position, data, edgeTangent_, 3 * sizeof(double));
  // lm_
  CORE::COMM::ParObject::ExtractfromPack(position, data, lm_, 3 * sizeof(double));
  // lmold_
  CORE::COMM::ParObject::ExtractfromPack(position, data, lmold_, 3 * sizeof(double));
  // lmuzawa_
  CORE::COMM::ParObject::ExtractfromPack(position, data, lmuzawa_, 3 * sizeof(double));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MORTAR::MortarNode::MortarNode(int id, const std::vector<double>& coords, const int owner,
    const std::vector<int>& dofs, const bool isslave)
    : DRT::Node(id, coords, owner),
      isslave_(isslave),
      istiedslave_(isslave),
      isonbound_(false),
      isonedge_(false),
      isoncorner_(false),
      isdbc_(false),
      dofs_(dofs),
      hasproj_(false),
      hassegment_(false),
      detected_(false),
      dentries_(0),
      modata_(Teuchos::null),
      nurbsw_(-1.0)
{
  for (std::size_t i = 0; i < coords.size(); ++i)
  {
    uold()[i] = 0.0;
    xspatial()[i] = X()[i];
    dbcdofs_[i] = false;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MORTAR::MortarNode::MortarNode(const MORTAR::MortarNode& old)
    : DRT::Node(old),
      isslave_(old.isslave_),
      istiedslave_(old.istiedslave_),
      isonbound_(old.isonbound_),
      isdbc_(old.isdbc_),
      dofs_(old.dofs_),
      hasproj_(old.hasproj_),
      hassegment_(old.hassegment_),
      detected_(false)
{
  for (int i = 0; i < 3; ++i)
  {
    uold()[i] = old.uold_[i];
    xspatial()[i] = old.xspatial_[i];
    dbcdofs_[i] = false;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MORTAR::MortarNode* MORTAR::MortarNode::Clone() const
{
  auto* newnode = new MORTAR::MortarNode(*this);
  return newnode;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const MORTAR::MortarNode& mrtrnode)
{
  mrtrnode.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::Print(std::ostream& os) const
{
  // Print id and coordinates
  os << "Mortar ";
  DRT::Node::Print(os);

  if (IsSlave())
    os << " Slave  ";
  else
    os << " Master ";

  if (IsOnBound())
    os << " Boundary ";
  else
    os << " Interior ";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class DRT::Node
  DRT::Node::Pack(data);
  // add isslave_
  AddtoPack(data, isslave_);
  // add istiedslave_
  AddtoPack(data, istiedslave_);
  // add isonbound_
  AddtoPack(data, isonbound_);
  // add isonbound_
  AddtoPack(data, isonedge_);
  // add isonbound_
  AddtoPack(data, isoncorner_);
  // add isdbc_
  AddtoPack(data, isdbc_);
  // add dbcdofs_
  AddtoPack(data, dbcdofs_[0]);
  AddtoPack(data, dbcdofs_[1]);
  AddtoPack(data, dbcdofs_[2]);
  // dentries_
  AddtoPack(data, dentries_);
  // add dofs_
  AddtoPack(data, dofs_);
  // add xspatial_
  AddtoPack(data, xspatial_, 3 * sizeof(double));
  // add uold_
  AddtoPack(data, uold_, 3 * sizeof(double));
  // add hasproj_
  AddtoPack(data, hasproj_);
  // add hassegment_
  AddtoPack(data, hassegment_);
  // add nurbsw_
  AddtoPack(data, nurbsw_);

  // add data_
  bool hasdata = (modata_ != Teuchos::null);
  AddtoPack(data, hasdata);
  if (hasdata) modata_->Pack(data);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class DRT::Node
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::Node::Unpack(basedata);
  // isslave_
  isslave_ = ExtractInt(position, data);
  // istiedslave_
  istiedslave_ = ExtractInt(position, data);
  // isonbound_
  isonbound_ = ExtractInt(position, data);
  // isonedge_
  isonedge_ = ExtractInt(position, data);
  // isoncorner_
  isoncorner_ = ExtractInt(position, data);
  // isdbc_
  isdbc_ = ExtractInt(position, data);
  // dbcdofs_
  dbcdofs_[0] = ExtractInt(position, data);
  dbcdofs_[1] = ExtractInt(position, data);
  dbcdofs_[2] = ExtractInt(position, data);
  // dentries_
  ExtractfromPack(position, data, dentries_);
  // dofs_
  ExtractfromPack(position, data, dofs_);
  // xspatial_
  ExtractfromPack(position, data, xspatial_, 3 * sizeof(double));
  // uold_
  ExtractfromPack(position, data, uold_, 3 * sizeof(double));
  // hasproj_
  hasproj_ = ExtractInt(position, data);
  // hassegment_
  hassegment_ = ExtractInt(position, data);
  // nurbsw_
  nurbsw_ = ExtractDouble(position, data);

  // data_
  bool hasdata = ExtractInt(position, data);
  if (hasdata)
  {
    modata_ = Teuchos::rcp(new MORTAR::MortarNodeDataContainer());
    modata_->Unpack(position, data);
  }
  else
  {
    modata_ = Teuchos::null;
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::AddDValue(const int& colnode, const double& val)
{
  // check if this is a master node or slave boundary node
  if (not IsSlave()) dserror("AddDValue: function called for master node %i", Id());
  if (IsOnBound()) dserror("AddDValue: function called for boundary node %i", Id());

  // check if this has been called before
  if ((int)MoData().GetD().size() == 0) MoData().GetD().resize(dentries_);

  // add the pair (col,val) to the given row
  MoData().GetD()[colnode] += val;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::AddDntsValue(const int& colnode, const double& val)
{
  // check if this has been called before
  if ((int)MoData().GetDnts().size() == 0) MoData().GetDnts().resize(dentries_);

  // add the pair (col,val) to the given row
  MoData().GetDnts()[colnode] += val;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::AddDltsValue(const int& colnode, const double& val)
{
  // check if this has been called before
  if ((int)MoData().GetDlts().size() == 0) MoData().GetDlts().resize(dentries_);

  // add the pair (col,val) to the given row
  MoData().GetDlts()[colnode] += val;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::AddDltlValue(const int& colnode, const double& val)
{
  // check if this is a master node or slave boundary node
  if (not IsSlave()) dserror("AddDValue: function called for master node %i", Id());
  if (not IsOnEdge()) dserror("function called for non-edge node %i", Id());

  // check if this has been called before
  if (static_cast<int>(MoData().GetDltl().size()) == 0) MoData().GetDltl().resize(dentries_);

  // add the pair (col,val) to the given row
  MoData().GetDltl()[colnode] += val;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::AddMValue(const int& colnode, const double& val)
{
  // check if this is a master node or slave boundary node
  if (not IsSlave()) dserror("AddMValue: function called for master node %i", Id());
  if (IsOnBoundorCE()) dserror("AddMValue: function called for boundary node %i", Id());

  // add the pair (col,val) to the given row
  MoData().GetM()[colnode] += val;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::AddMntsValue(const int& colnode, const double& val)
{
  // add the pair (col,val) to the given row
  MoData().GetMnts()[colnode] += val;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::AddMltsValue(const int& colnode, const double& val)
{
  // add the pair (col,val) to the given row
  MoData().GetMlts()[colnode] += val;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::AddMltlValue(const int& colnode, const double& val)
{
  // check if this is a master node or slave boundary node
  if (not IsSlave()) dserror("AddMValue: function called for master node %i", Id());
  if (not IsOnEdge()) dserror("function called for non-edge node %i", Id());

  // add the pair (col,val) to the given row
  MoData().GetMltl()[colnode] += val;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::AddMmodValue(const int& colnode, const double& val)
{
  // check if this is a master node or slave boundary node
  if (not IsSlave()) dserror("AddMmodValue: function called for master node %i", Id());
  if (IsOnBound()) dserror("AddMmodValue: function called for boundary node %i", Id());

  // add the pair (col,val) to the given row
  MoData().GetMmod()[colnode] += val;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::InitializeDataContainer()
{
  // get maximum size of nodal D-entries
  dentries_ = 0;
  std::set<int> sIdCheck;
  std::pair<std::set<int>::iterator, bool> check;
  for (int i = 0; i < NumElement(); ++i)
  {
    const int* snodeIds = Elements()[i]->NodeIds();
    for (int j = 0; j < Elements()[i]->NumNode(); ++j)
    {
      check = sIdCheck.insert(snodeIds[j]);
      if (check.second) dentries_ += Elements()[i]->NumDofPerNode(*(Elements()[i]->Nodes()[j]));
    }
  }

  // only initialize if not yet done
  if (modata_ == Teuchos::null) modata_ = Teuchos::rcp(new MORTAR::MortarNodeDataContainer());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::ResetDataContainer()
{
  // reset to Teuchos::null
  modata_ = Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MORTAR::MortarNode::BuildAveragedNormal()
{
  // reset normal and tangents when this method is called
  for (int j = 0; j < 3; ++j) MoData().n()[j] = 0.0;

  int nseg = NumElement();
  DRT::Element** adjeles = Elements();

  // we need to store some stuff here
  //**********************************************************************
  // elens(0,i): x-coord of element normal
  // elens(1,i): y-coord of element normal
  // elens(2,i): z-coord of element normal
  // elens(3,i): id of adjacent element i
  // elens(4,i): length of element normal
  // elens(5,i): length/area of element itself
  //**********************************************************************
  CORE::LINALG::SerialDenseMatrix elens(6, nseg);

  // loop over all adjacent elements
  for (int i = 0; i < nseg; ++i)
  {
    auto* adjmrtrele = dynamic_cast<MortarElement*>(adjeles[i]);

    // build element normal at current node
    // (we have to pass in the index i to be able to store the
    // normal and other information at the right place in elens)
    adjmrtrele->BuildNormalAtNode(Id(), i, elens);

    // add (weighted) element normal to nodal normal n
    for (int j = 0; j < 3; ++j) MoData().n()[j] += elens(j, i) / elens(4, i);
  }

  // create unit normal vector
  double length = sqrt(MoData().n()[0] * MoData().n()[0] + MoData().n()[1] * MoData().n()[1] +
                       MoData().n()[2] * MoData().n()[2]);
  if (length == 0.0)
    dserror("Nodal normal length 0, node ID %i", Id());
  else
    for (int j = 0; j < 3; ++j) MoData().n()[j] /= length;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MORTAR::MortarNode* MORTAR::MortarNode::FindClosestNode(
    const Teuchos::RCP<DRT::Discretization> intdis, const Teuchos::RCP<Epetra_Map> nodesearchmap,
    double& mindist)
{
  MortarNode* closestnode = nullptr;

  // loop over all nodes of the DRT::Discretization that are
  // included in the given Epetra_Map ("brute force" search)
  for (int i = 0; i < nodesearchmap->NumMyElements(); ++i)
  {
    int gid = nodesearchmap->GID(i);
    DRT::Node* node = intdis->gNode(gid);
    if (!node) dserror("FindClosestNode: Cannot find node with gid %", gid);
    auto* mrtrnode = dynamic_cast<MortarNode*>(node);

    // build distance between the two nodes
    double dist = 0.0;
    const double* p1 = xspatial();
    const double* p2 = mrtrnode->xspatial();

    for (int j = 0; j < 3; ++j) dist += (p1[j] - p2[j]) * (p1[j] - p2[j]);
    dist = sqrt(dist);

    // new closest node found, update
    if (dist <= mindist)
    {
      mindist = dist;
      closestnode = mrtrnode;
    }
  }

  if (!closestnode) dserror("FindClosestNode: No closest node found at all!");

  return closestnode;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool MORTAR::MortarNode::CheckMeshDistortion(double& relocation, const double& limit)
{
  // initialize return parameter
  bool ok = true;

  // loop over all adjacent elements of this node
  for (int i = 0; i < NumElement(); ++i)
  {
    // get the current element
    DRT::Element* ele = Elements()[i];
    if (!ele) dserror("Cannot find element with lid %", i);
    auto* mrtrele = dynamic_cast<MortarElement*>(ele);

    // minimal edge size of the current element
    const double minedgesize = mrtrele->MinEdgeSize();

    // check whether relocation is not too large
    if (relocation > limit * minedgesize)
    {
      // print information to screen
      std::cout << "\n*****************WARNING***********************" << '\n';
      std::cout << "Checking distortion for CNode:     " << Id() << '\n';
      std::cout << "Relocation distance:               " << relocation << '\n';
      std::cout << "AdjEle: " << mrtrele->Id() << "\tLimit*MinEdgeSize: " << limit * minedgesize
                << '\n';
      std::cout << "*****************WARNING***********************" << '\n';

      // set return parameter and stop
      ok = false;
      break;
    }
  }

  return ok;
}
