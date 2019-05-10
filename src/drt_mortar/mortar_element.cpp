/*----------------------------------------------------------------------*/
/*!
\file mortar_element.cpp

\brief A mortar coupling element

\level 2

\maintainer Matthias Mayr

*/
/*----------------------------------------------------------------------*/

#include "mortar_element.H"
#include "mortar_node.H"
#include "mortar_calc_utils.H"

#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"

#include "../drt_contact/contact_nitsche_utils.H"

#include "../drt_mat/material.H"

#include "../drt_so3/so_surface.H"
#include <Teuchos_RCP.hpp>



MORTAR::MortarElementType MORTAR::MortarElementType::instance_;

MORTAR::MortarElementType& MORTAR::MortarElementType::Instance() { return instance_; }

DRT::ParObject* MORTAR::MortarElementType::Create(const std::vector<char>& data)
{
  MORTAR::MortarElement* ele =
      new MORTAR::MortarElement(0, 0, DRT::Element::dis_none, 0, NULL, false);
  ele->Unpack(data);
  return ele;
}


Teuchos::RCP<DRT::Element> MORTAR::MortarElementType::Create(const int id, const int owner)
{
  // return Teuchos::rcp( new MortarElement( id, owner ) );
  return Teuchos::null;
}


void MORTAR::MortarElementType::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
}

void MORTAR::MortarElementType::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 12/10|
 *----------------------------------------------------------------------*/
MORTAR::MortarEleDataContainer::MortarEleDataContainer()
{
  // initialize area
  Area() = 0.0;
  dualshapecoeff_ = Teuchos::null;
  derivdualshapecoeff_ = Teuchos::null;
  trafocoeff_ = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            popp 12/10|
 *----------------------------------------------------------------------*/
void MORTAR::MortarEleDataContainer::Pack(DRT::PackBuffer& data) const
{
  // add area_
  DRT::ParObject::AddtoPack(data, area_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 02/10|
 *----------------------------------------------------------------------*/
void MORTAR::MortarEleDataContainer::Unpack(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  // area_
  DRT::ParObject::ExtractfromPack(position, data, area_);

  dualshapecoeff_ = Teuchos::null;
  derivdualshapecoeff_ = Teuchos::null;
  return;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
MORTAR::MortarElement::MortarElement(int id, int owner,
    const DRT::Element::DiscretizationType& shape, const int numnode, const int* nodeids,
    const bool isslave, bool isnurbs)
    : DRT::FaceElement(id, owner),
      shape_(shape),
      isslave_(isslave),
      attached_(false),
      nurbs_(isnurbs),
      normalfac_(1.0),    // normal factor for nurbs
      zero_sized_(false)  // information for nurbs integration
{
  SetNodeIds(numnode, nodeids);
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (protected)                                   kronbichler 03/15|
 *----------------------------------------------------------------------*/
MORTAR::MortarElement::MortarElement(int id, int owner)
    : DRT::FaceElement(id, owner),
      shape_(DRT::Element::dis_none),
      isslave_(false),
      attached_(false),
      nurbs_(false),
      normalfac_(1.0),    // normal factor for nurbs
      zero_sized_(false)  // information for nurbs integration
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 10/07|
 *----------------------------------------------------------------------*/
MORTAR::MortarElement::MortarElement(const MORTAR::MortarElement& old)
    : DRT::FaceElement(old), shape_(old.shape_), isslave_(old.isslave_)
{
  // not yet used and thus not necessarily consistent
  dserror("ERROR: MortarElement copy-ctor not yet implemented");

  return;
}

/*----------------------------------------------------------------------*
 |  clone-ctor (public)                                      mwgee 10/07|
 *----------------------------------------------------------------------*/
DRT::Element* MORTAR::MortarElement::Clone() const
{
  MORTAR::MortarElement* newele = new MORTAR::MortarElement(*this);
  return newele;
}


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const MORTAR::MortarElement& element)
{
  element.Print(os);
  return os;
}


/*----------------------------------------------------------------------*
 |  print element (public)                                   mwgee 10/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElement::Print(std::ostream& os) const
{
  os << "Mortar Element ";
  DRT::Element::Print(os);
  if (isslave_)
    os << " Slave  ";
  else
    os << " Master ";

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElement::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class DRT::FaceElement
  DRT::FaceElement::Pack(data);
  // add shape_
  AddtoPack(data, shape_);
  // add isslave_
  AddtoPack(data, isslave_);
  // add nurbs_
  AddtoPack(data, nurbs_);

  // for nurbs:
  if (nurbs_)
  {
    // add normalfac
    AddtoPack(data, normalfac_);
    // add zero_sized_
    AddtoPack(data, zero_sized_);
    // knots
    int nr = mortarknots_.size();
    DRT::ParObject::AddtoPack(data, nr);
    if (nr != 0)
    {
      for (int i = 0; i < nr; i++) DRT::ParObject::AddtoPack(data, (mortarknots_[i]));
    }
  }

  // add modata_
  bool hasdata = (modata_ != Teuchos::null);
  AddtoPack(data, hasdata);
  if (hasdata) modata_->Pack(data);

  // add physicaltype
  AddtoPack(data, static_cast<int>(physicaltype_));

  // mesh size
  AddtoPack(data, traceHE_);
  AddtoPack(data, traceHCond_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElement::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class DRT::FaceElement
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::FaceElement::Unpack(basedata);
  // shape_
  shape_ = static_cast<DiscretizationType>(ExtractInt(position, data));
  // isslave_
  isslave_ = ExtractInt(position, data);
  // nurbs_
  nurbs_ = ExtractInt(position, data);

  // for nurbs:
  if (nurbs_)
  {
    // normalfac_
    normalfac_ = ExtractDouble(position, data);
    // zero_sized_
    zero_sized_ = ExtractInt(position, data);
    // knots
    int nr;
    DRT::ParObject::ExtractfromPack(position, data, nr);

    if (nr != 0)
    {
      mortarknots_.resize(nr);
      for (int i = 0; i < nr; i++) DRT::ParObject::ExtractfromPack(position, data, mortarknots_[i]);
    }
  }

  // modata_
  bool hasdata = ExtractInt(position, data);
  if (hasdata)
  {
    modata_ = Teuchos::rcp(new MORTAR::MortarEleDataContainer());
    modata_->Unpack(position, data);
  }
  else
  {
    modata_ = Teuchos::null;
  }

  // physical type
  physicaltype_ = (PhysicalType)(ExtractInt(position, data));

  // mesh size
  traceHE_ = ExtractDouble(position, data);
  traceHCond_ = ExtractDouble(position, data);

  if (position != data.size())
    dserror(
        "Mismatch in size of available data (size %d) vs. position pointer "
        "of read data (size %d)",
        (int)data.size(), position);

  return;
}

/*----------------------------------------------------------------------*
 |  number of dofs per node (public)                         mwgee 10/07|
 *----------------------------------------------------------------------*/
int MORTAR::MortarElement::NumDofPerNode(const DRT::Node& node) const
{
  const MORTAR::MortarNode* mnode = dynamic_cast<const MORTAR::MortarNode*>(&node);
  if (!mnode) dserror("Node is not a MortarNode");
  return mnode->NumDof();
}

/*----------------------------------------------------------------------*
 |  evaluate element (public)                                mwgee 10/07|
 *----------------------------------------------------------------------*/
int MORTAR::MortarElement::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  dserror("MORTAR::MortarElement::Evaluate not implemented!");
  return -1;
}

/*----------------------------------------------------------------------*
 |  Get local coordinates for local node id                   popp 12/07|
 *----------------------------------------------------------------------*/
bool MORTAR::MortarElement::LocalCoordinatesOfNode(int lid, double* xi) const
{
  // 2D linear case (2noded line element)
  // 2D quadratic case (3noded line element)
  switch (Shape())
  {
    case line2:
    case line3:
    {
      switch (lid)
      {
        case 0:
          xi[0] = -1.0;
          break;
        case 1:
          xi[0] = 1.0;
          break;
        case 2:
          xi[0] = 0.0;
          break;
        default:
          dserror(
              "ERROR: LocalCoordinatesOfNode: Node number % in segment % out of range", lid, Id());
      }
      // we are in the 2D case here!
      xi[1] = 0.0;

      break;
    }

    // 3D linear case (2noded triangular element)
    // 3D quadratic case (3noded triangular element)
    case tri3:
    case tri6:
    {
      switch (lid)
      {
        case 0:
        {
          xi[0] = 0.0;
          xi[1] = 0.0;
          break;
        }
        case 1:
        {
          xi[0] = 1.0;
          xi[1] = 0.0;
          break;
        }
        case 2:
        {
          xi[0] = 0.0;
          xi[1] = 1.0;
          break;
        }
        case 3:
        {
          xi[0] = 0.5;
          xi[1] = 0.0;
          break;
        }
        case 4:
        {
          xi[0] = 0.5;
          xi[1] = 0.5;
          break;
        }
        case 5:
        {
          xi[0] = 0.0;
          xi[1] = 0.5;
          break;
        }
        default:
          dserror("ERROR: LocCoordsOfNode: Node number % in segment % out of range", lid, Id());
      }

      break;
    }

    // 3D bilinear case (4noded quadrilateral element)
    // 3D serendipity case (8noded quadrilateral element)
    // 3D biquadratic case (9noded quadrilateral element)
    case quad4:
    case quad8:
    case quad9:
    {
      switch (lid)
      {
        case 0:
        {
          xi[0] = -1.0;
          xi[1] = -1.0;
          break;
        }
        case 1:
        {
          xi[0] = 1.0;
          xi[1] = -1.0;
          break;
        }
        case 2:
        {
          xi[0] = 1.0;
          xi[1] = 1.0;
          break;
        }
        case 3:
        {
          xi[0] = -1.0;
          xi[1] = 1.0;
          break;
        }
        case 4:
        {
          xi[0] = 0.0;
          xi[1] = -1.0;
          break;
        }
        case 5:
        {
          xi[0] = 1.0;
          xi[1] = 0.0;
          break;
        }
        case 6:
        {
          xi[0] = 0.0;
          xi[1] = 1.0;
          break;
        }
        case 7:
        {
          xi[0] = -1.0;
          xi[1] = 0.0;
          break;
        }
        case 8:
        {
          xi[0] = 0.0;
          xi[1] = 0.0;
          break;
        }
        default:
          dserror("ERROR: LocCoordsOfNode: Node number % in segment % out of range", lid, Id());
      }

      break;
    }

    //==================================================
    //                     NURBS
    case nurbs2:
    {
      if (lid == 0)
        xi[0] = -1.0;
      else if (lid == 1)
        xi[0] = 1.0;
      else
        dserror(
            "ERROR: LocalCoordinatesOfNode: Node number % in segment % out of range", lid, Id());

      // we are in the 2D case here!
      xi[1] = 0.0;

      break;
    }
    case nurbs3:
    {
      if (lid == 0)
        xi[0] = -1.0;
      else if (lid == 1)
        xi[0] = 0.0;
      else if (lid == 2)
        xi[0] = 1.0;
      else
        dserror(
            "ERROR: LocalCoordinatesOfNode: Node number % in segment % out of range", lid, Id());

      // we are in the 2D case here!
      xi[1] = 0.0;

      break;
    }
    case nurbs9:
    {
      switch (lid)
      {
        case 0:
        {
          xi[0] = -1.0;
          xi[1] = -1.0;
          break;
        }
        case 1:
        {
          xi[0] = 0.0;
          xi[1] = -1.0;
          break;
        }
        case 2:
        {
          xi[0] = 1.0;
          xi[1] = -1.0;
          break;
        }
        case 3:
        {
          xi[0] = -1.0;
          xi[1] = 0.0;
          break;
        }
        case 4:
        {
          xi[0] = 0.0;
          xi[1] = 0.0;
          break;
        }
        case 5:
        {
          xi[0] = 1.0;
          xi[1] = 0.0;
          break;
        }
        case 6:
        {
          xi[0] = -1.0;
          xi[1] = 1.0;
          break;
        }
        case 7:
        {
          xi[0] = 0.0;
          xi[1] = 1.0;
          break;
        }
        case 8:
        {
          xi[0] = 1.0;
          xi[1] = 1.0;
          break;
        }
        default:
          dserror("ERROR: LocCoordsOfNode: Node number % in segment % out of range", lid, Id());
      }

      break;
    }
    // unknown case
    default:
      dserror("ERROR: LocalCoordinatesOfNode called for unknown element type");
      exit(EXIT_FAILURE);
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Get local numbering for global node id                    popp 12/07|
 *----------------------------------------------------------------------*/
int MORTAR::MortarElement::GetLocalNodeId(int nid) const
{
  int lid = -1;

  // look for global ID nid in element's nodes
  for (int i = 0; i < NumNode(); ++i)
    if (NodeIds()[i] == nid)
    {
      lid = i;
      break;
    }

  if (lid < 0) dserror("ERROR: Cannot find node % in segment %", nid, Id());

  return lid;
}

/*----------------------------------------------------------------------*
 |  Build element normal at node                              popp 12/07|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElement::BuildNormalAtNode(int nid, int& i, Epetra_SerialDenseMatrix& elens)
{
  // find this node in my list of nodes and get local numbering
  int lid = GetLocalNodeId(nid);

  // get local coordinates for this node
  double xi[2];
  LocalCoordinatesOfNode(lid, xi);

  // build an outward unit normal at xi and return it
  ComputeNormalAtXi(xi, i, elens);

  return;
}

/*----------------------------------------------------------------------*
 |  Compute element normal at loc. coord. xi                  popp 09/08|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElement::ComputeNormalAtXi(double* xi, int& i, Epetra_SerialDenseMatrix& elens)
{
  // empty local basis vectors
  double gxi[3];
  double geta[3];

  // metrics routine gives local basis vectors
  Metrics(xi, gxi, geta);

  // n is cross product of gxi and geta
  elens(0, i) = (gxi[1] * geta[2] - gxi[2] * geta[1]) * NormalFac();
  elens(1, i) = (gxi[2] * geta[0] - gxi[0] * geta[2]) * NormalFac();
  elens(2, i) = (gxi[0] * geta[1] - gxi[1] * geta[0]) * NormalFac();

  // store length of normal and other information into elens
  elens(4, i) =
      sqrt(elens(0, i) * elens(0, i) + elens(1, i) * elens(1, i) + elens(2, i) * elens(2, i));
  if (elens(4, i) < 1e-12) dserror("ERROR: ComputeNormalAtXi gives normal of length 0!");
  elens(3, i) = Id();
  elens(5, i) = MoData().Area();

  return;
}

/*----------------------------------------------------------------------*
 |  Compute element normal at loc. coord. xi                  popp 11/08|
 *----------------------------------------------------------------------*/
double MORTAR::MortarElement::ComputeUnitNormalAtXi(double* xi, double* n)
{
  // check input
  if (!xi) dserror("ERROR: ComputeUnitNormalAtXi called with xi=NULL");
  if (!n) dserror("ERROR: ComputeUnitNormalAtXi called with n=NULL");

  // empty local basis vectors
  double gxi[3];
  double geta[3];

  // metrics routine gives local basis vectors
  Metrics(xi, gxi, geta);

  // n is cross product of gxi and geta
  n[0] = (gxi[1] * geta[2] - gxi[2] * geta[1]) * NormalFac();
  n[1] = (gxi[2] * geta[0] - gxi[0] * geta[2]) * NormalFac();
  n[2] = (gxi[0] * geta[1] - gxi[1] * geta[0]) * NormalFac();

  // build unit normal
  const double length = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  if (length < 1e-12) dserror("ERROR: Normal of length zero!");
  for (int i = 0; i < 3; ++i) n[i] /= length;

  return length;
}


/*----------------------------------------------------------------------*
 |  Compute nodal averaged normal at xi                      farah 06/16|
 *----------------------------------------------------------------------*/
double MORTAR::MortarElement::ComputeAveragedUnitNormalAtXi(double* xi, double* n)
{
  // check input
  if (!xi) dserror("ERROR: ComputeUnitNormalAtXi called with xi=NULL");
  if (!n) dserror("ERROR: ComputeUnitNormalAtXi called with n=NULL");

  int nnodes = NumPoint();
  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes, 2, true);

  // get shape function values and derivatives at xi
  EvaluateShape(xi, val, deriv, nnodes, false);

  // initialize n
  n[0] = 0.0;
  n[1] = 0.0;
  n[2] = 0.0;

  // loop over all nodes of this element
  for (int i = 0; i < NumNode(); ++i)
  {
    MortarNode* mymrtrnode = dynamic_cast<MortarNode*>(Nodes()[i]);
    n[0] = val[i] * mymrtrnode->MoData().n()[0];
    n[1] = val[i] * mymrtrnode->MoData().n()[1];
    n[2] = val[i] * mymrtrnode->MoData().n()[2];
  }

  const double length = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  if (length < 1e-12) dserror("ERROR: Normal of length zero!");
  for (int i = 0; i < 3; ++i) n[i] /= length;

  return length;
}

/*----------------------------------------------------------------------*
 |  Compute unit normal derivative at loc. coord. xi          popp 03/09|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElement::DerivUnitNormalAtXi(
    double* xi, std::vector<GEN::pairedvector<int, double>>& derivn)
{
  // initialize variables
  const int nnodes = NumNode();
  DRT::Node** mynodes = Nodes();
  if (!mynodes) dserror("ERROR: DerivUnitNormalAtXi: Null pointer!");

  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes, 2, true);

  double gxi[3];
  double geta[3];

  // get shape function values and derivatives at xi
  EvaluateShape(xi, val, deriv, nnodes);

  // get local element basis vectors
  Metrics(xi, gxi, geta);

  // n is cross product of gxi and geta
  double n[3] = {0.0, 0.0, 0.0};
  n[0] = gxi[1] * geta[2] - gxi[2] * geta[1] * normalfac_;
  n[1] = gxi[2] * geta[0] - gxi[0] * geta[2] * normalfac_;
  n[2] = gxi[0] * geta[1] - gxi[1] * geta[0] * normalfac_;

  // build unit normal
  const double length = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  if (length < 1e-12) dserror("ERROR: Normal of length zero!");
  for (int i = 0; i < 3; ++i) n[i] /= length;

  // check if this mortar ele is an IntEle
  std::vector<std::vector<GEN::pairedvector<int, double>>> nodelin(0);
  NodeLinearization(nodelin);

  int nderiv = nnodes * 3;
  // to be safe if it is a IntEle for a nurbs9
  if (Shape() == DRT::Element::quad4) nderiv = 9 * 3;

  // resize derivn
  derivn.resize(3, nderiv);

  // non-unit normal derivative
  std::vector<GEN::pairedvector<int, double>> derivnnu(
      3, nderiv);  // assume that each node has 3 dofs...
  typedef GEN::pairedvector<int, double>::const_iterator CI;

  // now the derivative
  for (int n = 0; n < nnodes; ++n)
  {
    MortarNode* mymrtrnode = dynamic_cast<MortarNode*>(mynodes[n]);
    if (!mymrtrnode) dserror("ERROR: DerivUnitNormalAtXi: Null pointer!");
    int ndof = mymrtrnode->NumDof();

    // derivative weighting matrix for current node
    LINALG::Matrix<3, 3> F;
    F(0, 0) = 0.0;
    F(1, 1) = 0.0;
    F(2, 2) = 0.0;
    F(0, 1) = geta[2] * deriv(n, 0) - gxi[2] * deriv(n, 1);
    F(0, 2) = gxi[1] * deriv(n, 1) - geta[1] * deriv(n, 0);
    F(1, 0) = gxi[2] * deriv(n, 1) - geta[2] * deriv(n, 0);
    F(1, 2) = geta[0] * deriv(n, 0) - gxi[0] * deriv(n, 1);
    F(2, 0) = geta[1] * deriv(n, 0) - gxi[1] * deriv(n, 1);
    F(2, 1) = gxi[0] * deriv(n, 1) - geta[0] * deriv(n, 0);

    // create directional derivatives
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < ndof; ++k)
        for (CI p = nodelin[n][k].begin(); p != nodelin[n][k].end(); ++p)
          (derivnnu[j])[p->first] += F(j, k) * p->second;
  }

  const double ll = length * length;
  const double linv = 1.0 / length;
  const double lllinv = 1.0 / (length * length * length);
  const double sxsx = n[0] * n[0] * ll;
  const double sxsy = n[0] * n[1] * ll;
  const double sxsz = n[0] * n[2] * ll;
  const double sysy = n[1] * n[1] * ll;
  const double sysz = n[1] * n[2] * ll;
  const double szsz = n[2] * n[2] * ll;

  for (CI p = derivnnu[0].begin(); p != derivnnu[0].end(); ++p)
  {
    derivn[0][p->first] += linv * (p->second) * normalfac_;
    derivn[0][p->first] -= lllinv * sxsx * (p->second) * normalfac_;
    derivn[1][p->first] -= lllinv * sxsy * (p->second) * normalfac_;
    derivn[2][p->first] -= lllinv * sxsz * (p->second) * normalfac_;
  }

  for (CI p = derivnnu[1].begin(); p != derivnnu[1].end(); ++p)
  {
    derivn[1][p->first] += linv * (p->second) * normalfac_;
    derivn[1][p->first] -= lllinv * sysy * (p->second) * normalfac_;
    derivn[0][p->first] -= lllinv * sxsy * (p->second) * normalfac_;
    derivn[2][p->first] -= lllinv * sysz * (p->second) * normalfac_;
  }

  for (CI p = derivnnu[2].begin(); p != derivnnu[2].end(); ++p)
  {
    derivn[2][p->first] += linv * (p->second) * normalfac_;
    derivn[2][p->first] -= lllinv * szsz * (p->second) * normalfac_;
    derivn[0][p->first] -= lllinv * sxsz * (p->second) * normalfac_;
    derivn[1][p->first] -= lllinv * sysz * (p->second) * normalfac_;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Get nodal coordinates of the element                      popp 01/08|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElement::GetNodalCoords(LINALG::SerialDenseMatrix& coord)
{
  const int nnodes = NumPoint();
  DRT::Node** mynodes = Points();
  if (!mynodes) dserror("ERROR: GetNodalCoords: Null pointer!");
  if (coord.M() != 3 || coord.N() != nnodes) dserror("ERROR: GetNodalCoords: Dimensions!");

  for (int i = 0; i < nnodes; ++i)
  {
    MortarNode* mymrtrnode = dynamic_cast<MortarNode*>(mynodes[i]);
    if (!mymrtrnode) dserror("ERROR: GetNodalCoords: Null pointer!");

    const double* x = mymrtrnode->xspatial();
    std::copy(x, x + 3, &coord(0, i));
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Get old nodal coordinates of the element               gitterle 08/10|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElement::GetNodalCoordsOld(LINALG::SerialDenseMatrix& coord, bool isinit)
{
  const int nnodes = NumPoint();
  DRT::Node** mynodes = Points();
  if (!mynodes) dserror("ERROR: GetNodalCoordsOld: Null pointer!");
  if (coord.M() != 3 || coord.N() != nnodes) dserror("ERROR: GetNodalCoordsOld: Dimensions!");

  for (int i = 0; i < nnodes; ++i)
  {
    MortarNode* mymrtrnode = dynamic_cast<MortarNode*>(mynodes[i]);
    if (!mymrtrnode) dserror("ERROR: GetNodalCoordsOld: Null pointer!");

    coord(0, i) = mymrtrnode->X()[0] + mymrtrnode->uold()[0];
    coord(1, i) = mymrtrnode->X()[1] + mymrtrnode->uold()[1];
    coord(2, i) = mymrtrnode->X()[2] + mymrtrnode->uold()[2];
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Get lagrange multipliers of the element                gitterle 08/10|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElement::GetNodalLagMult(LINALG::SerialDenseMatrix& lagmult, bool isinit)
{
  int nnodes = NumNode();
  DRT::Node** mynodes = Nodes();
  if (!mynodes) dserror("ERROR: GetNodalLagMult: Null pointer!");
  if (lagmult.M() != 3 || lagmult.N() != nnodes) dserror("ERROR: GetNodalLagMult: Dimensions!");

  for (int i = 0; i < nnodes; ++i)
  {
    MortarNode* mymrtrnode = dynamic_cast<MortarNode*>(mynodes[i]);
    if (!mymrtrnode) dserror("ERROR: GetNodalCoords: Null pointer!");

    lagmult(0, i) = mymrtrnode->MoData().lm()[0];
    lagmult(1, i) = mymrtrnode->MoData().lm()[1];
    lagmult(2, i) = mymrtrnode->MoData().lm()[2];
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate element metrics (local basis vectors)            popp 08/08|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElement::Metrics(const double* xi, double* gxi, double* geta)
{
  std::fill(gxi, gxi + 3, 0.0);
  std::fill(geta, geta + 3, 0.0);

  int nnodes = NumPoint();

  int dim = 0;
  DRT::Element::DiscretizationType dt = Shape();
  switch (dt)
  {
    case line2:
    case line3:
    case nurbs2:
    case nurbs3:
    {
      dim = 2;
      break;
    }
    case tri3:
    case quad4:
    case tri6:
    case quad8:
    case quad9:
    case nurbs4:
    case nurbs8:
    case nurbs9:
    {
      dim = 3;
      break;
    }
    default:
      dserror("ERROR: Metrics called for unknown element type");
      exit(EXIT_FAILURE);
  }

  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes, 2, true);

  // get shape function values and derivatives at xi
  EvaluateShape(xi, val, deriv, nnodes, false);

  // get coordinates of element nodes
  LINALG::SerialDenseMatrix coord(3, nnodes);
  GetNodalCoords(coord);

  // build basis vectors gxi and geta
  for (int i = 0; i < nnodes; ++i)
  {
    // first local basis vector
    gxi[0] += deriv(i, 0) * coord(0, i);
    gxi[1] += deriv(i, 0) * coord(1, i);
    gxi[2] += deriv(i, 0) * coord(2, i);

    // second local basis vector
    geta[0] += deriv(i, 1) * coord(0, i);
    geta[1] += deriv(i, 1) * coord(1, i);
    geta[2] += deriv(i, 1) * coord(2, i);
  }

  // reset geta to (0,0,1) in 2D case
  if (dim == 2)
  {
    geta[0] = 0.0;
    geta[1] = 0.0;
    geta[2] = 1.0;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate Jacobian determinant                             popp 12/07|
 *----------------------------------------------------------------------*/
double MORTAR::MortarElement::Jacobian(const double* xi)
{
  double jac = 0.0;
  double gxi[3];
  double geta[3];
  DRT::Element::DiscretizationType dt = Shape();

  // 2D linear case (2noded line element)
  if (dt == line2) jac = MoData().Area() * 0.5;

  // 3D linear case (3noded triangular element)
  else if (dt == tri3)
    jac = MoData().Area() * 2.0;

  // 2D quadratic case (3noded line element)
  // 3D bilinear case (4noded quadrilateral element)
  // 3D quadratic case (6noded triangular element)
  // 3D serendipity case (8noded quadrilateral element)
  // 3D biquadratic case (9noded quadrilateral element)
  else if (dt == line3 || dt == quad4 || dt == tri6 || dt == quad8 || dt == quad9 || dt == nurbs2 ||
           dt == nurbs3 || dt == nurbs4 || dt == nurbs8 || dt == nurbs9)
  {
    // metrics routine gives local basis vectors
    Metrics(xi, gxi, geta);

    // cross product of gxi and geta
    double cross[3] = {0.0, 0.0, 0.0};
    cross[0] = gxi[1] * geta[2] - gxi[2] * geta[1];
    cross[1] = gxi[2] * geta[0] - gxi[0] * geta[2];
    cross[2] = gxi[0] * geta[1] - gxi[1] * geta[0];
    jac = sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);
  }

  // unknown case
  else
    dserror("ERROR: Jacobian called for unknown element type!");

  return jac;
}

/*----------------------------------------------------------------------*
 |  Evaluate directional deriv. of Jacobian det.              popp 05/08|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElement::DerivJacobian(
    const double* xi, GEN::pairedvector<int, double>& derivjac)
{
  // get element nodes
  int nnodes = NumNode();

  DRT::Node** mynodes = NULL;  // Nodes();
  mynodes = Nodes();

  if (!mynodes) dserror("ERROR: DerivJacobian: Null pointer!");

  // the inverse Jacobian
  double jacinv = 0.0;
  double gxi[3];
  double geta[3];

  // evaluate shape functions
  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes, 2, true);
  EvaluateShape(xi, val, deriv, nnodes, false);

  // metrics routine gives local basis vectors
  Metrics(xi, gxi, geta);

  // cross product of gxi and geta
  double cross[3] = {0.0, 0.0, 0.0};
  cross[0] = gxi[1] * geta[2] - gxi[2] * geta[1];
  cross[1] = gxi[2] * geta[0] - gxi[0] * geta[2];
  cross[2] = gxi[0] * geta[1] - gxi[1] * geta[0];

  DRT::Element::DiscretizationType dt = Shape();

  // 2D linear case (2noded line element)
  switch (dt)
  {
    // 3D linear case (3noded triangular element)
    case tri3:
    {
      jacinv = 1.0 / (MoData().Area() * 2.0);
      break;
    }
    // default 2-D case
    case line2:
    {
      jacinv = 2.0 / MoData().Area();
      break;
    }
    // 2D quadratic case (3noded line element)
    // 3D bilinear case (4noded quadrilateral element)
    // 3D quadratic case (6noded triangular element)
    // 3D serendipity case (8noded quadrilateral element)
    // 3D biquadratic case (9noded quadrilateral element)
    /* no break (upper case) */
    case line3:
    case quad4:
    case tri6:
    case quad8:
    case quad9:
    case nurbs2:
    case nurbs3:
    case nurbs4:
    case nurbs8:
    case nurbs9:
    {
      jacinv = 1.0 / sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);
      break;
    }
    default:
      dserror("ERROR: Jac. derivative not implemented for this type of CoElement");
      exit(EXIT_FAILURE);
  }

  // *********************************************************************
  // compute Jacobian derivative
  // *********************************************************************
  // (loop over all nodes and over all nodal dofs to capture all
  // potential dependencies of the Jacobian. Note that here we only
  // need to compute the DIRECT derivative of Lin(J), as the current
  // GP coordinate does not change! The derivative DJacDXi is done in
  // a special function (see above)!
  // *********************************************************************
  for (int i = 0; i < nnodes; ++i)
  {
    MORTAR::MortarNode* mymrtrnode = dynamic_cast<MORTAR::MortarNode*>(mynodes[i]);
    if (!mymrtrnode) dserror("ERROR: DerivJacobian: Null pointer!");

    derivjac[mymrtrnode->Dofs()[0]] +=
        jacinv * (cross[2] * geta[1] - cross[1] * geta[2]) * deriv(i, 0);
    derivjac[mymrtrnode->Dofs()[0]] +=
        jacinv * (cross[1] * gxi[2] - cross[2] * gxi[1]) * deriv(i, 1);
    derivjac[mymrtrnode->Dofs()[1]] +=
        jacinv * (cross[0] * geta[2] - cross[2] * geta[0]) * deriv(i, 0);
    derivjac[mymrtrnode->Dofs()[1]] +=
        jacinv * (cross[2] * gxi[0] - cross[0] * gxi[2]) * deriv(i, 1);

    if (mymrtrnode->NumDof() == 3)
    {
      derivjac[mymrtrnode->Dofs()[2]] +=
          jacinv * (cross[1] * geta[0] - cross[0] * geta[1]) * deriv(i, 0);
      derivjac[mymrtrnode->Dofs()[2]] +=
          jacinv * (cross[0] * gxi[1] - cross[1] * gxi[0]) * deriv(i, 1);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Compute length / area of the element                      popp 12/07|
 *----------------------------------------------------------------------*/
double MORTAR::MortarElement::ComputeArea()
{
  double area = 0.0;
  DRT::Element::DiscretizationType dt = Shape();

  // 2D linear case (2noded line element)
  if (dt == line2)
  {
    // no integration necessary (constant Jacobian)
    LINALG::SerialDenseMatrix coord(3, NumPoint());
    GetNodalCoords(coord);

    // build vector between the two nodes
    double tang[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < 3; ++k)
    {
      tang[k] = coord(k, 1) - coord(k, 0);
    }
    area = sqrt(tang[0] * tang[0] + tang[1] * tang[1] + tang[2] * tang[2]);
  }

  // 3D linear case (3noded triangular element)
  else if (dt == tri3)
  {
    // no integration necessary (constant Jacobian)
    LINALG::SerialDenseMatrix coord(3, NumPoint());
    GetNodalCoords(coord);

    // build vectors between the three nodes
    double t1[3] = {0.0, 0.0, 0.0};
    double t2[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < 3; ++k)
    {
      t1[k] = coord(k, 1) - coord(k, 0);
      t2[k] = coord(k, 2) - coord(k, 0);
    }

    // cross product of t1 and t2
    double t1xt2[3] = {0.0, 0.0, 0.0};
    t1xt2[0] = t1[1] * t2[2] - t1[2] * t2[1];
    t1xt2[1] = t1[2] * t2[0] - t1[0] * t2[2];
    t1xt2[2] = t1[0] * t2[1] - t1[1] * t2[0];
    area = 0.5 * sqrt(t1xt2[0] * t1xt2[0] + t1xt2[1] * t1xt2[1] + t1xt2[2] * t1xt2[2]);
  }

  // 2D quadratic case   (3noded line element)
  // 3D bilinear case    (4noded quadrilateral element)
  // 3D quadratic case   (6noded triangular element)
  // 3D serendipity case (8noded quadrilateral element)
  // 3D biquadratic case (9noded quadrilateral element)
  else if (dt == line3 || dt == quad4 || dt == tri6 || dt == quad8 || dt == quad9 || dt == nurbs2 ||
           dt == nurbs3 || dt == nurbs4 || dt == nurbs8 || dt == nurbs9)
  {
    // Gauss quadrature with correct NumGP and Dim
    MORTAR::ElementIntegrator integrator(dt);
    double detg = 0.0;

    // loop over all Gauss points, build Jacobian and compute area
    for (int j = 0; j < integrator.nGP(); ++j)
    {
      double gpc[2] = {integrator.Coordinate(j, 0), integrator.Coordinate(j, 1)};
      detg = Jacobian(gpc);
      area += integrator.Weight(j) * detg;
    }
  }

  // other cases not implemented yet
  else
    dserror("ERROR: Area computation not implemented for this type of MortarElement");

  return area;
}


/*----------------------------------------------------------------------*
 |  Compute length / area of the element                     seitz 09/17|
 *----------------------------------------------------------------------*/
double MORTAR::MortarElement::ComputeAreaDeriv(GEN::pairedvector<int, double>& area_deriv)
{
  double area = 0.0;
  DRT::Element::DiscretizationType dt = Shape();

  // 2D linear case (2noded line element)
  if (dt == line2)
  {
    // no integration necessary (constant Jacobian)
    LINALG::SerialDenseMatrix coord(3, NumPoint());
    GetNodalCoords(coord);

    // build vector between the two nodes
    double tang[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < 3; ++k)
    {
      tang[k] = coord(k, 1) - coord(k, 0);
    }
    area = sqrt(tang[0] * tang[0] + tang[1] * tang[1] + tang[2] * tang[2]);
  }

  // 3D linear case (3noded triangular element)
  else if (dt == tri3)
  {
    // no integration necessary (constant Jacobian)
    LINALG::SerialDenseMatrix coord(3, NumPoint());
    GetNodalCoords(coord);

    // build vectors between the three nodes
    double t1[3] = {0.0, 0.0, 0.0};
    double t2[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < 3; ++k)
    {
      t1[k] = coord(k, 1) - coord(k, 0);
      t2[k] = coord(k, 2) - coord(k, 0);
    }

    // cross product of t1 and t2
    double t1xt2[3] = {0.0, 0.0, 0.0};
    t1xt2[0] = t1[1] * t2[2] - t1[2] * t2[1];
    t1xt2[1] = t1[2] * t2[0] - t1[0] * t2[2];
    t1xt2[2] = t1[0] * t2[1] - t1[1] * t2[0];
    area = 0.5 * sqrt(t1xt2[0] * t1xt2[0] + t1xt2[1] * t1xt2[1] + t1xt2[2] * t1xt2[2]);
  }

  // 2D quadratic case   (3noded line element)
  // 3D bilinear case    (4noded quadrilateral element)
  // 3D quadratic case   (6noded triangular element)
  // 3D serendipity case (8noded quadrilateral element)
  // 3D biquadratic case (9noded quadrilateral element)
  else if (dt == line3 || dt == quad4 || dt == tri6 || dt == quad8 || dt == quad9 || dt == nurbs2 ||
           dt == nurbs3 || dt == nurbs4 || dt == nurbs8 || dt == nurbs9)
  {
    // Gauss quadrature with correct NumGP and Dim
    MORTAR::ElementIntegrator integrator(dt);
    double detg = 0.0;

    // loop over all Gauss points, build Jacobian and compute area
    for (int j = 0; j < integrator.nGP(); ++j)
    {
      double gpc[2] = {integrator.Coordinate(j, 0), integrator.Coordinate(j, 1)};
      detg = Jacobian(gpc);
      area += integrator.Weight(j) * detg;

      GEN::pairedvector<int, double> derivjac(NumNode() * Dim());
      DerivJacobian(gpc, derivjac);
      for (GEN::pairedvector<int, double>::const_iterator p = derivjac.begin(); p != derivjac.end();
           ++p)
        area_deriv[p->first] += integrator.Weight(j) * p->second;
    }
  }

  // other cases not implemented yet
  else
    dserror("ERROR: Area computation not implemented for this type of MortarElement");

  return area;
}


/*----------------------------------------------------------------------*
 |  Get global coords for given local coords                  popp 01/08|
 *----------------------------------------------------------------------*/
bool MORTAR::MortarElement::LocalToGlobal(const double* xi, double* globcoord, int inttype)
{
  // check input
  if (!xi) dserror("ERROR: LocalToGlobal called with xi=NULL");
  if (!globcoord) dserror("ERROR: LocalToGlobal called with globcoord=NULL");

  // collect fundamental data
  const int nnodes = NumNode();

  DRT::Node** mynodes = Nodes();
  if (!mynodes) dserror("ERROR: LocalToGlobal: Null pointer!");
  LINALG::SerialDenseMatrix coord(3, nnodes);
  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes, 2, true);

  // Evaluate shape, get nodal coords  and interpolate global coords
  EvaluateShape(xi, val, deriv, nnodes, false);
  GetNodalCoords(coord);

  // init globcoords
  for (int i = 0; i < 3; ++i) globcoord[i] = 0.0;

  for (int i = 0; i < nnodes; ++i)
  {
    if (inttype == 0)
    {
      // use shape function values for interpolation
      globcoord[0] += val[i] * coord(0, i);
      globcoord[1] += val[i] * coord(1, i);
      globcoord[2] += val[i] * coord(2, i);
    }
    else if (inttype == 1)
    {
      // use shape function derivatives xi for interpolation
      globcoord[0] += deriv(i, 0) * coord(0, i);
      globcoord[1] += deriv(i, 0) * coord(1, i);
      globcoord[2] += deriv(i, 0) * coord(2, i);
    }
    else if (inttype == 2)
    {
      // use shape function derivatives eta for interpolation
      globcoord[0] += deriv(i, 1) * coord(0, i);
      globcoord[1] += deriv(i, 1) * coord(1, i);
      globcoord[2] += deriv(i, 1) * coord(2, i);
    }
    else
      dserror("ERROR: Invalid interpolation type requested, only 0,1,2!");
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Compute minimal edge size of MortarElement                popp 11/08|
 *----------------------------------------------------------------------*/
double MORTAR::MortarElement::MinEdgeSize()
{
  double minedgesize = 1.0e12;
  DRT::Element::DiscretizationType dt = Shape();

  // get coordinates of element nodes
  LINALG::SerialDenseMatrix coord(3, NumPoint());
  GetNodalCoords(coord);

  // 2D case (2noded and 3noded line elements)
  if (dt == line2 || dt == line3)
  {
    // there is only one edge
    // (we approximate the quadratic case as linear)
    double diff[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < 3; ++k) diff[k] = coord(k, 1) - coord(k, 0);
    minedgesize = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
  }

  // 3D tri case (3noded and 6noded triangular elements)
  else if (dt == tri3 || dt == tri6)
  {
    // there are three edges
    // (we approximate the quadratic case as linear)
    for (int edge = 0; edge < 3; ++edge)
    {
      double diff[3] = {0.0, 0.0, 0.0};
      for (int k = 0; k < 3; ++k)
      {
        if (edge == 2)
          diff[k] = coord(k, 0) - coord(k, edge);
        else
          diff[k] = coord(k, edge + 1) - coord(k, edge);
      }
      double dist = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
      if (dist < minedgesize) minedgesize = dist;
    }
  }

  // 3D quad case (4noded, 8noded and 9noded quadrilateral elements)
  else if (dt == quad4 || dt == quad8 || dt == quad9)
  {
    // there are four edges
    // (we approximate the quadratic case as linear)
    for (int edge = 0; edge < 4; ++edge)
    {
      double diff[3] = {0.0, 0.0, 0.0};
      for (int k = 0; k < 3; ++k)
      {
        if (edge == 3)
          diff[k] = coord(k, 0) - coord(k, edge);
        else
          diff[k] = coord(k, edge + 1) - coord(k, edge);
      }
      double dist = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
      if (dist < minedgesize) minedgesize = dist;
    }
  }

  //==================================================
  //                     NURBS
  //==================================================
  else if (dt == nurbs3)
  {
    double sxi0[2] = {-1.0, 0.0};
    double sxi1[2] = {1.0, 0.0};
    int nrow = NumNode();
    LINALG::SerialDenseVector sval0(nrow);
    LINALG::SerialDenseVector sval1(nrow);
    LINALG::SerialDenseMatrix sderiv(nrow, 1);
    EvaluateShape(sxi0, sval0, sderiv, nrow);
    EvaluateShape(sxi1, sval1, sderiv, nrow);

    double gpx0[3] = {0.0, 0.0, 0.0};
    double gpx1[3] = {0.0, 0.0, 0.0};

    for (int j = 0; j < nrow; ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        gpx0[i] += sval0(j) * coord(i, j);
        gpx1[i] += sval1(j) * coord(i, j);
      }
    }

    double diff[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < 3; ++k) diff[k] = gpx1[k] - gpx0[k];
    minedgesize = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
  }
  else if (dt == nurbs9)
  {
    int nrow = NumNode();

    // get real point data
    LINALG::SerialDenseMatrix coordnurbs(3, nrow, true);

    // parameter space coordinates
    double sxi0[2] = {-1.0, -1.0};
    double sxi1[2] = {1.0, -1.0};
    double sxi2[2] = {1.0, 1.0};
    double sxi3[2] = {-1.0, 1.0};

    // evaluate shape functions at these coordinates
    LINALG::SerialDenseVector sval0(nrow);
    LINALG::SerialDenseVector sval1(nrow);
    LINALG::SerialDenseVector sval2(nrow);
    LINALG::SerialDenseVector sval3(nrow);
    LINALG::SerialDenseMatrix sderiv(nrow, 2);
    EvaluateShape(sxi0, sval0, sderiv, nrow);
    EvaluateShape(sxi1, sval1, sderiv, nrow);
    EvaluateShape(sxi2, sval2, sderiv, nrow);
    EvaluateShape(sxi3, sval3, sderiv, nrow);

    double gpx0[3] = {0.0, 0.0, 0.0};
    double gpx1[3] = {0.0, 0.0, 0.0};

    for (int j = 0; j < nrow; ++j)
    {
      for (int i = 0; i < 3; ++i)
      {
        gpx0[i] += sval0(j) * coord(i, j);
        gpx1[i] += sval1(j) * coord(i, j);
      }
    }

    // there are four edges
    // (we approximate the quadratic case as linear)
    for (int edge = 0; edge < 4; ++edge)
    {
      double diff[3] = {0.0, 0.0, 0.0};
      for (int k = 0; k < 3; ++k)
      {
        if (edge == 0)
          diff[k] = coord(k, 0) - coord(k, 2);
        else if (edge == 1)
          diff[k] = coord(k, 2) - coord(k, 8);
        else if (edge == 2)
          diff[k] = coord(k, 8) - coord(k, 6);
        else if (edge == 3)
          diff[k] = coord(k, 6) - coord(k, 0);
        else
          dserror("Wrong edge size!");
      }
      double dist = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
      if (dist < minedgesize) minedgesize = dist;
    }
  }
  // invalid case
  else
    dserror("ERROR: MinEdgeSize not implemented for this type of MortarElement");

  if (minedgesize == 1.0e12) dserror("ERROR: MinEdgeSize went wrong...!");
  return minedgesize;
}

/*----------------------------------------------------------------------*
 |  Compute maximal edge size of MortarElement                popp 11/08|
 *----------------------------------------------------------------------*/
double MORTAR::MortarElement::MaxEdgeSize()
{
  double maxedgesize = 0.0;
  DRT::Element::DiscretizationType dt = Shape();

  // get coordinates of element nodes
  LINALG::SerialDenseMatrix coord(3, NumPoint());
  GetNodalCoords(coord);

  // 2D case (2noded and 3noded line elements)
  if (dt == line2 || dt == line3)
  {
    // there is only one edge
    // (we approximate the quadratic case as linear)
    double diff[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < 3; ++k) diff[k] = coord(k, 1) - coord(k, 0);
    maxedgesize = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
  }

  // 3D tri case (3noded and 6noded triangular elements)
  else if (dt == tri3 || dt == tri6)
  {
    // there are three edges
    // (we approximate the quadratic case as linear)
    for (int edge = 0; edge < 3; ++edge)
    {
      double diff[3] = {0.0, 0.0, 0.0};
      for (int k = 0; k < 3; ++k)
      {
        if (edge == 2)
          diff[k] = coord(k, 0) - coord(k, edge);
        else
          diff[k] = coord(k, edge + 1) - coord(k, edge);
      }
      double dist = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
      if (dist > maxedgesize) maxedgesize = dist;
    }
  }

  // 3D quad case (4noded, 8noded and 9noded quadrilateral elements)
  else if (dt == quad4 || dt == quad8 || dt == quad9)
  {
    // there are four edges
    // (we approximate the quadratic case as linear)
    for (int edge = 0; edge < 4; ++edge)
    {
      double diff[3] = {0.0, 0.0, 0.0};
      for (int k = 0; k < 3; ++k)
      {
        if (edge == 3)
          diff[k] = coord(k, 0) - coord(k, edge);
        else
          diff[k] = coord(k, edge + 1) - coord(k, edge);
      }
      double dist = sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]);
      if (dist > maxedgesize) maxedgesize = dist;
    }
  }

  // invalid case
  else
    dserror("ERROR: MaxEdgeSize not implemented for this type of MortarElement");

  if (maxedgesize < 1e-12) dserror("ERROR: MaxEdgeSize went wrong...!");
  return maxedgesize;
}

/*----------------------------------------------------------------------*
 |  Initialize data container                                 popp 12/10|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElement::InitializeDataContainer()
{
  // only initialize if not yet done
  if (modata_ == Teuchos::null) modata_ = Teuchos::rcp(new MORTAR::MortarEleDataContainer());

  if (ParentElement() != NULL)
  {
    int numdof =
        ParentElement()->NumNode() * ParentElement()->NumDofPerNode(*ParentElement()->Nodes()[0]);
    MoData().ParentDisp() = std::vector<double>(numdof);
    for (int i = 0; i < numdof; ++i) MoData().ParentDisp()[i] = 0.0;
  }
}

/*----------------------------------------------------------------------*
 |  Reset data container                                      popp 12/10|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElement::ResetDataContainer()
{
  // reset to Teuchos::null
  modata_ = Teuchos::null;

  return;
}

/*----------------------------------------------------------------------*
 |  Add one MortarElement to potential contact partners       popp 01/08|
 *----------------------------------------------------------------------*/
bool MORTAR::MortarElement::AddSearchElements(const int& gid)
{
  // check calling element type
  if (!IsSlave()) dserror("ERROR: AddSearchElements called for infeasible MortarElement!");

  // add new gid to vector of search candidates
  MoData().SearchElements().push_back(gid);

  return true;
}

/*----------------------------------------------------------------------*
 |  reset found search elements                              farah 10/13|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElement::DeleteSearchElements()
{
  // check calling element type
  if (!IsSlave()) dserror("ERROR: DeleteSearchElements called for infeasible MortarElement!");

  // add new gid to vector of search candidates
  MoData().SearchElements().clear();

  return;
}

/*----------------------------------------------------------------------*
 |  Derivatives of nodal spatial coords                      seitz 03/15|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElement::NodeLinearization(
    std::vector<std::vector<GEN::pairedvector<int, double>>>& nodelin)
{
  // resize the linearizations
  nodelin.resize(NumNode(), std::vector<GEN::pairedvector<int, double>>(3, 1));

  // loop over all intEle nodes
  for (int in = 0; in < NumNode(); ++in)
  {
    MORTAR::MortarNode* mrtrnode = dynamic_cast<MORTAR::MortarNode*>(Nodes()[in]);
    for (int dim = 0; dim < Dim(); ++dim) nodelin[in][dim][mrtrnode->Dofs()[dim]] += 1.;
  }
}

/*----------------------------------------------------------------------*
 |                                                           seitz 11/16|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElement::EstimateNitscheTraceMaxEigenvalueCombined()
{
  if (Dim() != 3)
    dserror(
        "Contact using Nitsche's method is only supported for 3D problems."
        "We do not intend to support 2D problems.");

  Teuchos::RCP<DRT::Element> surf_ele = ParentElement()->Surfaces()[FaceParentNumber()];
  DRT::ELEMENTS::StructuralSurface* surf =
      dynamic_cast<DRT::ELEMENTS::StructuralSurface*>(surf_ele.get());

  traceHE_ = 1. / surf->EstimateNitscheTraceMaxEigenvalueCombined(MoData().ParentDisp());

  if (ParentElement()->NumMaterial() > 1)
    if (ParentElement()->Material(1)->MaterialType() == INPAR::MAT::m_th_fourier_iso)
      traceHCond_ = 1. / surf->EstimateNitscheTraceMaxEigenvalueTSI(MoData().ParentDisp());
}


/*----------------------------------------------------------------------*
 |                                                           seitz 10/16|
 *----------------------------------------------------------------------*/
MORTAR::MortarElementNitscheContainer& MORTAR::MortarElement::GetNitscheContainer()
{
  if (!ParentElement()) dserror("parent element pointer not set");
  if (nitsche_container_ == Teuchos::null) switch (ParentElement()->Shape())
    {
      case DRT::Element::hex8:
        nitsche_container_ =
            Teuchos::rcp(new MORTAR::MortarElementNitscheData<DRT::Element::hex8>());
        break;
      case DRT::Element::tet4:
        nitsche_container_ =
            Teuchos::rcp(new MORTAR::MortarElementNitscheData<DRT::Element::tet4>());
        break;
      case DRT::Element::hex27:
        nitsche_container_ =
            Teuchos::rcp(new MORTAR::MortarElementNitscheData<DRT::Element::hex27>());
        break;
      case DRT::Element::nurbs27:
        nitsche_container_ =
            Teuchos::rcp(new MORTAR::MortarElementNitscheData<DRT::Element::nurbs27>());
        break;
      default:
        dserror("Nitsche data container not ready. Just add it here...");
    }
  return *nitsche_container_;
}
