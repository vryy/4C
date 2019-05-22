/*----------------------------------------------------------------------------*/
/*!
\brief three dimensional total Lagrange truss element

\level 3

\maintainer Maximilian Grill

*/
/*---------------------------------------------------------------------------*/

#include "truss3.H"
#include "../drt_beam3/beam3eb.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils_nullspace.H"
//#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_structure_new/str_elements_paramsinterface.H"

DRT::ELEMENTS::Truss3Type DRT::ELEMENTS::Truss3Type::instance_;

DRT::ELEMENTS::Truss3Type& DRT::ELEMENTS::Truss3Type::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::Truss3Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Truss3* object = new DRT::ELEMENTS::Truss3(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Truss3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "TRUSS3")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Truss3(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Truss3Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Truss3(id, owner));
  return ele;
}


void DRT::ELEMENTS::Truss3Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

void DRT::ELEMENTS::Truss3Type::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  DRT::UTILS::ComputeStructure3DNullSpace(dis, ns, x0, numdf, dimns);
}

void DRT::ELEMENTS::Truss3Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["TRUSS3"];

  defs["LINE2"]
      .AddIntVector("LINE2", 2)
      .AddNamedInt("MAT")
      .AddNamedDouble("CROSS")
      .AddNamedString("KINEM");

  defs["LIN2"].AddIntVector("LIN2", 2).AddNamedInt("MAT").AddNamedDouble("CROSS").AddNamedString(
      "KINEM");
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 08/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss3::Truss3(int id, int owner)
    : DRT::Element(id, owner),
      data_(),
      isinit_(false),
      diff_disp_ref_(LINALG::Matrix<1, 3>(true)),
      deltatheta_(LINALG::Matrix<1, 3>(true)),
      material_(0),
      lrefe_(0),
      lcurr_(0),
      crosssec_(0),
      NormMoment(0),
      NormForce(0),
      RatioNormForceMoment(0),
      Theta0_(LINALG::Matrix<3, 1>(true)),
      Theta_(LINALG::Matrix<3, 1>(true)),
      kintype_(tr3_totlag),
      // note: for corotational approach integration for Neumann conditions only
      // hence enough to integrate 3rd order polynomials exactly
      gaussrule_(DRT::UTILS::intrule_line_2point)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 08/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss3::Truss3(const DRT::ELEMENTS::Truss3& old)
    : DRT::Element(old),
      data_(old.data_),
      isinit_(old.isinit_),
      X_(old.X_),
      trefNode_(old.trefNode_),
      ThetaRef_(old.ThetaRef_),
      diff_disp_ref_(old.diff_disp_ref_),
      deltatheta_(old.deltatheta_),
      material_(old.material_),
      lrefe_(old.lrefe_),
      lcurr_(old.lcurr_),
      jacobimass_(old.jacobimass_),
      jacobinode_(old.jacobinode_),
      crosssec_(old.crosssec_),
      NormMoment(old.NormMoment),
      NormForce(old.NormForce),
      RatioNormForceMoment(old.RatioNormForceMoment),
      Theta0_(old.Theta0_),
      Theta_(old.Theta_),
      kintype_(old.kintype_),
      gaussrule_(old.gaussrule_)
{
  return;
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Truss3 and return pointer to it (public) |
 |                                                            cyron 08/08|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Truss3::Clone() const
{
  DRT::ELEMENTS::Truss3* newelement = new DRT::ELEMENTS::Truss3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            cyron 08/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss3::~Truss3() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::Print(std::ostream& os) const
{
  os << "Truss3 ";
  Element::Print(os);
  os << " gaussrule_: " << gaussrule_ << " ";
  return;
}



/*----------------------------------------------------------------------*
 | Print the change in angle of this element            mukherjee 10/14 |
 *----------------------------------------------------------------------*/
LINALG::Matrix<1, 3> DRT::ELEMENTS::Truss3::DeltaTheta() const
{
  // for now constant, since we only implemented 4-noded interpolated element with linear shape
  // functions
  return deltatheta_;
}

/*----------------------------------------------------------------------*
 |(public)                                                   cyron 08/08|
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Truss3::Shape() const { return line2; }


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);
  AddtoPack(data, isinit_);
  AddtoPack<6, 1>(data, X_);
  AddtoPack(data, trefNode_);
  AddtoPack(data, ThetaRef_);
  AddtoPack<1, 3>(data, diff_disp_ref_);
  AddtoPack<1, 3>(data, deltatheta_);
  AddtoPack(data, material_);
  AddtoPack(data, lrefe_);
  AddtoPack(data, lcurr_);
  AddtoPack(data, jacobimass_);
  AddtoPack(data, jacobinode_);
  AddtoPack(data, crosssec_);
  AddtoPack(data, NormMoment);
  AddtoPack(data, NormForce);
  AddtoPack(data, RatioNormForceMoment);
  AddtoPack<3, 1>(data, Theta0_);
  AddtoPack<3, 1>(data, Theta_);
  AddtoPack(data, gaussrule_);  // implicit conversion from enum to integer
  AddtoPack(data, kintype_);
  AddtoPack(data, data_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);
  isinit_ = ExtractInt(position, data);
  ExtractfromPack<6, 1>(position, data, X_);
  ExtractfromPack(position, data, trefNode_);
  ExtractfromPack(position, data, ThetaRef_);
  ExtractfromPack<1, 3>(position, data, diff_disp_ref_);
  ExtractfromPack<1, 3>(position, data, deltatheta_);
  ExtractfromPack(position, data, material_);
  ExtractfromPack(position, data, lrefe_);
  ExtractfromPack(position, data, lcurr_);
  ExtractfromPack(position, data, jacobimass_);
  ExtractfromPack(position, data, jacobinode_);
  ExtractfromPack(position, data, crosssec_);
  ExtractfromPack(position, data, NormMoment);
  ExtractfromPack(position, data, NormForce);
  ExtractfromPack(position, data, RatioNormForceMoment);
  ExtractfromPack<3, 1>(position, data, Theta0_);
  ExtractfromPack<3, 1>(position, data, Theta_);
  // gaussrule_
  int gausrule_integer;
  ExtractfromPack(position, data, gausrule_integer);
  gaussrule_ =
      DRT::UTILS::GaussRule1D(gausrule_integer);  // explicit conversion from integer to enum
  // kinematic type
  kintype_ = static_cast<KinematicType>(ExtractInt(position, data));
  std::vector<char> tmp(0);
  ExtractfromPack(position, data, tmp);
  data_.Unpack(tmp);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                              cyron 08/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Truss3::Lines()
{
  std::vector<Teuchos::RCP<Element>> lines(1);
  lines[0] = Teuchos::rcp(this, false);
  return lines;
}

/*----------------------------------------------------------------------*
 |determine Gauss rule from required type of integration                |
 |                                                   (public)cyron 09/09|
 *----------------------------------------------------------------------*/
DRT::UTILS::GaussRule1D DRT::ELEMENTS::Truss3::MyGaussRule(
    int nnode, IntegrationType integrationtype)
{
  DRT::UTILS::GaussRule1D gaussrule = DRT::UTILS::intrule1D_undefined;

  switch (nnode)
  {
    case 2:
    {
      switch (integrationtype)
      {
        case gaussexactintegration:
        {
          gaussrule = DRT::UTILS::intrule_line_2point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule = DRT::UTILS::intrule_line_1point;
          break;
        }
        case lobattointegration:
        {
          gaussrule = DRT::UTILS::intrule_line_lobatto2point;
          break;
        }
        default:
          dserror("unknown type of integration");
      }
      break;
    }
    case 3:
    {
      switch (integrationtype)
      {
        case gaussexactintegration:
        {
          gaussrule = DRT::UTILS::intrule_line_3point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule = DRT::UTILS::intrule_line_2point;
          break;
        }
        case lobattointegration:
        {
          gaussrule = DRT::UTILS::intrule_line_lobatto3point;
          break;
        }
        default:
          dserror("unknown type of integration");
      }
      break;
    }
    case 4:
    {
      switch (integrationtype)
      {
        case gaussexactintegration:
        {
          gaussrule = DRT::UTILS::intrule_line_4point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule = DRT::UTILS::intrule_line_3point;
          break;
        }
        default:
          dserror("unknown type of integration");
      }
      break;
    }
    case 5:
    {
      switch (integrationtype)
      {
        case gaussexactintegration:
        {
          gaussrule = DRT::UTILS::intrule_line_5point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule = DRT::UTILS::intrule_line_4point;
          break;
        }
        default:
          dserror("unknown type of integration");
      }
      break;
    }
    default:
      dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
  }

  return gaussrule;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::SetUpReferenceGeometry(
    const std::vector<double>& xrefe, const std::vector<double>& rotrefe)
{
  if (isinit_)
    dserror(
        " You are about to setup an element that has already been set up. That does not make "
        "sense.");

  isinit_ = true;

  // setting reference coordinates
  for (int i = 0; i < 6; ++i) X_(i) = xrefe[i];

  // length in reference configuration
  lrefe_ = std::pow((X_(3) - X_(0)) * (X_(3) - X_(0)) + (X_(4) - X_(1)) * (X_(4) - X_(1)) +
                        (X_(5) - X_(2)) * (X_(5) - X_(2)),
      0.5);

  lcurr_ = lrefe_;

  // set jacobi determinants for integration of mass matrix and at nodes
  jacobimass_.resize(2);
  jacobimass_[0] = lrefe_ / 2.0;
  jacobimass_[1] = lrefe_ / 2.0;
  jacobinode_.resize(2);
  jacobinode_[0] = lrefe_ / 2.0;
  jacobinode_[1] = lrefe_ / 2.0;

  double abs_rotrefe = 0;
  for (int i = 0; i < 6; ++i) abs_rotrefe += rotrefe[i] * rotrefe[i];

  if (abs_rotrefe != 0)
  {
    // assign size to the vector
    ThetaRef_.resize(3);
    trefNode_.resize(3);

    for (int node = 0; node < 2; node++)
    {
      trefNode_[node].Clear();
      for (int dof = 0; dof < 3; dof++) trefNode_[node](dof) = rotrefe[3 * node + dof];
    }
    diff_disp_ref_.Clear();

    // Calculate reference directional vector of the truss element
    for (int j = 0; j < 3; ++j)
    {
      diff_disp_ref_(j) = Nodes()[1]->X()[j] - Nodes()[0]->X()[j];
    }

    for (int location = 0; location < 3;
         location++)  // Location of torsional spring. There are three locations
    {
      double dotprod = 0.0;
      LINALG::Matrix<1, 3> crossprod(true);
      double CosTheta = 0.0;
      double SinTheta = 0.0;

      if (location == 0)
      {
        double norm_v_ref = diff_disp_ref_.Norm2();
        double norm_t1_ref = trefNode_[location].Norm2();
        for (int j = 0; j < 3; ++j) dotprod += trefNode_[location](j) * diff_disp_ref_(j);

        CosTheta = dotprod / (norm_v_ref * norm_t1_ref);

        // Cross Product
        crossprod(0) =
            trefNode_[location](1) * diff_disp_ref_(2) - trefNode_[location](2) * diff_disp_ref_(1);
        crossprod(1) =
            trefNode_[location](2) * diff_disp_ref_(0) - trefNode_[location](0) * diff_disp_ref_(2);
        crossprod(2) =
            trefNode_[location](0) * diff_disp_ref_(1) - trefNode_[location](1) * diff_disp_ref_(0);

        double norm = crossprod.Norm2();
        SinTheta = norm / (norm_v_ref * norm_t1_ref);
      }
      else if (location == 1)
      {
        double norm_v_ref = diff_disp_ref_.Norm2();
        double norm_t2_ref = trefNode_[location].Norm2();
        for (int j = 0; j < 3; ++j)
          dotprod +=
              trefNode_[location](j) * diff_disp_ref_(j);  // From the opposite direction v_2 =-v_1

        CosTheta = dotprod / (norm_v_ref * norm_t2_ref);

        // cross product
        crossprod(0) =
            trefNode_[location](1) * diff_disp_ref_(2) - trefNode_[location](2) * diff_disp_ref_(1);
        crossprod(1) =
            trefNode_[location](2) * diff_disp_ref_(0) - trefNode_[location](0) * diff_disp_ref_(2);
        crossprod(2) =
            trefNode_[location](0) * diff_disp_ref_(1) - trefNode_[location](1) * diff_disp_ref_(0);
        double norm = crossprod.Norm2();
        SinTheta = norm / (norm_v_ref * norm_t2_ref);
      }
      else  // i.e. for calculation of reference angle between t1 & t2
      {
        double norm_t1_ref = trefNode_[location - 2].Norm2();
        double norm_t2_ref = trefNode_[location - 1].Norm2();
        for (int j = 0; j < 3; ++j)
          dotprod += trefNode_[location - 1](j) * trefNode_[location - 2](j);

        CosTheta = dotprod / (norm_t1_ref * norm_t2_ref);

        // cross product
        crossprod(0) = trefNode_[location - 2](1) * trefNode_[location - 1](2) -
                       trefNode_[location - 2](2) * trefNode_[location - 1](1);
        crossprod(1) = trefNode_[location - 2](2) * trefNode_[location - 1](0) -
                       trefNode_[location - 2](0) * trefNode_[location - 1](2);
        crossprod(2) = trefNode_[location - 2](0) * trefNode_[location - 1](1) -
                       trefNode_[location - 2](1) * trefNode_[location - 1](0);

        double norm = crossprod.Norm2();
        SinTheta = norm / (norm_t1_ref * norm_t2_ref);
      }

      double ThetaBoundary1 = M_PI / 4;
      double ThetaBoundary2 = 3 * M_PI / 4;

      ThetaRef_[location] = 0;
      if (SinTheta >= 0)
      {
        if (CosTheta >= cos(ThetaBoundary1))
          ThetaRef_[location] = asin(SinTheta);
        else if (CosTheta <= cos(ThetaBoundary2))
          ThetaRef_[location] = M_PI - asin(SinTheta);
        else
          ThetaRef_[location] = acos(CosTheta);
      }
      else
        dserror("Angle more than 180 degrees!");

      Theta0_(location) = ThetaRef_[location];
    }
  }


  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::ScaleReferenceLength(double scalefac)
{
  // scale length in reference configuration
  X_(3) = X_(0) + (scalefac * (X_(3) - X_(0)));
  X_(4) = X_(1) + (scalefac * (X_(4) - X_(1)));
  X_(5) = X_(2) + (scalefac * (X_(5) - X_(2)));

  lrefe_ *= scalefac;

  // set jacobi determinants for integration of mass matrix and at nodes
  jacobimass_.resize(2);
  jacobimass_[0] = jacobimass_[1] = lrefe_ * 0.5;
  jacobinode_.resize(2);
  jacobinode_[0] = jacobinode_[1] = lrefe_ * 0.5;
}

/*----------------------------------------------------------------------*
 | Get degrees of freedom used by this element               fang 05/16 |
 *----------------------------------------------------------------------*/
// TODO: remove once truss3 element is fixed and no longer expects more dofs (6) than it can
// inherently handle (3)...
void DRT::ELEMENTS::Truss3::LocationVector(
    const Discretization& dis, LocationArray& la, bool doDirichlet) const
{
  const int numnode = NumNode();
  const DRT::Node* const* nodes = Nodes();

  la.Clear();

  // we need to look at all DofSets of our Discretization
  for (int dofset = 0; dofset < la.Size(); ++dofset)
  {
    std::vector<int>& lm = la[dofset].lm_;
    std::vector<int>& lmdirich = la[dofset].lmdirich_;
    std::vector<int>& lmowner = la[dofset].lmowner_;
    std::vector<int>& lmstride = la[dofset].stride_;

    // fill the vector with nodal dofs
    if (nodes)
    {
      for (int i = 0; i < numnode; ++i)
      {
        const DRT::Node* node = nodes[i];

        const int owner = node->Owner();
        std::vector<int> dof;
        dis.Dof(dof, node, dofset, 0, this);
        const unsigned size = dof.size();

        if (size) lmstride.push_back(size);
        for (unsigned j = 0; j < size; ++j)
        {
          lmowner.push_back(owner);
          lm.push_back(dof[j]);
        }

        if (doDirichlet)
        {
          const std::vector<int>* flag = NULL;
          DRT::Condition* dirich = node->GetCondition("Dirichlet");
          if (dirich)
          {
            if (dirich->Type() != DRT::Condition::PointDirichlet &&
                dirich->Type() != DRT::Condition::LineDirichlet &&
                dirich->Type() != DRT::Condition::SurfaceDirichlet &&
                dirich->Type() != DRT::Condition::VolumeDirichlet)
              dserror("condition with name Dirichlet is not of type Dirichlet");
            flag = dirich->Get<std::vector<int>>("onoff");
          }
          for (unsigned j = 0; j < size; ++j)
          {
            if (flag && (*flag)[j])
              lmdirich.push_back(1);
            else
              lmdirich.push_back(0);
          }
        }
      }
    }

    // fill the vector with element dofs
    const int owner = Owner();
    std::vector<int> dof = dis.Dof(dofset, this);
    if (dof.size()) lmstride.push_back(dof.size());
    for (unsigned j = 0; j < dof.size(); ++j)
    {
      lmowner.push_back(owner);
      lm.push_back(dof[j]);
    }

    if (doDirichlet)
    {
      const std::vector<int>* flag = NULL;
      DRT::Condition* dirich = GetCondition("Dirichlet");
      if (dirich)
      {
        if (dirich->Type() != DRT::Condition::PointDirichlet &&
            dirich->Type() != DRT::Condition::LineDirichlet &&
            dirich->Type() != DRT::Condition::SurfaceDirichlet &&
            dirich->Type() != DRT::Condition::VolumeDirichlet)
          dserror("condition with name Dirichlet is not of type Dirichlet");
        flag = dirich->Get<std::vector<int>>("onoff");
      }
      for (unsigned j = 0; j < dof.size(); ++j)
      {
        if (flag && (*flag)[j])
          lmdirich.push_back(1);
        else
          lmdirich.push_back(0);
      }
    }

  }  // for (int dofset=0; dofset<la.Size(); ++dofset)

  return;
}


int DRT::ELEMENTS::Truss3Type::Initialize(DRT::Discretization& dis)
{
  //  dserror("stop. Please modify the way reference tangent is calculated!");
  // reference node positions
  std::vector<double> xrefe;

  // reference nodal tangent positions
  std::vector<double> rotrefe;
  LINALG::Matrix<3, 1> trefNodeAux(true);
  // resize vectors for the number of coordinates we need to store
  xrefe.resize(3 * 2);
  rotrefe.resize(3 * 2);
  for (int i = 0; i < 6; i++) rotrefe[i] = 0;

  // setting beam reference director correctly
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    // in case that current element is not a truss3 element there is nothing to do and we go back
    // to the head of the loop
    if (dis.lColElement(i)->ElementType() != *this) continue;

    // if we get so far current element is a truss3 element and  we get a pointer at it
    DRT::ELEMENTS::Truss3* currele = dynamic_cast<DRT::ELEMENTS::Truss3*>(dis.lColElement(i));
    if (!currele) dserror("cast to Truss3* failed");

    // getting element's nodal coordinates and treating them as reference configuration
    if (currele->Nodes()[0] == NULL || currele->Nodes()[1] == NULL)
      dserror("Cannot get nodes in order to compute reference configuration'");
    else
    {
      for (int k = 0; k < 2; k++)  // element has two nodes
        for (int l = 0; l < 3; l++) xrefe[k * 3 + l] = currele->Nodes()[k]->X()[l];
    }

    // ask the truss element about the first element the first node is connected to
    DRT::Element* Element = currele->Nodes()[0]->Elements()[0];
    // Check via dynamic cast, if it's a beam3eb element
    DRT::ELEMENTS::Beam3eb* BeamElement = dynamic_cast<DRT::ELEMENTS::Beam3eb*>(Element);
    if (BeamElement != NULL)
    {
      for (int k = 0; k < 2; k++)  // element has two nodes
        for (int l = 0; l < 3; l++)
        {
          trefNodeAux = BeamElement->Tref()[k];
          rotrefe[k * 3 + l] = trefNodeAux(l);
        }
    }

    currele->SetUpReferenceGeometry(xrefe, rotrefe);


  }  // for (int i=0; i<dis_.NumMyColElements(); ++i)


  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::SetParamsInterfacePtr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
    interface_ptr_ = Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>(
        p.get<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface"));
  else
    interface_ptr_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> DRT::ELEMENTS::Truss3::ParamsInterfacePtr()
{
  return interface_ptr_;
}
