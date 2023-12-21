/*----------------------------------------------------------------------------*/
/*! \file
\brief three dimensional total Lagrange truss element

\level 3


*/
/*---------------------------------------------------------------------------*/

#include "baci_truss3.H"

#include "baci_io_linedefinition.H"
#include "baci_lib_discret.H"
#include "baci_lib_node.H"
#include "baci_so3_nullspace.H"
#include "baci_structure_new_elements_paramsinterface.H"

BACI_NAMESPACE_OPEN

DRT::ELEMENTS::Truss3Type DRT::ELEMENTS::Truss3Type::instance_;

DRT::ELEMENTS::Truss3Type& DRT::ELEMENTS::Truss3Type::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::Truss3Type::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::Truss3(-1, -1);
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

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::Truss3Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void DRT::ELEMENTS::Truss3Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions["TRUSS3"];

  defs["LINE2"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("LINE2", 2)
                      .AddNamedInt("MAT")
                      .AddNamedDouble("CROSS")
                      .AddNamedString("KINEM")
                      .Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 08/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss3::Truss3(int id, int owner)
    : DRT::Element(id, owner),
      crosssec_(0.0),
      eint_(0.0),
      lrefe_(0.0),
      gaussrule_(CORE::FE::GaussRule1D::line_2point),
      data_(),
      diff_disp_ref_(CORE::LINALG::Matrix<1, 3>(true)),
      interface_ptr_(Teuchos::null),
      isinit_(false),
      jacobimass_(),
      jacobinode_(),
      kintype_(KinematicType::tr3_totlag),
      material_(0),
      X_()
{
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 08/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss3::Truss3(const DRT::ELEMENTS::Truss3& old)
    : DRT::Element(old),
      crosssec_(old.crosssec_),
      eint_(old.eint_),
      lrefe_(old.lrefe_),
      gaussrule_(old.gaussrule_),
      data_(old.data_),
      diff_disp_ref_(old.diff_disp_ref_),
      interface_ptr_(old.interface_ptr_),
      isinit_(old.isinit_),
      jacobimass_(old.jacobimass_),
      jacobinode_(old.jacobinode_),
      kintype_(old.kintype_),
      material_(old.material_),
      X_(old.X_)

{
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Truss3 and return pointer to it (public) |
 |                                                            cyron 08/08|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Truss3::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::Truss3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |(public)                                                   cyron 08/08|
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::Truss3::Shape() const { return CORE::FE::CellType::line2; }


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);
  AddtoPack(data, isinit_);
  AddtoPack<6, 1>(data, X_);
  AddtoPack<1, 3>(data, diff_disp_ref_);
  AddtoPack(data, material_);
  AddtoPack(data, lrefe_);
  AddtoPack(data, jacobimass_);
  AddtoPack(data, jacobinode_);
  AddtoPack(data, crosssec_);
  AddtoPack(data, gaussrule_);
  AddtoPack(data, kintype_);
  AddtoPack(data, data_);
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);
  isinit_ = ExtractInt(position, data);
  ExtractfromPack<6, 1>(position, data, X_);
  ExtractfromPack<1, 3>(position, data, diff_disp_ref_);
  ExtractfromPack(position, data, material_);
  ExtractfromPack(position, data, lrefe_);
  ExtractfromPack(position, data, jacobimass_);
  ExtractfromPack(position, data, jacobinode_);
  ExtractfromPack(position, data, crosssec_);
  ExtractfromPack(position, data, gaussrule_);
  // kinematic type
  kintype_ = static_cast<KinematicType>(ExtractInt(position, data));
  std::vector<char> tmp(0);
  ExtractfromPack(position, data, tmp);
  data_.Unpack(tmp);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                              cyron 08/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Truss3::Lines()
{
  return {Teuchos::rcpFromRef(*this)};
}

/*----------------------------------------------------------------------*
 |determine Gauss rule from required type of integration                |
 |                                                   (public)cyron 09/09|
 *----------------------------------------------------------------------*/
CORE::FE::GaussRule1D DRT::ELEMENTS::Truss3::MyGaussRule(int nnode, IntegrationType integrationtype)
{
  CORE::FE::GaussRule1D gaussrule = CORE::FE::GaussRule1D::undefined;

  switch (nnode)
  {
    case 2:
    {
      switch (integrationtype)
      {
        case gaussexactintegration:
        {
          gaussrule = CORE::FE::GaussRule1D::line_2point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule = CORE::FE::GaussRule1D::line_1point;
          break;
        }
        case lobattointegration:
        {
          gaussrule = CORE::FE::GaussRule1D::line_lobatto2point;
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
          gaussrule = CORE::FE::GaussRule1D::line_3point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule = CORE::FE::GaussRule1D::line_2point;
          break;
        }
        case lobattointegration:
        {
          gaussrule = CORE::FE::GaussRule1D::line_lobatto3point;
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
          gaussrule = CORE::FE::GaussRule1D::line_4point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule = CORE::FE::GaussRule1D::line_3point;
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
          gaussrule = CORE::FE::GaussRule1D::line_5point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule = CORE::FE::GaussRule1D::line_4point;
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
void DRT::ELEMENTS::Truss3::SetUpReferenceGeometry(const std::vector<double>& xrefe)
{
  if (!isinit_)
  {
    // setting reference coordinates
    for (int i = 0; i < 6; ++i) X_(i) = xrefe[i];

    // length in reference configuration
    lrefe_ = std::sqrt((X_(3) - X_(0)) * (X_(3) - X_(0)) + (X_(4) - X_(1)) * (X_(4) - X_(1)) +
                       (X_(5) - X_(2)) * (X_(5) - X_(2)));

    // set jacobi determinants for integration of mass matrix and at nodes
    jacobimass_.resize(2);
    jacobimass_[0] = lrefe_ / 2.0;
    jacobimass_[1] = lrefe_ / 2.0;
    jacobinode_.resize(2);
    jacobinode_[0] = lrefe_ / 2.0;
    jacobinode_[1] = lrefe_ / 2.0;

    isinit_ = true;
  }
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
          const std::vector<int>* flag = nullptr;
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
    std::vector<int> dofs = dis.Dof(dofset, this);
    if (dofs.size()) lmstride.push_back(dofs.size());
    for (int& dof : dofs)
    {
      lmowner.push_back(owner);
      lm.push_back(dof);
    }

    if (doDirichlet)
    {
      const std::vector<int>* flag = nullptr;
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
      for (unsigned j = 0; j < dofs.size(); ++j)
      {
        if (flag && (*flag)[j])
          lmdirich.push_back(1);
        else
          lmdirich.push_back(0);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Truss3Type::Initialize(DRT::Discretization& dis)
{
  // reference node positions
  std::vector<double> xrefe;

  // reference nodal tangent positions
  CORE::LINALG::Matrix<3, 1> trefNodeAux(true);
  // resize vectors for the number of coordinates we need to store
  xrefe.resize(3 * 2);

  // setting beam reference director correctly
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    // in case that current element is not a truss3 element there is nothing to do and we go back
    // to the head of the loop
    if (dis.lColElement(i)->ElementType() != *this) continue;

    // if we get so far current element is a truss3 element and  we get a pointer at it
    auto* currele = dynamic_cast<DRT::ELEMENTS::Truss3*>(dis.lColElement(i));
    if (!currele) dserror("cast to Truss3* failed");

    // getting element's nodal coordinates and treating them as reference configuration
    if (currele->Nodes()[0] == nullptr || currele->Nodes()[1] == nullptr)
      dserror("Cannot get nodes in order to compute reference configuration'");
    else
    {
      for (int k = 0; k < 2; k++)  // element has two nodes
        for (int l = 0; l < 3; l++) xrefe[k * 3 + l] = currele->Nodes()[k]->X()[l];
    }

    currele->SetUpReferenceGeometry(xrefe);
  }

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::SetParamsInterfacePtr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
  {
    interface_ptr_ = Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>(
        p.get<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface"));
  }
  else
    interface_ptr_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> DRT::ELEMENTS::Truss3::ParamsInterfacePtr()
{
  return interface_ptr_;
}

BACI_NAMESPACE_CLOSE
