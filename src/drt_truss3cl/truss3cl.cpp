/*!----------------------------------------------------------------------
\file truss3cl.cpp
 \brief three dimensional interpolated total Lagrange hybrid beam-truss element
 (can be connected to beam3 elements)

\level 3

<pre>
\maintainer Dhrubajyoti Mukherjee
            mukherjee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15270
</pre>

*----------------------------------------------------------------------*/

#include "truss3cl.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

DRT::ELEMENTS::Truss3CLType DRT::ELEMENTS::Truss3CLType::instance_;

DRT::ELEMENTS::Truss3CLType& DRT::ELEMENTS::Truss3CLType::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::Truss3CLType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Truss3CL* object = new DRT::ELEMENTS::Truss3CL(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Truss3CLType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "Truss3CL")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Truss3CL(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Truss3CLType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Truss3CL(id, owner));
  return ele;
}


void DRT::ELEMENTS::Truss3CLType::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

void DRT::ELEMENTS::Truss3CLType::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  DRT::UTILS::ComputeStructure3DNullSpace(dis, ns, x0, numdf, dimns);
}

void DRT::ELEMENTS::Truss3CLType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["Truss3CL"];

  defs["LINE2"]
      .AddIntVector("LINE2", 2)
      .AddNamedInt("MAT")
      .AddNamedDouble("CROSS")
      .AddNamedString("KINEM");

  defs["LIN2"].AddIntVector("LIN2", 2).AddNamedInt("MAT").AddNamedDouble("CROSS").AddNamedString(
      "KINEM");
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                        mukherjee 01/14|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss3CL::Truss3CL(int id, int owner)
    : DRT::Element(id, owner),
      data_(),
      isinit_(false),
      material_(0),
      lrefe_(0),
      crosssec_(0),
      kintype_(tr3_totlag),
      // note: for corotational approach integration for Neumann conditions only
      // hence enough to integrate 3rd order polynomials exactly
      gaussrule_(DRT::UTILS::intrule_line_2point),
      mybindingposition_(0)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                   mukherjee 01/14|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss3CL::Truss3CL(const DRT::ELEMENTS::Truss3CL& old)
    : DRT::Element(old),
      data_(old.data_),
      isinit_(old.isinit_),
      X_(old.X_),
      material_(old.material_),
      lrefe_(old.lrefe_),
      xrefe_(old.xrefe_),
      jacobimass_(old.jacobimass_),
      jacobinode_(old.jacobinode_),
      crosssec_(old.crosssec_),
      kintype_(old.kintype_),
      gaussrule_(old.gaussrule_),
      mybindingposition_(old.mybindingposition_)
{
  return;
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Truss3CL and return pointer to it (public) |
 |                                                        mukherjee 01/14|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Truss3CL::Clone() const
{
  DRT::ELEMENTS::Truss3CL* newelement = new DRT::ELEMENTS::Truss3CL(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                        mukherjee 01/14|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Truss3CL::~Truss3CL() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                          mukherjee 01/14|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3CL::Print(std::ostream& os) const
{
  os << "Truss3CL ";
  Element::Print(os);
  os << " gaussrule_: " << gaussrule_ << " ";
  return;
}


/*----------------------------------------------------------------------*
 |(public)                                               mukherjee 01/14|
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Truss3CL::Shape() const { return line2; }


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                       mukherjee 01/14|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3CL::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);
  AddtoPack(data, isinit_);
  AddtoPack(data, X_);
  AddtoPack(data, material_);
  AddtoPack(data, lrefe_);
  AddtoPack(data, jacobimass_);
  AddtoPack(data, jacobinode_);
  AddtoPack(data, crosssec_);
  AddtoPack(data, gaussrule_);  // implicit conversion from enum to integer
  AddtoPack(data, kintype_);
  AddtoPack(data, data_);
  AddtoPack(data, xrefe_);
  AddtoPack(data, mybindingposition_);
  AddtoPack(data, xrefe_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                       mukherjee 01/14|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3CL::Unpack(const std::vector<char>& data)
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
  ExtractfromPack(position, data, X_);
  ExtractfromPack(position, data, material_);
  ExtractfromPack(position, data, lrefe_);
  ExtractfromPack(position, data, jacobimass_);
  ExtractfromPack(position, data, jacobinode_);
  ExtractfromPack(position, data, crosssec_);
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
  ExtractfromPack(position, data, xrefe_);
  ExtractfromPack(position, data, mybindingposition_);
  ExtractfromPack(position, data, xrefe_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                          mukherjee 01/14|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Truss3CL::Lines()
{
  std::vector<Teuchos::RCP<Element>> lines(1);
  lines[0] = Teuchos::rcp(this, false);
  return lines;
}

/*----------------------------------------------------------------------*
 |determine Gauss rule from required type of integration                |
 |                                               (public)mukherjee 01/14|
 *----------------------------------------------------------------------*/
DRT::UTILS::GaussRule1D DRT::ELEMENTS::Truss3CL::MyGaussRule(
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

void DRT::ELEMENTS::Truss3CL::SetUpReferenceGeometry(
    const std::vector<double>& xrefe, const bool secondinit)
{
  /*this method initializes geometric variables of the element; the initilization can usually be
   *applied to elements only once; therefore after the first initilization the flag isinit is set to
   *true and from then on this method does not take any action when called again unless it is called
   *on purpose with the additional parameter secondinit. If this parameter is passed into the method
   *and is true the element is initialized another time with respective xrefe; note: the isinit_
   *flag is important for avoiding reinitialization upon restart. However, it should be possible to
   *conduct a second initilization in principle (e.g. for periodic boundary conditions*/
  if (!isinit_ || secondinit)
  {
    isinit_ = true;

    // setting reference coordinates
    for (int i = 0; i < 6; i++) X_(i) = xrefe[i];

    // length in reference configuration
    lrefe_ = std::sqrt(pow(X_(3) - X_(0), 2) + pow(X_(4) - X_(1), 2) + pow(X_(5) - X_(2), 2));

    // set jacobi determinants for integration of mass matrix and at nodes
    jacobimass_.resize(2);
    jacobimass_[0] = lrefe_ / 2.0;
    jacobimass_[1] = lrefe_ / 2.0;
    jacobinode_.resize(2);
    jacobinode_[0] = lrefe_ / 2.0;
    jacobinode_[1] = lrefe_ / 2.0;
  }

  return;
}


/*----------------------------------------------------------------------*
 | Get degrees of freedom used by this element               fang 05/16 |
 *----------------------------------------------------------------------*/
// TODO: remove once truss3cl element is fixed and no longer expects more dofs (6) than it can
// inherently handle (3)...
void DRT::ELEMENTS::Truss3CL::LocationVector(
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


int DRT::ELEMENTS::Truss3CLType::Initialize(DRT::Discretization& dis)
{
  const int nnode = 4;
  const int fnnode = 2;
  // setting truss reference director correctly
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    // in case that current element is not a Truss3CL element there is nothing to do and we go back
    // to the head of the loop
    if (dis.lColElement(i)->ElementType() != *this) continue;

    // if we get so far current element is a Truss3CL element and  we get a pointer at it
    DRT::ELEMENTS::Truss3CL* currele = dynamic_cast<DRT::ELEMENTS::Truss3CL*>(dis.lColElement(i));
    if (!currele) dserror("cast to Truss3CL* failed");

    // reference node position of real nodes
    LINALG::Matrix<12, 1> xrefe;
    xrefe.Clear();
    currele->xrefe_.resize(3 * fnnode);
    currele->xrefe_.clear();

    // getting element's nodal coordinates and treating them as reference configuration
    if (currele->Nodes()[0] == NULL || currele->Nodes()[1] == NULL)
      dserror("Cannot get nodes in order to compute reference configuration'");
    else
    {
      for (int node = 0; node < nnode; node++)  // element has four nodes
        for (int dof = 0; dof < 3; dof++) xrefe(node * 3 + dof) = currele->Nodes()[node]->X()[dof];

      // since filaments are straight at initial configuration, only lagrange polynomials are enough

      std::vector<LINALG::Matrix<1, 2>> Ibp(2);
      for (int filament = 0; filament < 2; filament++)
      {
        DRT::UTILS::shape_function_1D(
            Ibp[filament], currele->mybindingposition_[filament], currele->Shape());
      }

      for (int filament = 0; filament < 2; filament++)
        for (int node = 0; node < fnnode; node++)
          for (int dof = 0; i < 3; i++)
          {
            currele->xrefe_[i + 3 * filament] +=
                Ibp[filament](node) * xrefe(dof + 3 * node + 6 * filament);
          }
    }
    currele->SetUpReferenceGeometry(currele->xrefe_);
  }  // for (int i=0; i<dis_.NumMyColElements(); ++i)
  return 0;
}
