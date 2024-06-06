/*----------------------------------------------------------------------------*/
/*! \file
\brief three dimensional total Lagrange truss element

\level 3


*/
/*---------------------------------------------------------------------------*/

#include "4C_truss3.hpp"

#include "4C_discretization_fem_general_node.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::Truss3Type Discret::ELEMENTS::Truss3Type::instance_;

Discret::ELEMENTS::Truss3Type& Discret::ELEMENTS::Truss3Type::Instance() { return instance_; }

Core::Communication::ParObject* Discret::ELEMENTS::Truss3Type::Create(const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::Truss3(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Truss3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "TRUSS3")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Truss3(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Truss3Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Truss3(id, owner));
  return ele;
}


void Discret::ELEMENTS::Truss3Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::Truss3Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void Discret::ELEMENTS::Truss3Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["TRUSS3"];

  defs["LINE2"] = Input::LineDefinition::Builder()
                      .AddIntVector("LINE2", 2)
                      .AddNamedInt("MAT")
                      .AddNamedDouble("CROSS")
                      .AddNamedString("KINEM")
                      .Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 08/08|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Truss3::Truss3(int id, int owner)
    : Core::Elements::Element(id, owner),
      crosssec_(0.0),
      eint_(0.0),
      lrefe_(0.0),
      gaussrule_(Core::FE::GaussRule1D::line_2point),
      diff_disp_ref_(Core::LinAlg::Matrix<1, 3>(true)),
      interface_ptr_(Teuchos::null),
      isinit_(false),
      jacobimass_(),
      jacobinode_(),
      kintype_(KinematicType::tr3_totlag),
      material_(0),
      x_()
{
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 08/08|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Truss3::Truss3(const Discret::ELEMENTS::Truss3& old)
    : Core::Elements::Element(old),
      crosssec_(old.crosssec_),
      eint_(old.eint_),
      lrefe_(old.lrefe_),
      gaussrule_(old.gaussrule_),
      diff_disp_ref_(old.diff_disp_ref_),
      interface_ptr_(old.interface_ptr_),
      isinit_(old.isinit_),
      jacobimass_(old.jacobimass_),
      jacobinode_(old.jacobinode_),
      kintype_(old.kintype_),
      material_(old.material_),
      x_(old.x_)

{
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Truss3 and return pointer to it (public) |
 |                                                            cyron 08/08|
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Truss3::Clone() const
{
  auto* newelement = new Discret::ELEMENTS::Truss3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |(public)                                                   cyron 08/08|
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Truss3::Shape() const { return Core::FE::CellType::line2; }


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 08/08|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Truss3::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);
  AddtoPack(data, isinit_);
  AddtoPack<6, 1>(data, x_);
  AddtoPack<1, 3>(data, diff_disp_ref_);
  AddtoPack(data, material_);
  AddtoPack(data, lrefe_);
  AddtoPack(data, jacobimass_);
  AddtoPack(data, jacobinode_);
  AddtoPack(data, crosssec_);
  AddtoPack(data, gaussrule_);
  AddtoPack(data, kintype_);
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 08/08|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Truss3::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);
  isinit_ = ExtractInt(position, data);
  ExtractfromPack<6, 1>(position, data, x_);
  ExtractfromPack<1, 3>(position, data, diff_disp_ref_);
  ExtractfromPack(position, data, material_);
  ExtractfromPack(position, data, lrefe_);
  ExtractfromPack(position, data, jacobimass_);
  ExtractfromPack(position, data, jacobinode_);
  ExtractfromPack(position, data, crosssec_);
  ExtractfromPack(position, data, gaussrule_);
  // kinematic type
  kintype_ = static_cast<KinematicType>(ExtractInt(position, data));

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                              cyron 08/08|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Truss3::Lines()
{
  return {Teuchos::rcpFromRef(*this)};
}

/*----------------------------------------------------------------------*
 |determine Gauss rule from required type of integration                |
 |                                                   (public)cyron 09/09|
 *----------------------------------------------------------------------*/
Core::FE::GaussRule1D Discret::ELEMENTS::Truss3::my_gauss_rule(
    int nnode, IntegrationType integrationtype)
{
  Core::FE::GaussRule1D gaussrule = Core::FE::GaussRule1D::undefined;

  switch (nnode)
  {
    case 2:
    {
      switch (integrationtype)
      {
        case gaussexactintegration:
        {
          gaussrule = Core::FE::GaussRule1D::line_2point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule = Core::FE::GaussRule1D::line_1point;
          break;
        }
        case lobattointegration:
        {
          gaussrule = Core::FE::GaussRule1D::line_lobatto2point;
          break;
        }
        default:
          FOUR_C_THROW("unknown type of integration");
      }
      break;
    }
    case 3:
    {
      switch (integrationtype)
      {
        case gaussexactintegration:
        {
          gaussrule = Core::FE::GaussRule1D::line_3point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule = Core::FE::GaussRule1D::line_2point;
          break;
        }
        case lobattointegration:
        {
          gaussrule = Core::FE::GaussRule1D::line_lobatto3point;
          break;
        }
        default:
          FOUR_C_THROW("unknown type of integration");
      }
      break;
    }
    case 4:
    {
      switch (integrationtype)
      {
        case gaussexactintegration:
        {
          gaussrule = Core::FE::GaussRule1D::line_4point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule = Core::FE::GaussRule1D::line_3point;
          break;
        }
        default:
          FOUR_C_THROW("unknown type of integration");
      }
      break;
    }
    case 5:
    {
      switch (integrationtype)
      {
        case gaussexactintegration:
        {
          gaussrule = Core::FE::GaussRule1D::line_5point;
          break;
        }
        case gaussunderintegration:
        {
          gaussrule = Core::FE::GaussRule1D::line_4point;
          break;
        }
        default:
          FOUR_C_THROW("unknown type of integration");
      }
      break;
    }
    default:
      FOUR_C_THROW("Only Line2, Line3, Line4 and Line5 Elements implemented.");
  }

  return gaussrule;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Truss3::set_up_reference_geometry(const std::vector<double>& xrefe)
{
  if (!isinit_)
  {
    // setting reference coordinates
    for (int i = 0; i < 6; ++i) x_(i) = xrefe[i];

    // length in reference configuration
    lrefe_ = std::sqrt((x_(3) - x_(0)) * (x_(3) - x_(0)) + (x_(4) - x_(1)) * (x_(4) - x_(1)) +
                       (x_(5) - x_(2)) * (x_(5) - x_(2)));

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
void Discret::ELEMENTS::Truss3::scale_reference_length(double scalefac)
{
  // scale length in reference configuration
  x_(3) = x_(0) + (scalefac * (x_(3) - x_(0)));
  x_(4) = x_(1) + (scalefac * (x_(4) - x_(1)));
  x_(5) = x_(2) + (scalefac * (x_(5) - x_(2)));

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
void Discret::ELEMENTS::Truss3::LocationVector(
    const Discretization& dis, LocationArray& la, bool doDirichlet) const
{
  const int numnode = num_node();
  const Core::Nodes::Node* const* nodes = Nodes();

  la.Clear();

  // we need to look at all DofSets of our discretization
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
        const Core::Nodes::Node* node = nodes[i];

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
          Core::Conditions::Condition* dirich = node->GetCondition("Dirichlet");
          if (dirich)
          {
            if (dirich->Type() != Core::Conditions::PointDirichlet &&
                dirich->Type() != Core::Conditions::LineDirichlet &&
                dirich->Type() != Core::Conditions::SurfaceDirichlet &&
                dirich->Type() != Core::Conditions::VolumeDirichlet)
              FOUR_C_THROW("condition with name Dirichlet is not of type Dirichlet");
            flag = &dirich->parameters().Get<std::vector<int>>("onoff");
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
      Core::Conditions::Condition* dirich = GetCondition("Dirichlet");
      if (dirich)
      {
        if (dirich->Type() != Core::Conditions::PointDirichlet &&
            dirich->Type() != Core::Conditions::LineDirichlet &&
            dirich->Type() != Core::Conditions::SurfaceDirichlet &&
            dirich->Type() != Core::Conditions::VolumeDirichlet)
          FOUR_C_THROW("condition with name Dirichlet is not of type Dirichlet");
        flag = &dirich->parameters().Get<std::vector<int>>("onoff");
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
int Discret::ELEMENTS::Truss3Type::Initialize(Discret::Discretization& dis)
{
  // reference node positions
  std::vector<double> xrefe;

  // reference nodal tangent positions
  Core::LinAlg::Matrix<3, 1> trefNodeAux(true);
  // resize vectors for the number of coordinates we need to store
  xrefe.resize(3 * 2);

  // setting beam reference director correctly
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    // in case that current element is not a truss3 element there is nothing to do and we go back
    // to the head of the loop
    if (dis.lColElement(i)->ElementType() != *this) continue;

    // if we get so far current element is a truss3 element and  we get a pointer at it
    auto* currele = dynamic_cast<Discret::ELEMENTS::Truss3*>(dis.lColElement(i));
    if (!currele) FOUR_C_THROW("cast to Truss3* failed");

    // getting element's nodal coordinates and treating them as reference configuration
    if (currele->Nodes()[0] == nullptr || currele->Nodes()[1] == nullptr)
      FOUR_C_THROW("Cannot get nodes in order to compute reference configuration'");
    else
    {
      for (int k = 0; k < 2; k++)  // element has two nodes
        for (int l = 0; l < 3; l++) xrefe[k * 3 + l] = currele->Nodes()[k]->X()[l];
    }

    currele->set_up_reference_geometry(xrefe);
  }

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Truss3::set_params_interface_ptr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
  {
    interface_ptr_ = Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>(
        p.get<Teuchos::RCP<Core::Elements::ParamsInterface>>("interface"));
  }
  else
    interface_ptr_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::ParamsInterface> Discret::ELEMENTS::Truss3::ParamsInterfacePtr()
{
  return interface_ptr_;
}

FOUR_C_NAMESPACE_CLOSE
