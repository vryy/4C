/*----------------------------------------------------------------------*/
/*! \file

\brief Solid Hex8 element

\level 1


*----------------------------------------------------------------------*/

#include "4C_so3_hex8.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_fiber_node.hpp"
#include "4C_fem_general_fiber_node_holder.hpp"
#include "4C_fem_general_fiber_node_utils.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_element_service.hpp"
#include "4C_so3_hex8fbar.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_prestress.hpp"
#include "4C_so3_prestress_service.hpp"
#include "4C_so3_surface.hpp"
#include "4C_so3_utils.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::SoHex8Type Discret::ELEMENTS::SoHex8Type::instance_;

Discret::ELEMENTS::SoHex8Type& Discret::ELEMENTS::SoHex8Type::Instance() { return instance_; }

namespace
{
  const std::string name = Discret::ELEMENTS::SoHex8Type::Instance().Name();
}

Core::Communication::ParObject* Discret::ELEMENTS::SoHex8Type::Create(const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::SoHex8(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::SoHex8(id, owner));
    return ele;
  }

  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::SoHex8Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::SoHex8(id, owner));
  return ele;
}


void Discret::ELEMENTS::SoHex8Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::SoHex8Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void Discret::ELEMENTS::SoHex8Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX8"] = Input::LineDefinition::Builder()
                     .AddIntVector("HEX8", 8)
                     .AddNamedInt("MAT")
                     .AddNamedString("KINEM")
                     .AddNamedString("EAS")
                     .add_optional_named_double_vector("RAD", 3)
                     .add_optional_named_double_vector("AXI", 3)
                     .add_optional_named_double_vector("CIR", 3)
                     .add_optional_named_double_vector("FIBER1", 3)
                     .add_optional_named_double_vector("FIBER2", 3)
                     .add_optional_named_double_vector("FIBER3", 3)
                     .add_optional_named_double("STRENGTH")
                     .add_optional_named_double("GROWTHTRIG")
                     .Build();
}

// initialization of static gauss point rule for the so_hex8 element
const Core::FE::IntPointsAndWeights<NUMDIM_SOH8> Discret::ELEMENTS::SoHex8::gp_rule_(
    Core::FE::IntPointsAndWeights<NUMDIM_SOH8>(
        static_cast<enum Core::FE::GaussRule3D>(GpRuleSoH8::rule)));

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoHex8::SoHex8(int id, int owner)
    : SoBase(id, owner),
      easdata_(EASData()),
      analyticalmaterialtangent_(true),
      pstype_(Inpar::STR::PreStress::none),
      pstime_(0.0),
      time_(0.0),
      old_step_length_(0.0)
{
  eastype_ = soh8_easnone;
  neas_ = 0;
  invJ_.resize(NUMGPT_SOH8, Core::LinAlg::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>(true));
  detJ_.resize(NUMGPT_SOH8, 0.0);

  Teuchos::RCP<const Teuchos::ParameterList> params =
      Global::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    const Teuchos::ParameterList& sdyn = Global::Problem::Instance()->structural_dynamic_params();

    pstype_ = Prestress::GetType();
    pstime_ = Prestress::GetPrestressTime();
    if (Core::UTILS::IntegralValue<int>(sdyn, "MATERIALTANGENT"))
      analyticalmaterialtangent_ = false;
  }
  if (Prestress::IsMulf(pstype_))
    prestress_ = Teuchos::rcp(new Discret::ELEMENTS::PreStress(NUMNOD_SOH8, NUMGPT_SOH8));

  if (Global::Problem::Instance()->GetProblemType() == Core::ProblemType::struct_ale)
  {
    if (kintype_ == Inpar::STR::KinemType::linear)
      FOUR_C_THROW("Structure-Ale approach only for nonlinear kinematics !!!");

    structale_ = true;
  }
  else
    structale_ = false;


  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::SoHex8::SoHex8(const Discret::ELEMENTS::SoHex8& old)
    : SoBase(old),
      eastype_(old.eastype_),
      neas_(old.neas_),
      easdata_(old.easdata_),
      detJ_(old.detJ_),
      analyticalmaterialtangent_(old.analyticalmaterialtangent_),
      pstype_(old.pstype_),
      pstime_(old.pstime_),
      time_(old.time_),
      old_step_length_(old.old_step_length_)
{
  invJ_.resize(old.invJ_.size());
  for (int i = 0; i < (int)invJ_.size(); ++i)
  {
    // can this size be anything but NUMDIM_SOH8 x NUMDIM_SOH8?
    // invJ_[i].Shape(old.invJ_[i].numRows(),old.invJ_[i].numCols());
    invJ_[i] = old.invJ_[i];
  }

  if (Prestress::IsMulf(pstype_))
    prestress_ = Teuchos::rcp(new Discret::ELEMENTS::PreStress(*(old.prestress_)));

  if (Global::Problem::Instance()->GetProblemType() == Core::ProblemType::struct_ale)
  {
    if (kintype_ == Inpar::STR::KinemType::linear)
      FOUR_C_THROW("Structure-Ale approach only for nonlinear kinematics !!!");

    structale_ = true;
  }
  else
    structale_ = false;

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::SoHex8::Clone() const
{
  auto* newelement = new Discret::ELEMENTS::SoHex8(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::SoHex8::Shape() const { return Core::FE::CellType::hex8; }

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  SoBase::Pack(data);
  // eastype_
  AddtoPack(data, eastype_);
  // neas_
  AddtoPack(data, neas_);
  // analyticalmaterialtangent_
  AddtoPack(data, analyticalmaterialtangent_);
  // eas data
  pack_eas_data(data);
  // line search
  AddtoPack(data, old_step_length_);
  // Pack prestress type
  AddtoPack(data, static_cast<int>(pstype_));
  AddtoPack(data, pstime_);
  AddtoPack(data, time_);
  if (Prestress::IsMulf(pstype_))
  {
    Core::Communication::ParObject::AddtoPack(data, *prestress_);
  }

  // detJ_
  AddtoPack(data, detJ_);

  // invJ_
  const auto size = (int)invJ_.size();
  AddtoPack(data, size);
  for (int i = 0; i < size; ++i) AddtoPack(data, invJ_[i]);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  SoBase::Unpack(basedata);
  // eastype_
  eastype_ = static_cast<EASType>(ExtractInt(position, data));
  // neas_
  ExtractfromPack(position, data, neas_);
  // analyticalmaterialtangent_
  analyticalmaterialtangent_ = ExtractInt(position, data);
  // eas data
  unpack_eas_data(position, data);
  // line search
  ExtractfromPack(position, data, old_step_length_);
  // Extract prestress
  pstype_ = static_cast<Inpar::STR::PreStress>(ExtractInt(position, data));
  ExtractfromPack(position, data, pstime_);
  ExtractfromPack(position, data, time_);
  if (Prestress::IsMulf(pstype_))
  {
    std::vector<char> tmpprestress(0);
    ExtractfromPack(position, data, tmpprestress);
    if (prestress_ == Teuchos::null)
    {
      int numgpt = NUMGPT_SOH8;
      // see whether I am actually a So_hex8fbar element
      auto* me = dynamic_cast<Discret::ELEMENTS::SoHex8fbar*>(this);
      if (me) numgpt += 1;  // one more history entry for centroid data in hex8fbar
      prestress_ = Teuchos::rcp(new Discret::ELEMENTS::PreStress(NUMNOD_SOH8, numgpt));
    }
    prestress_->Unpack(tmpprestress);
  }

  // detJ_
  ExtractfromPack(position, data, detJ_);
  // invJ_
  int size = 0;
  ExtractfromPack(position, data, size);
  invJ_.resize(size, Core::LinAlg::Matrix<NUMDIM_SOH8, NUMDIM_SOH8>(true));
  for (int i = 0; i < size; ++i) ExtractfromPack(position, data, invJ_[i]);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8::Print(std::ostream& os) const
{
  os << "So_hex8 ";
  Element::Print(os);
  // std::cout << std::endl;
  return;
}

/*====================================================================*/
/* 8-node hexhedra node topology*/
/*--------------------------------------------------------------------*/
/* parameter coordinates (r,s,t) of nodes
 * of biunit cube [-1,1]x[-1,1]x[-1,1]
 *  8-node hexahedron: node 0,1,...,7
 *                      t
 *                      |
 *             4========|================7
 *           //|        |               /||
 *          // |        |              //||
 *         //  |        |             // ||
 *        //   |        |            //  ||
 *       //    |        |           //   ||
 *      //     |        |          //    ||
 *     //      |        |         //     ||
 *     5=========================6       ||
 *    ||       |        |        ||      ||
 *    ||       |        o--------||---------s
 *    ||       |       /         ||      ||
 *    ||       0------/----------||------3
 *    ||      /      /           ||     //
 *    ||     /      /            ||    //
 *    ||    /      /             ||   //
 *    ||   /      /              ||  //
 *    ||  /      /               || //
 *    || /      r                ||//
 *    ||/                        ||/
 *     1=========================2
 *
 */
/*====================================================================*/

/*----------------------------------------------------------------------*
|  get vector of surfaces (public)                             maf 04/07|
|  surface normals always point outward                                 |
*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::SoHex8::Surfaces()
{
  return Core::Communication::ElementBoundaryFactory<StructuralSurface, Core::Elements::Element>(
      Core::Communication::buildSurfaces, *this);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               maf 04/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::SoHex8::Lines()
{
  return Core::Communication::ElementBoundaryFactory<StructuralLine, Core::Elements::Element>(
      Core::Communication::buildLines, *this);
}

/*----------------------------------------------------------------------*
 |  get location of element center                              jb 08/11|
 *----------------------------------------------------------------------*/
std::vector<double> Discret::ELEMENTS::SoHex8::element_center_refe_coords()
{
  // update element geometry
  Core::Nodes::Node** nodes = Nodes();
  Core::LinAlg::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;  // material coord. of element
  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    const auto& x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }
  const Core::FE::CellType distype = Shape();
  Core::LinAlg::Matrix<NUMNOD_SOH8, 1> funct;
  // Element midpoint at r=s=t=0.0
  Core::FE::shape_function_3D(funct, 0.0, 0.0, 0.0, distype);
  Core::LinAlg::Matrix<1, NUMDIM_SOH8> midpoint;
  // midpoint.Multiply('T','N',1.0,funct,xrefe,0.0);
  midpoint.MultiplyTN(funct, xrefe);
  std::vector<double> centercoords(3);
  centercoords[0] = midpoint(0, 0);
  centercoords[1] = midpoint(0, 1);
  centercoords[2] = midpoint(0, 2);
  return centercoords;
}
/*----------------------------------------------------------------------*
 |  Return names of visualization data (public)                maf 01/08|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoHex8::VisNames(std::map<std::string, int>& names)
{
  SolidMaterial()->VisNames(names);

  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                         maf 01/08|
 *----------------------------------------------------------------------*/
bool Discret::ELEMENTS::SoHex8::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (Core::Elements::Element::VisData(name, data)) return true;

  return SolidMaterial()->VisData(name, data, NUMGPT_SOH8, this->Id());
}

// Compute nodal fibers and call post setup routine of the materials
void Discret::ELEMENTS::SoHex8::material_post_setup(Teuchos::ParameterList& params)
{
  if (Core::Nodes::HaveNodalFibers<Core::FE::CellType::hex8>(Nodes()))
  {
    // This element has fiber nodes.
    // Interpolate fibers to the Gauss points and pass them to the material

    // Get shape functions
    const std::vector<Core::LinAlg::Matrix<NUMNOD_SOH8, 1>> shapefcts = soh8_shapefcts();

    // add fibers to the ParameterList
    // ParameterList does not allow to store a std::vector, so we have to add every gp fiber
    // with a separate key. To keep it clean, It is added to a sublist.
    Core::Nodes::NodalFiberHolder fiberHolder;

    // Do the interpolation
    Core::Nodes::ProjectFibersToGaussPoints<Core::FE::CellType::hex8>(
        Nodes(), shapefcts, fiberHolder);

    params.set("fiberholder", fiberHolder);
  }

  // Call super post setup
  SoBase::material_post_setup(params);

  // Cleanup ParameterList to not carry all fibers the whole simulation
  // do not throw an error if key does not exist.
  params.remove("fiberholder", false);
}

FOUR_C_NAMESPACE_CLOSE
