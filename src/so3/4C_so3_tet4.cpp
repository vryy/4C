/*----------------------------------------------------------------------*/
/*! \file

\brief Solid Tet4 element

\level 1


*----------------------------------------------------------------------*/

#include "4C_so3_tet4.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_discretization_fem_general_fiber_node.hpp"
#include "4C_discretization_fem_general_fiber_node_holder.hpp"
#include "4C_discretization_fem_general_fiber_node_utils.hpp"
#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_prestress.hpp"
#include "4C_so3_prestress_service.hpp"
#include "4C_so3_surface.hpp"
#include "4C_so3_utils.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN


DRT::ELEMENTS::SoTet4Type DRT::ELEMENTS::SoTet4Type::instance_;

DRT::ELEMENTS::SoTet4Type& DRT::ELEMENTS::SoTet4Type::Instance() { return instance_; }

//------------------------------------------------------------------------
CORE::COMM::ParObject* DRT::ELEMENTS::SoTet4Type::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::SoTet4(-1, -1);
  object->Unpack(data);
  return object;
}


//------------------------------------------------------------------------
Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoTet4Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::SoTet4(id, owner));
    return ele;
  }
  return Teuchos::null;
}


//------------------------------------------------------------------------
Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::SoTet4Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::SoTet4(id, owner));
  return ele;
}


//------------------------------------------------------------------------
void DRT::ELEMENTS::SoTet4Type::nodal_block_information(
    CORE::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

//------------------------------------------------------------------------
CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::SoTet4Type::ComputeNullSpace(
    CORE::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

//------------------------------------------------------------------------
void DRT::ELEMENTS::SoTet4Type::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["TET4"] = INPUT::LineDefinition::Builder()
                     .AddIntVector("TET4", 4)
                     .AddNamedInt("MAT")
                     .AddNamedString("KINEM")
                     .add_optional_named_double_vector("RAD", 3)
                     .add_optional_named_double_vector("AXI", 3)
                     .add_optional_named_double_vector("CIR", 3)
                     .add_optional_named_double_vector("FIBER1", 3)
                     .add_optional_named_double_vector("FIBER2", 3)
                     .add_optional_named_double_vector("FIBER3", 3)
                     .add_optional_named_double("GROWTHTRIG")
                     .Build();
}

/*----------------------------------------------------------------------***
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoTet4::SoTet4(int id, int owner)
    : SoBase(id, owner),
      // material_(0),
      V_(-1.0),
      pstype_(INPAR::STR::PreStress::none),
      pstime_(0.0),
      time_(0.0)
{
  Teuchos::RCP<const Teuchos::ParameterList> params =
      GLOBAL::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    pstype_ = PRESTRESS::GetType();
    pstime_ = PRESTRESS::GetPrestressTime();

    DRT::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
        GLOBAL::Problem::Instance()->structural_dynamic_params(), get_element_type_string());
  }
  if (PRESTRESS::IsMulf(pstype_))
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOTET4, NUMGPT_SOTET4, true));
}

/*----------------------------------------------------------------------***
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoTet4::SoTet4(const DRT::ELEMENTS::SoTet4& old)
    : SoBase(old),
      // material_(old.material_),
      V_(old.V_),
      pstype_(old.pstype_),
      pstime_(old.pstime_),
      time_(old.time_)
{
  if (PRESTRESS::IsMulf(pstype_))
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(*(old.prestress_)));
}

/*----------------------------------------------------------------------***
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
CORE::Elements::Element* DRT::ELEMENTS::SoTet4::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::SoTet4(*this);
  return newelement;
}

/*----------------------------------------------------------------------***
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::SoTet4::Shape() const { return CORE::FE::CellType::tet4; }

/*----------------------------------------------------------------------***
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoTet4::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  SoBase::Pack(data);
  // ngp_
  // AddtoPack(data,ngp_,3*sizeof(int));
  // material_
  // AddtoPack(data,material_);

  // V_
  AddtoPack(data, V_);

  // Pack prestress
  AddtoPack(data, static_cast<int>(pstype_));
  AddtoPack(data, pstime_);
  AddtoPack(data, time_);
  if (PRESTRESS::IsMulf(pstype_))
  {
    CORE::COMM::ParObject::AddtoPack(data, *prestress_);
  }
}


/*----------------------------------------------------------------------***
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoTet4::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  SoBase::Unpack(basedata);
  // ngp_
  // ExtractfromPack(position,data,ngp_,3*sizeof(int));
  // material_
  // ExtractfromPack(position,data,material_);
  // V_
  ExtractfromPack(position, data, V_);

  // Extract prestress
  pstype_ = static_cast<INPAR::STR::PreStress>(ExtractInt(position, data));
  ExtractfromPack(position, data, pstime_);
  ExtractfromPack(position, data, time_);
  if (PRESTRESS::IsMulf(pstype_))
  {
    std::vector<char> tmpprestress(0);
    ExtractfromPack(position, data, tmpprestress);
    if (prestress_ == Teuchos::null)
      prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOTET4, NUMGPT_SOTET4, true));
    prestress_->Unpack(tmpprestress);
  }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}


/*----------------------------------------------------------------------***
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoTet4::Print(std::ostream& os) const
{
  os << "So_tet4 ";
  Element::Print(os);
  std::cout << std::endl;
  return;
}

/*====================================================================*/
/* 4-node tetrahedra node topology*/
/*--------------------------------------------------------------------*/
/* parameter coordinates (ksi1, ksi2, ksi3) of nodes
 * of a common tetrahedron [0,1]x[0,1]x[0,1]
 *  4-node hexahedron: node 0,1,...,3
 *
 * -----------------------
 *- this is the numbering used in GiD & EXODUS!!
 *      3-
 *      |\ ---
 *      |  \    ---
 *      |    \      ---
 *      |      \        -2
 *      |        \       /\
 *      |          \   /   \
 *      |            X      \
 *      |          /   \     \
 *      |        /       \    \
 *      |      /           \   \
 *      |    /               \  \
 *      |  /                   \ \
 *      |/                       \\
 *      0--------------------------1
 */
/*====================================================================*/


/*----------------------------------------------------------------------*
|  get vector of surfaces (public)                             maf 04/07|
|  surface normals always point outward                                 |
*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::SoTet4::Surfaces()
{
  return CORE::COMM::ElementBoundaryFactory<StructuralSurface, CORE::Elements::Element>(
      CORE::COMM::buildSurfaces, *this);
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
std::vector<double> DRT::ELEMENTS::SoTet4::element_center_refe_coords()
{
  // update element geometry
  CORE::Nodes::Node** nodes = Nodes();
  CORE::LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xrefe;  // material coord. of element
  for (int i = 0; i < NUMNOD_SOTET4; ++i)
  {
    const auto& x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }
  const CORE::FE::CellType distype = Shape();
  CORE::LINALG::Matrix<NUMNOD_SOTET4, 1> funct;
  // Centroid of a tet with (0,1)(0,1)(0,1) is (0.25, 0.25, 0.25)
  CORE::FE::shape_function_3D(funct, 0.25, 0.25, 0.25, distype);
  CORE::LINALG::Matrix<1, NUMDIM_SOTET4> midpoint;
  // midpoint.Multiply('T','N',1.0,funct,xrefe,0.0);
  midpoint.MultiplyTN(funct, xrefe);
  std::vector<double> centercoords(3);
  centercoords[0] = midpoint(0, 0);
  centercoords[1] = midpoint(0, 1);
  centercoords[2] = midpoint(0, 2);
  return centercoords;
}

/*----------------------------------------------------------------------***++
 |  get vector of lines (public)                               maf 04/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::SoTet4::Lines()
{
  return CORE::COMM::ElementBoundaryFactory<StructuralLine, CORE::Elements::Element>(
      CORE::COMM::buildLines, *this);
}

/*----------------------------------------------------------------------*
 |  Return names of visualization data (public)                 st 01/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoTet4::VisNames(std::map<std::string, int>& names)
{
  SolidMaterial()->VisNames(names);

  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                          st 01/10|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::SoTet4::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (CORE::Elements::Element::VisData(name, data)) return true;

  return SolidMaterial()->VisData(name, data, NUMGPT_SOTET4, this->Id());
}

/*----------------------------------------------------------------------*
 |  Call post setup routine of the materials                            |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoTet4::material_post_setup(Teuchos::ParameterList& params)
{
  if (CORE::Nodes::HaveNodalFibers<CORE::FE::CellType::tet4>(Nodes()))
  {
    // This element has fiber nodes.
    // Interpolate fibers to the Gauss points and pass them to the material

    // Get shape functions
    static const std::vector<CORE::LINALG::Matrix<NUMNOD_SOTET4, 1>> shapefcts =
        so_tet4_1gp_shapefcts();

    // add fibers to the ParameterList
    // ParameterList does not allow to store a std::vector, so we have to add every gp fiber
    // with a separate key. To keep it clean, It is added to a sublist.
    CORE::Nodes::NodalFiberHolder fiberHolder;

    // Do the interpolation
    CORE::Nodes::ProjectFibersToGaussPoints<CORE::FE::CellType::tet4>(
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
