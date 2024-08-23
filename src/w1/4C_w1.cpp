/*----------------------------------------------------------------------------*/
/*! \file
\brief wall1 element.

\level 1


*/
/*---------------------------------------------------------------------------*/

#include "4C_w1.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::Wall1Type Discret::ELEMENTS::Wall1Type::instance_;

Discret::ELEMENTS::Wall1Type& Discret::ELEMENTS::Wall1Type::instance() { return instance_; }

Core::Communication::ParObject* Discret::ELEMENTS::Wall1Type::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Discret::ELEMENTS::Wall1* object = new Discret::ELEMENTS::Wall1(-1, -1);
  object->unpack(buffer);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Wall1Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "WALL")
  {
    if (eledistype != "NURBS4" and eledistype != "NURBS9")
    {
      return Teuchos::rcp(new Discret::ELEMENTS::Wall1(id, owner));
    }
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Wall1Type::create(
    const int id, const int owner)
{
  return Teuchos::rcp(new Discret::ELEMENTS::Wall1(id, owner));
}


void Discret::ELEMENTS::Wall1Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 2;
  dimns = 3;
  nv = 2;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::Wall1Type::compute_null_space(
    Core::Nodes::Node& node, const double* x0, int const numdof, int const dimnsp)
{
  return compute_solid_2d_null_space(node, x0);
}

void Discret::ELEMENTS::Wall1Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["WALL"];

  defs["QUAD4"] = Input::LineDefinition::Builder()
                      .add_int_vector("QUAD4", 4)
                      .add_named_int("MAT")
                      .add_named_string("KINEM")
                      .add_named_string("EAS")
                      .add_named_double("THICK")
                      .add_named_string("STRESS_STRAIN")
                      .add_named_int_vector("GP", 2)
                      .build();

  defs["QUAD8"] = Input::LineDefinition::Builder()
                      .add_int_vector("QUAD8", 8)
                      .add_named_int("MAT")
                      .add_named_string("KINEM")
                      .add_named_string("EAS")
                      .add_named_double("THICK")
                      .add_named_string("STRESS_STRAIN")
                      .add_named_int_vector("GP", 2)
                      .build();

  defs["QUAD9"] = Input::LineDefinition::Builder()
                      .add_int_vector("QUAD9", 9)
                      .add_named_int("MAT")
                      .add_named_string("KINEM")
                      .add_named_string("EAS")
                      .add_named_double("THICK")
                      .add_named_string("STRESS_STRAIN")
                      .add_named_int_vector("GP", 2)
                      .build();

  defs["TRI3"] = Input::LineDefinition::Builder()
                     .add_int_vector("TRI3", 3)
                     .add_named_int("MAT")
                     .add_named_string("KINEM")
                     .add_named_string("EAS")
                     .add_named_double("THICK")
                     .add_named_string("STRESS_STRAIN")
                     .add_named_int_vector("GP", 2)
                     .build();

  defs["TRI6"] = Input::LineDefinition::Builder()
                     .add_int_vector("TRI6", 6)
                     .add_named_int("MAT")
                     .add_named_string("KINEM")
                     .add_named_string("EAS")
                     .add_named_double("THICK")
                     .add_named_string("STRESS_STRAIN")
                     .add_named_int_vector("GP", 2)
                     .build();

  defs["NURBS4"] = Input::LineDefinition::Builder()
                       .add_int_vector("NURBS4", 4)
                       .add_named_int("MAT")
                       .add_named_string("KINEM")
                       .add_named_string("EAS")
                       .add_named_double("THICK")
                       .add_named_string("STRESS_STRAIN")
                       .add_named_int_vector("GP", 2)
                       .build();

  defs["NURBS9"] = Input::LineDefinition::Builder()
                       .add_int_vector("NURBS9", 9)
                       .add_named_int("MAT")
                       .add_named_string("KINEM")
                       .add_named_string("EAS")
                       .add_named_double("THICK")
                       .add_named_string("STRESS_STRAIN")
                       .add_named_int_vector("GP", 2)
                       .build();
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Wall1LineType::create(
    const int id, const int owner)
{
  // return Teuchos::rcp( new Wall1Line( id, owner ) );
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mgit 01/08/|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Wall1::Wall1(int id, int owner)
    : SoBase(id, owner),
      material_(0),
      thickness_(0.0),
      old_step_length_(0.0),
      gaussrule_(Core::FE::GaussRule2D::undefined),
      wtype_(plane_none),
      stresstype_(w1_none),
      iseas_(false),
      eastype_(eas_vague),
      easdata_(EASData()),
      distype_(Core::FE::CellType::dis_none)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mgit 01/08|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Wall1::Wall1(const Discret::ELEMENTS::Wall1& old)
    : SoBase(old),
      material_(old.material_),
      thickness_(old.thickness_),
      old_step_length_(old.old_step_length_),
      gaussrule_(old.gaussrule_),
      wtype_(old.wtype_),
      stresstype_(old.stresstype_),
      iseas_(old.iseas_),
      eastype_(old.eas_vague),
      easdata_(old.easdata_),
      distype_(old.distype_)
{
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Wall1 and return pointer to it (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Wall1::clone() const
{
  Discret::ELEMENTS::Wall1* newelement = new Discret::ELEMENTS::Wall1(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          mgit 04/07 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Wall1::shape() const { return distype_; }


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Wall1::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class Element
  SoBase::pack(data);
  // material_
  add_to_pack(data, material_);
  // thickness
  add_to_pack(data, thickness_);
  // plane strain or plane stress information
  add_to_pack(data, wtype_);
  // gaussrule_
  add_to_pack(data, gaussrule_);
  // stresstype
  add_to_pack(data, stresstype_);
  // eas
  add_to_pack(data, iseas_);
  // eas type
  add_to_pack(data, eastype_);
  // eas data
  pack_eas_data(data);
  // distype
  add_to_pack(data, distype_);
  // line search
  add_to_pack(data, old_step_length_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Wall1::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer basedata_buffer(basedata);
  SoBase::unpack(basedata_buffer);
  // material_
  extract_from_pack(buffer, material_);
  // thickness_
  extract_from_pack(buffer, thickness_);
  // plane strain or plane stress information_
  wtype_ = static_cast<DimensionalReduction>(extract_int(buffer));
  // gaussrule_
  extract_from_pack(buffer, gaussrule_);
  // stresstype_
  stresstype_ = static_cast<StressType>(extract_int(buffer));
  // iseas_
  iseas_ = extract_int(buffer);
  // eastype_
  eastype_ = static_cast<EasType>(extract_int(buffer));
  // easdata_
  unpack_eas_data(buffer);
  // distype_
  distype_ = static_cast<Core::FE::CellType>(extract_int(buffer));
  // line search
  extract_from_pack(buffer, old_step_length_);
  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
  return;
}



/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             mgit 07/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Wall1::lines()
{
  return Core::Communication::element_boundary_factory<Wall1Line, Wall1>(
      Core::Communication::buildLines, *this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          mgit 03/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Wall1::surfaces()
{
  return {Teuchos::rcpFromRef(*this)};
}

/*-----------------------------------------------------------------------------*
| Map plane Green-Lagrange strains to 3d                       mayr.mt 05/2014 |
*-----------------------------------------------------------------------------*/
void Discret::ELEMENTS::Wall1::green_lagrange_plane3d(
    const Core::LinAlg::SerialDenseVector& glplane, Core::LinAlg::Matrix<6, 1>& gl3d)
{
  gl3d(0) = glplane(0);               // E_{11}
  gl3d(1) = glplane(1);               // E_{22}
  gl3d(2) = 0.0;                      // E_{33}
  gl3d(3) = glplane(2) + glplane(3);  // 2*E_{12}=E_{12}+E_{21}
  gl3d(4) = 0.0;                      // 2*E_{23}
  gl3d(5) = 0.0;                      // 2*E_{31}

  return;
}

FOUR_C_NAMESPACE_CLOSE
