// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fluid_ele.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fluid_ele_nullspace.hpp"
#include "4C_fluid_ele_tds.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::Elements::FluidType Discret::Elements::FluidType::instance_;

Discret::Elements::FluidType& Discret::Elements::FluidType::instance() { return instance_; }

Core::Communication::ParObject* Discret::Elements::FluidType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Discret::Elements::Fluid* object = new Discret::Elements::Fluid(-1, -1);
  object->unpack(buffer);
  return object;
}


std::shared_ptr<Core::Elements::Element> Discret::Elements::FluidType::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "FLUID")
  {
    return std::make_shared<Discret::Elements::Fluid>(id, owner);
  }
  else if (eletype == "FLUID2" || eletype == "FLUID3")
  {
    FOUR_C_THROW("Fluid element types FLUID2 and FLUID3 are no longer in use. Switch to FLUID.");
  }
  return nullptr;
}


std::shared_ptr<Core::Elements::Element> Discret::Elements::FluidType::create(
    const int id, const int owner)
{
  return std::make_shared<Discret::Elements::Fluid>(id, owner);
}


void Discret::Elements::FluidType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = dwele->num_dof_per_node(*(dwele->nodes()[0]));
  dimns = numdf;
  nv = numdf - 1;
  np = 1;
}


Core::LinAlg::SerialDenseMatrix Discret::Elements::FluidType::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return FLD::compute_fluid_null_space(node, numdof, dimnsp);
}

void Discret::Elements::FluidType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defsgeneral = definitions["FLUID"];

  defsgeneral["HEX8"] = Input::LineDefinition::Builder()
                            .add_int_vector("HEX8", 8)
                            .add_named_int("MAT")
                            .add_named_string("NA")
                            .build();

  defsgeneral["HEX20"] = Input::LineDefinition::Builder()
                             .add_int_vector("HEX20", 20)
                             .add_named_int("MAT")
                             .add_named_string("NA")
                             .build();

  defsgeneral["HEX27"] = Input::LineDefinition::Builder()
                             .add_int_vector("HEX27", 27)
                             .add_named_int("MAT")
                             .add_named_string("NA")
                             .build();

  defsgeneral["TET4"] = Input::LineDefinition::Builder()
                            .add_int_vector("TET4", 4)
                            .add_named_int("MAT")
                            .add_named_string("NA")
                            .build();

  defsgeneral["TET10"] = Input::LineDefinition::Builder()
                             .add_int_vector("TET10", 10)
                             .add_named_int("MAT")
                             .add_named_string("NA")
                             .build();

  defsgeneral["WEDGE6"] = Input::LineDefinition::Builder()
                              .add_int_vector("WEDGE6", 6)
                              .add_named_int("MAT")
                              .add_named_string("NA")
                              .build();

  defsgeneral["WEDGE15"] = Input::LineDefinition::Builder()
                               .add_int_vector("WEDGE15", 15)
                               .add_named_int("MAT")
                               .add_named_string("NA")
                               .build();

  defsgeneral["PYRAMID5"] = Input::LineDefinition::Builder()
                                .add_int_vector("PYRAMID5", 5)
                                .add_named_int("MAT")
                                .add_named_string("NA")
                                .build();

  defsgeneral["NURBS8"] = Input::LineDefinition::Builder()
                              .add_int_vector("NURBS8", 8)
                              .add_named_int("MAT")
                              .add_named_string("NA")
                              .build();

  defsgeneral["NURBS27"] = Input::LineDefinition::Builder()
                               .add_int_vector("NURBS27", 27)
                               .add_named_int("MAT")
                               .add_named_string("NA")
                               .build();

  // 2D elements
  defsgeneral["QUAD4"] = Input::LineDefinition::Builder()
                             .add_int_vector("QUAD4", 4)
                             .add_named_int("MAT")
                             .add_named_string("NA")
                             .build();

  defsgeneral["QUAD8"] = Input::LineDefinition::Builder()
                             .add_int_vector("QUAD8", 8)
                             .add_named_int("MAT")
                             .add_named_string("NA")
                             .build();

  defsgeneral["QUAD9"] = Input::LineDefinition::Builder()
                             .add_int_vector("QUAD9", 9)
                             .add_named_int("MAT")
                             .add_named_string("NA")
                             .build();

  defsgeneral["TRI3"] = Input::LineDefinition::Builder()
                            .add_int_vector("TRI3", 3)
                            .add_named_int("MAT")
                            .add_named_string("NA")
                            .build();

  defsgeneral["TRI6"] = Input::LineDefinition::Builder()
                            .add_int_vector("TRI6", 6)
                            .add_named_int("MAT")
                            .add_named_string("NA")
                            .build();

  defsgeneral["NURBS4"] = Input::LineDefinition::Builder()
                              .add_int_vector("NURBS4", 4)
                              .add_named_int("MAT")
                              .add_named_string("NA")
                              .build();

  defsgeneral["NURBS9"] = Input::LineDefinition::Builder()
                              .add_int_vector("NURBS9", 9)
                              .add_named_int("MAT")
                              .add_named_string("NA")
                              .build();
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 02/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::Elements::Fluid::Fluid(int id, int owner)
    : Core::Elements::Element(id, owner), is_ale_(false)
{
  distype_ = Core::FE::CellType::dis_none;
  tds_ = nullptr;
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 02/08|
 *----------------------------------------------------------------------*/
Discret::Elements::Fluid::Fluid(const Discret::Elements::Fluid& old)
    : Core::Elements::Element(old), distype_(old.distype_), is_ale_(old.is_ale_)
{
  tds_ = nullptr;
  if (old.tds_ != nullptr)
    FOUR_C_THROW("clone() method for deep copying tds_ not yet implemented!");
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Fluid and return pointer to it (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::Elements::Fluid::clone() const
{
  Discret::Elements::Fluid* newelement = new Discret::Elements::Fluid(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
void Discret::Elements::Fluid::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class Element
  Element::pack(data);
  // is_ale_
  add_to_pack(data, is_ale_);
  // Discretisation type
  add_to_pack(data, distype_);

  // time-dependent subgrid scales
  bool is_tds(false);
  if (tds_ != nullptr)
  {
    is_tds = true;
    add_to_pack(data, is_tds);
    tds_->pack(data);
  }
  else
  {
    add_to_pack(data, is_tds);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
void Discret::Elements::Fluid::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer base_buffer(basedata);
  Element::unpack(base_buffer);
  // is_ale_
  extract_from_pack(buffer, is_ale_);
  // distype
  extract_from_pack(buffer, distype_);

  // time-dependent subgrid scales
  bool is_tds;
  extract_from_pack(buffer, is_tds);
  if (is_tds)
  {
    tds_ = std::make_shared<FLD::TDSEleData>();
    std::vector<char> pbtest;
    extract_from_pack(buffer, pbtest);
    if (pbtest.size() == 0) FOUR_C_THROW("Seems no TDS data available");
    Core::Communication::UnpackBuffer pbtest_buffer(pbtest);
    tds_->unpack(pbtest_buffer);
  }
  else
    tds_ = nullptr;

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 02/08|
 *----------------------------------------------------------------------*/
void Discret::Elements::Fluid::print(std::ostream& os) const
{
  os << "Fluid ";
  Element::print(os);
  // cout << endl;
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of lines              (public)                 ae  02/010|
 *----------------------------------------------------------------------*/
std::vector<std::shared_ptr<Core::Elements::Element>> Discret::Elements::Fluid::lines()
{
  return Core::Communication::get_element_lines<FluidBoundary, Fluid>(*this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          ehrl  02/10|
 *----------------------------------------------------------------------*/
std::vector<std::shared_ptr<Core::Elements::Element>> Discret::Elements::Fluid::surfaces()
{
  return Core::Communication::get_element_surfaces<FluidBoundary, Fluid>(*this);
}


/*----------------------------------------------------------------------*
 |  get face element (public)                               schott 03/12|
 *----------------------------------------------------------------------*/
std::shared_ptr<Core::Elements::Element> Discret::Elements::Fluid::create_face_element(
    Core::Elements::Element* parent_slave,  //!< parent slave fluid3 element
    int nnode,                              //!< number of surface nodes
    const int* nodeids,                     //!< node ids of surface element
    Core::Nodes::Node** nodes,              //!< nodes of surface element
    const int lsurface_master,              //!< local surface number w.r.t master parent element
    const int lsurface_slave,               //!< local surface number w.r.t slave parent element
    const std::vector<int>& localtrafomap   //! local trafo map
)
{
  // dynamic cast for slave parent element
  Discret::Elements::Fluid* slave_pele = dynamic_cast<Discret::Elements::Fluid*>(parent_slave);


  // insert both parent elements
  return Core::Communication::element_int_face_factory<FluidIntFace, Fluid>(
      -1,               //!< internal face element id
      -1,               //!< owner of internal face element
      nnode,            //!< number of surface nodes
      nodeids,          //!< node ids of surface element
      nodes,            //!< nodes of surface element
      this,             //!< master parent element
      slave_pele,       //!< slave parent element
      lsurface_master,  //!< local surface number w.r.t master parent element
      lsurface_slave,   //!< local surface number w.r.t slave parent element
      localtrafomap     //!< local trafo map
  );
}


/*----------------------------------------------------------------------*
 |  activate time dependent subgrid scales (public)      gamnitzer 05/10|
 *----------------------------------------------------------------------*/
void Discret::Elements::Fluid::activate_tds(
    int nquad, int nsd, double** saccn, double** sveln, double** svelnp)
{
  if (tds_ == nullptr) tds_ = std::make_shared<FLD::TDSEleData>();

  tds_->activate_tds(nquad, nsd, saccn, sveln, svelnp);
}

FOUR_C_NAMESPACE_CLOSE
