// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_so3_hex20.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
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

Discret::Elements::SoHex20Type Discret::Elements::SoHex20Type::instance_;

Discret::Elements::SoHex20Type& Discret::Elements::SoHex20Type::instance() { return instance_; }

Core::Communication::ParObject* Discret::Elements::SoHex20Type::create(
    Core::Communication::UnpackBuffer& buffer)
{
  auto* object = new Discret::Elements::SoHex20(-1, -1);
  object->unpack(buffer);
  return object;
}


std::shared_ptr<Core::Elements::Element> Discret::Elements::SoHex20Type::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    std::shared_ptr<Core::Elements::Element> ele =
        std::make_shared<Discret::Elements::SoHex20>(id, owner);
    return ele;
  }
  return nullptr;
}


std::shared_ptr<Core::Elements::Element> Discret::Elements::SoHex20Type::create(
    const int id, const int owner)
{
  std::shared_ptr<Core::Elements::Element> ele =
      std::make_shared<Discret::Elements::SoHex20>(id, owner);
  return ele;
}


void Discret::Elements::SoHex20Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

Core::LinAlg::SerialDenseMatrix Discret::Elements::SoHex20Type::compute_null_space(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return compute_solid_3d_null_space(node, x0);
}

void Discret::Elements::SoHex20Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["HEX20"] = Input::LineDefinition::Builder()
                      .add_int_vector("HEX20", 20)
                      .add_named_int("MAT")
                      .add_named_string("KINEM")
                      .add_optional_named_double_vector("RAD", 3)
                      .add_optional_named_double_vector("AXI", 3)
                      .add_optional_named_double_vector("CIR", 3)
                      .add_optional_named_double_vector("FIBER1", 3)
                      .add_optional_named_double_vector("FIBER2", 3)
                      .add_optional_named_double_vector("FIBER3", 3)
                      .add_optional_named_double("STRENGTH")
                      .add_optional_named_double("GROWTHTRIG")
                      .build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                                       |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::Elements::SoHex20::SoHex20(int id, int owner)
    : SoBase(id, owner), pstype_(Inpar::Solid::PreStress::none), pstime_(0.0), time_(0.0)
{
  invJ_.resize(NUMGPT_SOH20, Core::LinAlg::Matrix<NUMDIM_SOH20, NUMDIM_SOH20>(true));
  detJ_.resize(NUMGPT_SOH20, 0.0);

  std::shared_ptr<const Teuchos::ParameterList> params =
      Global::Problem::instance()->get_parameter_list();
  if (params != nullptr)
  {
    pstype_ = Prestress::get_type();
    pstime_ = Prestress::get_prestress_time();

    Discret::Elements::Utils::throw_error_fd_material_tangent(
        Global::Problem::instance()->structural_dynamic_params(), get_element_type_string());
  }
  if (Prestress::is_mulf(pstype_))
    prestress_ = std::make_shared<Discret::Elements::PreStress>(NUMNOD_SOH20, NUMGPT_SOH20);

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                                  |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::Elements::SoHex20::SoHex20(const Discret::Elements::SoHex20& old)
    : SoBase(old), detJ_(old.detJ_), pstype_(old.pstype_), pstime_(old.pstime_), time_(old.time_)
{
  invJ_.resize(old.invJ_.size());
  for (int i = 0; i < (int)invJ_.size(); ++i)
  {
    // can this size be anything but NUMDIM_SOH20 x NUMDIM_SOH20?
    // invJ_[i].Shape(old.invJ_[i].numRows(),old.invJ_[i].numCols());
    invJ_[i] = old.invJ_[i];
  }

  if (Prestress::is_mulf(pstype_))
    prestress_ = std::make_shared<Discret::Elements::PreStress>(*(old.prestress_));

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::Elements::SoHex20::clone() const
{
  auto* newelement = new Discret::Elements::SoHex20(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::Elements::SoHex20::shape() const { return Core::FE::CellType::hex20; }

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 *----------------------------------------------------------------------*/
void Discret::Elements::SoHex20::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class Element
  SoBase::pack(data);

  // detJ_
  add_to_pack(data, detJ_);

  // invJ_
  const auto size = (int)invJ_.size();
  add_to_pack(data, size);
  for (int i = 0; i < size; ++i) add_to_pack(data, invJ_[i]);
  // Pack prestress
  add_to_pack(data, pstype_);
  add_to_pack(data, pstime_);
  add_to_pack(data, time_);
  if (Prestress::is_mulf(pstype_))
  {
    add_to_pack(data, *prestress_);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 *----------------------------------------------------------------------*/
void Discret::Elements::SoHex20::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer basedata_buffer(basedata);
  SoBase::unpack(basedata_buffer);

  // detJ_
  extract_from_pack(buffer, detJ_);
  // invJ_
  int size = 0;
  extract_from_pack(buffer, size);
  invJ_.resize(size, Core::LinAlg::Matrix<NUMDIM_SOH20, NUMDIM_SOH20>(true));
  for (int i = 0; i < size; ++i) extract_from_pack(buffer, invJ_[i]);

  // Extract prestress
  extract_from_pack(buffer, pstype_);
  extract_from_pack(buffer, pstime_);
  extract_from_pack(buffer, time_);
  if (Prestress::is_mulf(pstype_))
  {
    std::vector<char> tmpprestress(0);
    extract_from_pack(buffer, tmpprestress);
    if (prestress_ == nullptr)
      prestress_ = std::make_shared<Discret::Elements::PreStress>(NUMNOD_SOH20, NUMGPT_SOH20);
    Core::Communication::UnpackBuffer tmpprestress_buffer(tmpprestress);
    prestress_->unpack(tmpprestress_buffer);
  }

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                                         |
 *----------------------------------------------------------------------*/
void Discret::Elements::SoHex20::print(std::ostream& os) const
{
  os << "So_hex20 ";
  Element::print(os);
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
|  get vector of surfaces (public)                                      |
|  surface normals always point outward                                 |
*----------------------------------------------------------------------*/
std::vector<std::shared_ptr<Core::Elements::Element>> Discret::Elements::SoHex20::surfaces()
{
  return Core::Communication::element_boundary_factory<StructuralSurface, Core::Elements::Element>(
      Core::Communication::buildSurfaces, *this);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                                        |
 *----------------------------------------------------------------------*/
std::vector<std::shared_ptr<Core::Elements::Element>> Discret::Elements::SoHex20::lines()
{
  return Core::Communication::element_boundary_factory<StructuralLine, Core::Elements::Element>(
      Core::Communication::buildLines, *this);
}

/*----------------------------------------------------------------------*
 |  Return names of visualization data (public)                         |
 *----------------------------------------------------------------------*/
void Discret::Elements::SoHex20::vis_names(std::map<std::string, int>& names)
{
  solid_material()->vis_names(names);
  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                                  |
 *----------------------------------------------------------------------*/
bool Discret::Elements::SoHex20::vis_data(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (Core::Elements::Element::vis_data(name, data)) return true;

  return solid_material()->vis_data(name, data, NUMGPT_SOH20, this->id());
}

FOUR_C_NAMESPACE_CLOSE
