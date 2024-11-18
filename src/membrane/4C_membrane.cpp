// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_membrane.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN


std::shared_ptr<Core::Elements::Element> Discret::Elements::MembraneLine2Type::create(
    const int id, const int owner)
{
  // return Teuchos::rcp( new MembraneLine( id, owner ) );
  return nullptr;
}

std::shared_ptr<Core::Elements::Element> Discret::Elements::MembraneLine3Type::create(
    const int id, const int owner)
{
  // return Teuchos::rcp( new MembraneLine( id, owner ) );
  return nullptr;
}

/*-----------------------------------------------------------------------------*
 |  constructor (public)                                          fbraeu 06/16 |
 |  id          (in)  this element's global id                                 |
 *-----------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::Elements::Membrane<distype>::Membrane(int id, int owner)
    : Core::Elements::Element(id, owner),
      thickness_(0.0),
      cur_thickness_(0),
      planetype_(plane_stress),
      intpoints_(Core::FE::GaussRule2D::tri_3point)
{
  switch (distype)
  {
    case Core::FE::CellType::tri3:
    {
      Core::FE::GaussRule2D gaussrule = Core::FE::GaussRule2D::tri_3point;
      // get gauss integration points
      intpoints_ = Core::FE::IntegrationPoints2D(gaussrule);
      break;
    }
    case Core::FE::CellType::tri6:
    {
      Core::FE::GaussRule2D gaussrule = Core::FE::GaussRule2D::tri_6point;
      // get gauss integration points
      intpoints_ = Core::FE::IntegrationPoints2D(gaussrule);
      break;
    }
    case Core::FE::CellType::quad4:
    {
      Core::FE::GaussRule2D gaussrule = Core::FE::GaussRule2D::quad_4point;
      // get gauss integration points
      intpoints_ = Core::FE::IntegrationPoints2D(gaussrule);
      break;
    }
    case Core::FE::CellType::quad9:
    {
      Core::FE::GaussRule2D gaussrule = Core::FE::GaussRule2D::quad_9point;
      // get gauss integration points
      intpoints_ = Core::FE::IntegrationPoints2D(gaussrule);
      break;
    }
    default:
      FOUR_C_THROW("shape type unknown!\n");
      break;
  }
  cur_thickness_.resize(intpoints_.nquad, thickness_);
  return;
}

/*-----------------------------------------------------------------------------*
 |  copy-constructor (public)                                     fbraeu 06/16 |
 |  id               (in)  this element's global id                            |
 *-----------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::Elements::Membrane<distype>::Membrane(const Discret::Elements::Membrane<distype>& old)
    : Core::Elements::Element(old),
      thickness_(old.thickness_),
      cur_thickness_(old.cur_thickness_),
      planetype_(old.planetype_),
      intpoints_(old.intpoints_)
{
  cur_thickness_.resize(intpoints_.nquad, thickness_);
  return;
}

/*------------------------------------------------------------------------*
 |  Deep copy this instance of Membrane and return pointer to it (public) |
 |                                                           fbraeu 06/16 |
 *------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Core::Elements::Element* Discret::Elements::Membrane<distype>::clone() const
{
  Discret::Elements::Membrane<distype>* newelement =
      new Discret::Elements::Membrane<distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Core::FE::CellType Discret::Elements::Membrane<distype>::shape() const
{
  return distype;
}

/*----------------------------------------------------------------------*
 |  Return number of lines of this element                     (public) |
 |                                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::Elements::Membrane<distype>::num_line() const
{
  return Core::FE::get_number_of_element_lines(distype);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::Elements::Membrane<distype>::pack(Core::Communication::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // add base class Element
  Element::pack(data);

  // thickness_
  add_to_pack(data, thickness_);

  // current thickness_
  add_to_pack(data, cur_thickness_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::Elements::Membrane<distype>::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class Element
  Element::unpack(buffer);
  // thickness_
  extract_from_pack(buffer, thickness_);
  // current thickness_
  extract_from_pack(buffer, cur_thickness_);


  return;
}


/*----------------------------------------------------------------------*
 |  return solid material (public)                         sfuchs 05/17 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
std::shared_ptr<Mat::So3Material> Discret::Elements::Membrane<distype>::solid_material(
    int nummat) const
{
  return std::dynamic_pointer_cast<Mat::So3Material>(Core::Elements::Element::material(nummat));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::Elements::Membrane<distype>::set_params_interface_ptr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
    interface_ptr_ = p.get<std::shared_ptr<Core::Elements::ParamsInterface>>("interface");
  else
    interface_ptr_ = nullptr;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
std::shared_ptr<Core::Elements::ParamsInterface>
Discret::Elements::Membrane<distype>::params_interface_ptr()
{
  return interface_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Solid::Elements::ParamsInterface& Discret::Elements::Membrane<distype>::str_params_interface()
{
  if (not is_params_interface()) FOUR_C_THROW("The interface ptr is not set!");
  return *(std::dynamic_pointer_cast<Solid::Elements::ParamsInterface>(interface_ptr_));
}

/*----------------------------------------------------------------------*
 |  print this element (public)                            fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::Elements::Membrane<distype>::print(std::ostream& os) const
{
  os << "Membrane ";
  os << " discretization type: " << Core::FE::cell_type_to_string(distype).c_str();
  Element::print(os);
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                           fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
std::vector<std::shared_ptr<Core::Elements::Element>> Discret::Elements::Membrane<distype>::lines()
{
  return Core::Communication::element_boundary_factory<MembraneLine<distype>, Membrane<distype>>(
      Core::Communication::buildLines, *this);
}

/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                        fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
std::vector<std::shared_ptr<Core::Elements::Element>>
Discret::Elements::Membrane<distype>::surfaces()
{
  return {Core::Utils::shared_ptr_from_ref(*this)};
}

template class Discret::Elements::Membrane<Core::FE::CellType::tri3>;
template class Discret::Elements::Membrane<Core::FE::CellType::tri6>;
template class Discret::Elements::Membrane<Core::FE::CellType::quad4>;
template class Discret::Elements::Membrane<Core::FE::CellType::quad9>;

FOUR_C_NAMESPACE_CLOSE
