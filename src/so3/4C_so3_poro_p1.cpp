// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_so3_poro_p1.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_poro_p1_eletypes.hpp"
#include "4C_so3_surface.hpp"

FOUR_C_NAMESPACE_OPEN

template <class So3Ele, Core::FE::CellType distype>
Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::So3PoroP1(int id, int owner)
    : So3Poro<So3Ele, distype>(id, owner), init_porosity_(Teuchos::null), is_init_porosity_(false)
{
}

template <class So3Ele, Core::FE::CellType distype>
Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::So3PoroP1(
    const Discret::ELEMENTS::So3PoroP1<So3Ele, distype>& old)
    : So3Poro<So3Ele, distype>(old),
      init_porosity_(old.init_porosity_),
      is_init_porosity_(old.is_init_porosity_)
{
}

template <class So3Ele, Core::FE::CellType distype>
Core::Elements::Element* Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::clone() const
{
  auto* newelement = new Discret::ELEMENTS::So3PoroP1<So3Ele, distype>(*this);
  return newelement;
}

template <class So3Ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::pack(
    Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  data.add_to_pack<int>(is_init_porosity_);

  if (is_init_porosity_) add_to_pack(data, *init_porosity_);

  // add base class Element
  Base::pack(data);
}

template <class So3Ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::unpack(
    Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  is_init_porosity_ = extract_int(buffer);

  if (is_init_porosity_)
  {
    init_porosity_ = Teuchos::make_rcp<Core::LinAlg::Matrix<Base::numnod_, 1>>(true);
    extract_from_pack(buffer, *init_porosity_);
  }


  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer basedata_buffer(basedata);
  Base::unpack(basedata_buffer);

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
}

template <class So3Ele, Core::FE::CellType distype>
std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::surfaces()
{
  return Core::Communication::element_boundary_factory<StructuralSurface, Core::Elements::Element>(
      Core::Communication::buildSurfaces, *this);
}

template <class So3Ele, Core::FE::CellType distype>
std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::lines()
{
  return Core::Communication::element_boundary_factory<StructuralLine, Core::Elements::Element>(
      Core::Communication::buildLines, *this);
}

template <class So3Ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::print(std::ostream& os) const
{
  os << "So3_Poro_P1 ";
  os << Core::FE::cell_type_to_string(distype).c_str() << " ";
  Core::Elements::Element::print(os);
}

template <class So3Ele, Core::FE::CellType distype>
int Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::unique_par_object_id() const
{
  switch (distype)
  {
    case Core::FE::CellType::hex8:
      return SoHex8PoroP1Type::instance().unique_par_object_id();
      break;
    case Core::FE::CellType::tet4:
      return SoTet4PoroP1Type::instance().unique_par_object_id();
      break;
    default:
      FOUR_C_THROW("unknown element type!");
      break;
  }
  return -1;
}

template <class So3Ele, Core::FE::CellType distype>
Core::Elements::ElementType& Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::element_type() const
{
  switch (distype)
  {
    case Core::FE::CellType::tet4:
      return SoTet4PoroP1Type::instance();
    case Core::FE::CellType::hex8:
      return SoHex8PoroP1Type::instance();
    default:
      FOUR_C_THROW("unknown element type!");
      break;
  }
  return SoHex8PoroP1Type::instance();
}

FOUR_C_NAMESPACE_CLOSE

#include "4C_so3_poro_p1_fwd.hpp"
