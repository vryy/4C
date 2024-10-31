// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_so3_thermo.hpp"

#include "4C_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | ctor (public)                                             dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
Discret::Elements::So3Thermo<So3Ele, distype>::So3Thermo(int id, int owner)
    : So3Ele(id, owner), intpoints_(distype)
{
  numgpt_ = intpoints_.num_points();
  return;
}


/*----------------------------------------------------------------------*
 | copy-ctor (public)                                        dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
Discret::Elements::So3Thermo<So3Ele, distype>::So3Thermo(
    const Discret::Elements::So3Thermo<So3Ele, distype>& old)
    : So3Ele(old), intpoints_(distype)
{
  numgpt_ = intpoints_.num_points();
  return;
}


/*----------------------------------------------------------------------*
 | deep copy this instance of Solid3 and return pointer to   dano 08/12 |
 | it (public)                                                          |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
Core::Elements::Element* Discret::Elements::So3Thermo<So3Ele, distype>::clone() const
{
  auto* newelement = new Discret::Elements::So3Thermo<So3Ele, distype>(*this);

  return newelement;
}


/*----------------------------------------------------------------------*
 | pack data (public)                                        dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::Elements::So3Thermo<So3Ele, distype>::pack(
    Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // detJ_
  add_to_pack(data, detJ_);

  // invJ_
  const auto size = (int)invJ_.size();
  add_to_pack(data, size);
  for (int i = 0; i < size; ++i) add_to_pack(data, invJ_[i]);

  // add base class Element
  So3Ele::pack(data);

  return;

}  // pack()


/*----------------------------------------------------------------------*
 | unpack data (public)                                      dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::Elements::So3Thermo<So3Ele, distype>::unpack(
    Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // detJ_
  extract_from_pack(buffer, detJ_);
  // invJ_
  int size = 0;
  extract_from_pack(buffer, size);
  invJ_.resize(size, Core::LinAlg::Matrix<nsd_, nsd_>(true));
  for (int i = 0; i < size; ++i) extract_from_pack(buffer, invJ_[i]);

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer basedata_buffer(basedata);
  So3Ele::unpack(basedata_buffer);

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
  return;

}  // unpack()


/*----------------------------------------------------------------------*
 | print this element (public)                               dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::Elements::So3Thermo<So3Ele, distype>::print(std::ostream& os) const
{
  os << "So3_Thermo ";
  Core::Elements::Element::print(os);
  return;
}


/*----------------------------------------------------------------------*
 | read this element, get the material (public)              dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
bool Discret::Elements::So3Thermo<So3Ele, distype>::read_element(const std::string& eletype,
    const std::string& eledistype, const Core::IO::InputParameterContainer& container)
{
  So3Ele::read_element(eletype, eledistype, container);

  return true;

}  // read_element()

/*----------------------------------------------------------------------*
 | get the nodes from so3 (public)                           dano 05/13 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
int Discret::Elements::So3Thermo<So3Ele, distype>::unique_par_object_id() const
{
  switch (distype)
  {
    case Core::FE::CellType::hex8:
      return SoHex8ThermoType::instance().unique_par_object_id();
    default:
      FOUR_C_THROW("unknown element type!");
      break;
  }
}  // unique_par_object_id()


/*----------------------------------------------------------------------*
 | get the nodes from so3 (public)                           dano 05/13 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
Core::Elements::ElementType& Discret::Elements::So3Thermo<So3Ele, distype>::element_type() const
{
  switch (distype)
  {
    case Core::FE::CellType::hex8:
      return SoHex8ThermoType::instance();
    default:
      FOUR_C_THROW("unknown element type!");
      break;
  }
  // Intel compiler needs a return
  return SoHex8ThermoType::instance();

};  // element_type()


/*----------------------------------------------------------------------*
 | get the nodes from so3 (public)                           dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
inline Core::Nodes::Node** Discret::Elements::So3Thermo<So3Ele, distype>::nodes()
{
  return So3Ele::nodes();
}


/*----------------------------------------------------------------------*
 | get the material from so3 (public)                        dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
inline Teuchos::RCP<Core::Mat::Material> Discret::Elements::So3Thermo<So3Ele, distype>::material()
    const
{
  return So3Ele::material();
}


/*----------------------------------------------------------------------*
 | get the node Ids from so3 (public)                        dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
inline int Discret::Elements::So3Thermo<So3Ele, distype>::id() const
{
  return So3Ele::id();
}


/*----------------------------------------------------------------------*
 | return names of visualization data (public)               dano 04/13 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::Elements::So3Thermo<So3Ele, distype>::vis_names(std::map<std::string, int>& names)
{
  So3Ele::vis_names(names);

  return;
}  // vis_names()

/*----------------------------------------------------------------------*
 | return visualization data (public)                        dano 04/13 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
bool Discret::Elements::So3Thermo<So3Ele, distype>::vis_data(
    const std::string& name, std::vector<double>& data)
{
  return So3Ele::vis_data(name, data);

}  // vis_data()

FOUR_C_NAMESPACE_CLOSE

/*----------------------------------------------------------------------*/
// include the file at the end of so3_thermo.cpp
#include "4C_so3_thermo_fwd.hpp"
