/*----------------------------------------------------------------------*/
/*! \file
\brief 3d TSI solid element
\level 1
*/


/*----------------------------------------------------------------------*
 | headers                                                   dano 11/12 |
 *----------------------------------------------------------------------*/
#include "4C_so3_thermo.hpp"

#include "4C_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | ctor (public)                                             dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
Discret::ELEMENTS::So3Thermo<So3Ele, distype>::So3Thermo(int id, int owner)
    : So3Ele(id, owner), intpoints_(distype)
{
  numgpt_ = intpoints_.num_points();
  return;
}


/*----------------------------------------------------------------------*
 | copy-ctor (public)                                        dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
Discret::ELEMENTS::So3Thermo<So3Ele, distype>::So3Thermo(
    const Discret::ELEMENTS::So3Thermo<So3Ele, distype>& old)
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
Core::Elements::Element* Discret::ELEMENTS::So3Thermo<So3Ele, distype>::clone() const
{
  auto* newelement = new Discret::ELEMENTS::So3Thermo<So3Ele, distype>(*this);

  return newelement;
}


/*----------------------------------------------------------------------*
 | pack data (public)                                        dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3Thermo<So3Ele, distype>::pack(
    Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  So3Ele::add_to_pack(data, type);
  // detJ_
  So3Ele::add_to_pack(data, detJ_);

  // invJ_
  const auto size = (int)invJ_.size();
  So3Ele::add_to_pack(data, size);
  for (int i = 0; i < size; ++i) So3Ele::add_to_pack(data, invJ_[i]);

  // add base class Element
  So3Ele::pack(data);

  return;

}  // pack()


/*----------------------------------------------------------------------*
 | unpack data (public)                                      dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3Thermo<So3Ele, distype>::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // detJ_
  So3Ele::extract_from_pack(position, data, detJ_);
  // invJ_
  int size = 0;
  So3Ele::extract_from_pack(position, data, size);
  invJ_.resize(size, Core::LinAlg::Matrix<nsd_, nsd_>(true));
  for (int i = 0; i < size; ++i) So3Ele::extract_from_pack(position, data, invJ_[i]);

  // extract base class Element
  std::vector<char> basedata(0);
  So3Ele::extract_from_pack(position, data, basedata);
  So3Ele::unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;

}  // unpack()


/*----------------------------------------------------------------------*
 | print this element (public)                               dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3Thermo<So3Ele, distype>::print(std::ostream& os) const
{
  os << "So3_Thermo ";
  Core::Elements::Element::print(os);
  return;
}


/*----------------------------------------------------------------------*
 | read this element, get the material (public)              dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
bool Discret::ELEMENTS::So3Thermo<So3Ele, distype>::read_element(const std::string& eletype,
    const std::string& eledistype, const Core::IO::InputParameterContainer& container)
{
  So3Ele::read_element(eletype, eledistype, container);

  return true;

}  // read_element()

/*----------------------------------------------------------------------*
 | get the nodes from so3 (public)                           dano 05/13 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
int Discret::ELEMENTS::So3Thermo<So3Ele, distype>::unique_par_object_id() const
{
  switch (distype)
  {
    case Core::FE::CellType::hex8:
    {
      // cast the most specialised element
      // otherwise cast fails, because hex8fbar == hex8
      const auto* ele = dynamic_cast<const Discret::ELEMENTS::SoHex8fbar*>(this);
      if (ele != nullptr)
        return SoHex8fbarThermoType::instance().unique_par_object_id();
      else
        return SoHex8ThermoType::instance().unique_par_object_id();
      break;
    }  // hex8
    case Core::FE::CellType::tet4:
      return SoTet4ThermoType::instance().unique_par_object_id();
      break;
    case Core::FE::CellType::tet10:
      return SoTet10ThermoType::instance().unique_par_object_id();
      break;
    case Core::FE::CellType::hex27:
      return SoHex27ThermoType::instance().unique_par_object_id();
      break;
    case Core::FE::CellType::hex20:
      return SoHex20ThermoType::instance().unique_par_object_id();
      break;
    case Core::FE::CellType::nurbs27:
      return SoNurbs27ThermoType::instance().unique_par_object_id();
      break;
    default:
      FOUR_C_THROW("unknown element type!");
      break;
  }
  // Intel compiler needs a return
  return -1;

}  // unique_par_object_id()


/*----------------------------------------------------------------------*
 | get the nodes from so3 (public)                           dano 05/13 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
Core::Elements::ElementType& Discret::ELEMENTS::So3Thermo<So3Ele, distype>::element_type() const
{
  switch (distype)
  {
    case Core::FE::CellType::hex8:
    {
      // cast the most specialised element
      // caution: otherwise does not work, because hex8fbar == hex8
      const auto* ele = dynamic_cast<const Discret::ELEMENTS::SoHex8fbar*>(this);
      if (ele != nullptr)
        return SoHex8fbarThermoType::instance();
      else
        return SoHex8ThermoType::instance();
      break;
    }
    case Core::FE::CellType::tet4:
      return SoTet4ThermoType::instance();
      break;
    case Core::FE::CellType::tet10:
      return SoTet10ThermoType::instance();
      break;
    case Core::FE::CellType::hex27:
      return SoHex27ThermoType::instance();
      break;
    case Core::FE::CellType::hex20:
      return SoHex20ThermoType::instance();
      break;
    case Core::FE::CellType::nurbs27:
      return SoNurbs27ThermoType::instance();
      break;
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
inline Core::Nodes::Node** Discret::ELEMENTS::So3Thermo<So3Ele, distype>::nodes()
{
  return So3Ele::nodes();
}


/*----------------------------------------------------------------------*
 | get the material from so3 (public)                        dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
inline Teuchos::RCP<Core::Mat::Material> Discret::ELEMENTS::So3Thermo<So3Ele, distype>::material()
    const
{
  return So3Ele::material();
}


/*----------------------------------------------------------------------*
 | get the node Ids from so3 (public)                        dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
inline int Discret::ELEMENTS::So3Thermo<So3Ele, distype>::id() const
{
  return So3Ele::id();
}


/*----------------------------------------------------------------------*
 | return names of visualization data (public)               dano 04/13 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3Thermo<So3Ele, distype>::vis_names(std::map<std::string, int>& names)
{
  So3Ele::vis_names(names);

  return;
}  // vis_names()

/*----------------------------------------------------------------------*
 | return visualization data (public)                        dano 04/13 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
bool Discret::ELEMENTS::So3Thermo<So3Ele, distype>::vis_data(
    const std::string& name, std::vector<double>& data)
{
  return So3Ele::vis_data(name, data);

}  // vis_data()

FOUR_C_NAMESPACE_CLOSE

/*----------------------------------------------------------------------*/
// include the file at the end of so3_thermo.cpp
#include "4C_so3_thermo_fwd.hpp"
