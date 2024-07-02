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
  numgpt_ = intpoints_.NumPoints();
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
  numgpt_ = intpoints_.NumPoints();
  return;
}


/*----------------------------------------------------------------------*
 | deep copy this instance of Solid3 and return pointer to   dano 08/12 |
 | it (public)                                                          |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
Core::Elements::Element* Discret::ELEMENTS::So3Thermo<So3Ele, distype>::Clone() const
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
  int type = UniqueParObjectId();
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

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

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
bool Discret::ELEMENTS::So3Thermo<So3Ele, distype>::ReadElement(
    const std::string& eletype, const std::string& eledistype, Input::LineDefinition* linedef)
{
  So3Ele::ReadElement(eletype, eledistype, linedef);

  return true;

}  // ReadElement()

/*----------------------------------------------------------------------*
 | get the nodes from so3 (public)                           dano 05/13 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
int Discret::ELEMENTS::So3Thermo<So3Ele, distype>::UniqueParObjectId() const
{
  switch (distype)
  {
    case Core::FE::CellType::hex8:
    {
      // cast the most specialised element
      // otherwise cast fails, because hex8fbar == hex8
      const auto* ele = dynamic_cast<const Discret::ELEMENTS::SoHex8fbar*>(this);
      if (ele != nullptr)
        return SoHex8fbarThermoType::Instance().UniqueParObjectId();
      else
        return SoHex8ThermoType::Instance().UniqueParObjectId();
      break;
    }  // hex8
    case Core::FE::CellType::tet4:
      return SoTet4ThermoType::Instance().UniqueParObjectId();
      break;
    case Core::FE::CellType::tet10:
      return SoTet10ThermoType::Instance().UniqueParObjectId();
      break;
    case Core::FE::CellType::hex27:
      return SoHex27ThermoType::Instance().UniqueParObjectId();
      break;
    case Core::FE::CellType::hex20:
      return SoHex20ThermoType::Instance().UniqueParObjectId();
      break;
    case Core::FE::CellType::nurbs27:
      return SoNurbs27ThermoType::Instance().UniqueParObjectId();
      break;
    default:
      FOUR_C_THROW("unknown element type!");
      break;
  }
  // Intel compiler needs a return
  return -1;

}  // UniqueParObjectId()


/*----------------------------------------------------------------------*
 | get the nodes from so3 (public)                           dano 05/13 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
Core::Elements::ElementType& Discret::ELEMENTS::So3Thermo<So3Ele, distype>::ElementType() const
{
  switch (distype)
  {
    case Core::FE::CellType::hex8:
    {
      // cast the most specialised element
      // caution: otherwise does not work, because hex8fbar == hex8
      const auto* ele = dynamic_cast<const Discret::ELEMENTS::SoHex8fbar*>(this);
      if (ele != nullptr)
        return SoHex8fbarThermoType::Instance();
      else
        return SoHex8ThermoType::Instance();
      break;
    }
    case Core::FE::CellType::tet4:
      return SoTet4ThermoType::Instance();
      break;
    case Core::FE::CellType::tet10:
      return SoTet10ThermoType::Instance();
      break;
    case Core::FE::CellType::hex27:
      return SoHex27ThermoType::Instance();
      break;
    case Core::FE::CellType::hex20:
      return SoHex20ThermoType::Instance();
      break;
    case Core::FE::CellType::nurbs27:
      return SoNurbs27ThermoType::Instance();
      break;
    default:
      FOUR_C_THROW("unknown element type!");
      break;
  }
  // Intel compiler needs a return
  return SoHex8ThermoType::Instance();

};  // ElementType()


/*----------------------------------------------------------------------*
 | get the nodes from so3 (public)                           dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
inline Core::Nodes::Node** Discret::ELEMENTS::So3Thermo<So3Ele, distype>::Nodes()
{
  return So3Ele::Nodes();
}


/*----------------------------------------------------------------------*
 | get the material from so3 (public)                        dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
inline Teuchos::RCP<Core::Mat::Material> Discret::ELEMENTS::So3Thermo<So3Ele, distype>::material()
    const
{
  return So3Ele::Material();
}


/*----------------------------------------------------------------------*
 | get the node Ids from so3 (public)                        dano 08/12 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
inline int Discret::ELEMENTS::So3Thermo<So3Ele, distype>::id() const
{
  return So3Ele::Id();
}


/*----------------------------------------------------------------------*
 | return names of visualization data (public)               dano 04/13 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3Thermo<So3Ele, distype>::VisNames(std::map<std::string, int>& names)
{
  So3Ele::VisNames(names);

  return;
}  // VisNames()

/*----------------------------------------------------------------------*
 | return visualization data (public)                        dano 04/13 |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
bool Discret::ELEMENTS::So3Thermo<So3Ele, distype>::VisData(
    const std::string& name, std::vector<double>& data)
{
  return So3Ele::VisData(name, data);

}  // VisData()

FOUR_C_NAMESPACE_CLOSE

/*----------------------------------------------------------------------*/
// include the file at the end of so3_thermo.cpp
#include "4C_so3_thermo_fwd.hpp"
