/*----------------------------------------------------------------------*/
/*! \file
\brief 3d TSI solid element
\level 1
*/


/*----------------------------------------------------------------------*
 | headers                                                   dano 11/12 |
 *----------------------------------------------------------------------*/
#include "baci_so3_thermo.hpp"

#include "baci_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | ctor (public)                                             dano 08/12 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
DRT::ELEMENTS::So3_Thermo<so3_ele, distype>::So3_Thermo(int id, int owner)
    : so3_ele(id, owner), intpoints_(distype)
{
  numgpt_ = intpoints_.NumPoints();
  return;
}


/*----------------------------------------------------------------------*
 | copy-ctor (public)                                        dano 08/12 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
DRT::ELEMENTS::So3_Thermo<so3_ele, distype>::So3_Thermo(
    const DRT::ELEMENTS::So3_Thermo<so3_ele, distype>& old)
    : so3_ele(old), intpoints_(distype)
{
  numgpt_ = intpoints_.NumPoints();
  return;
}


/*----------------------------------------------------------------------*
 | deep copy this instance of Solid3 and return pointer to   dano 08/12 |
 | it (public)                                                          |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
DRT::Element* DRT::ELEMENTS::So3_Thermo<so3_ele, distype>::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::So3_Thermo<so3_ele, distype>(*this);

  return newelement;
}


/*----------------------------------------------------------------------*
 | pack data (public)                                        dano 08/12 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Thermo<so3_ele, distype>::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  so3_ele::AddtoPack(data, type);
  // detJ_
  so3_ele::AddtoPack(data, detJ_);

  // invJ_
  const auto size = (int)invJ_.size();
  so3_ele::AddtoPack(data, size);
  for (int i = 0; i < size; ++i) so3_ele::AddtoPack(data, invJ_[i]);

  // add base class Element
  so3_ele::Pack(data);

  return;

}  // Pack()


/*----------------------------------------------------------------------*
 | unpack data (public)                                      dano 08/12 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Thermo<so3_ele, distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // detJ_
  so3_ele::ExtractfromPack(position, data, detJ_);
  // invJ_
  int size = 0;
  so3_ele::ExtractfromPack(position, data, size);
  invJ_.resize(size, CORE::LINALG::Matrix<nsd_, nsd_>(true));
  for (int i = 0; i < size; ++i) so3_ele::ExtractfromPack(position, data, invJ_[i]);

  // extract base class Element
  std::vector<char> basedata(0);
  so3_ele::ExtractfromPack(position, data, basedata);
  so3_ele::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;

}  // Unpack()


/*----------------------------------------------------------------------*
 | print this element (public)                               dano 08/12 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Thermo<so3_ele, distype>::Print(std::ostream& os) const
{
  os << "So3_Thermo ";
  Element::Print(os);
  return;
}


/*----------------------------------------------------------------------*
 | read this element, get the material (public)              dano 08/12 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
bool DRT::ELEMENTS::So3_Thermo<so3_ele, distype>::ReadElement(
    const std::string& eletype, const std::string& eledistype, INPUT::LineDefinition* linedef)
{
  so3_ele::ReadElement(eletype, eledistype, linedef);

  return true;

}  // ReadElement()

/*----------------------------------------------------------------------*
 | get the nodes from so3 (public)                           dano 05/13 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
int DRT::ELEMENTS::So3_Thermo<so3_ele, distype>::UniqueParObjectId() const
{
  switch (distype)
  {
    case CORE::FE::CellType::hex8:
    {
      // cast the most specialised element
      // otherwise cast fails, because hex8fbar == hex8
      const auto* ele = dynamic_cast<const DRT::ELEMENTS::So_hex8fbar*>(this);
      if (ele != nullptr)
        return So_hex8fbarThermoType::Instance().UniqueParObjectId();
      else
        return So_hex8ThermoType::Instance().UniqueParObjectId();
      break;
    }  // hex8
    case CORE::FE::CellType::tet4:
      return So_tet4ThermoType::Instance().UniqueParObjectId();
      break;
    case CORE::FE::CellType::tet10:
      return So_tet10ThermoType::Instance().UniqueParObjectId();
      break;
    case CORE::FE::CellType::hex27:
      return So_hex27ThermoType::Instance().UniqueParObjectId();
      break;
    case CORE::FE::CellType::hex20:
      return So_hex20ThermoType::Instance().UniqueParObjectId();
      break;
    case CORE::FE::CellType::nurbs27:
      return So_nurbs27ThermoType::Instance().UniqueParObjectId();
      break;
    default:
      dserror("unknown element type!");
      break;
  }
  // Intel compiler needs a return
  return -1;

}  // UniqueParObjectId()


/*----------------------------------------------------------------------*
 | get the nodes from so3 (public)                           dano 05/13 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
DRT::ElementType& DRT::ELEMENTS::So3_Thermo<so3_ele, distype>::ElementType() const
{
  switch (distype)
  {
    case CORE::FE::CellType::hex8:
    {
      // cast the most specialised element
      // caution: otherwise does not work, because hex8fbar == hex8
      const auto* ele = dynamic_cast<const DRT::ELEMENTS::So_hex8fbar*>(this);
      if (ele != nullptr)
        return So_hex8fbarThermoType::Instance();
      else
        return So_hex8ThermoType::Instance();
      break;
    }
    case CORE::FE::CellType::tet4:
      return So_tet4ThermoType::Instance();
      break;
    case CORE::FE::CellType::tet10:
      return So_tet10ThermoType::Instance();
      break;
    case CORE::FE::CellType::hex27:
      return So_hex27ThermoType::Instance();
      break;
    case CORE::FE::CellType::hex20:
      return So_hex20ThermoType::Instance();
      break;
    case CORE::FE::CellType::nurbs27:
      return So_nurbs27ThermoType::Instance();
      break;
    default:
      dserror("unknown element type!");
      break;
  }
  // Intel compiler needs a return
  return So_hex8ThermoType::Instance();

};  // ElementType()


/*----------------------------------------------------------------------*
 | get the nodes from so3 (public)                           dano 08/12 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
inline DRT::Node** DRT::ELEMENTS::So3_Thermo<so3_ele, distype>::Nodes()
{
  return so3_ele::Nodes();
}


/*----------------------------------------------------------------------*
 | get the material from so3 (public)                        dano 08/12 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
inline Teuchos::RCP<MAT::Material> DRT::ELEMENTS::So3_Thermo<so3_ele, distype>::Material() const
{
  return so3_ele::Material();
}


/*----------------------------------------------------------------------*
 | get the node Ids from so3 (public)                        dano 08/12 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
inline int DRT::ELEMENTS::So3_Thermo<so3_ele, distype>::Id() const
{
  return so3_ele::Id();
}


/*----------------------------------------------------------------------*
 | return names of visualization data (public)               dano 04/13 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Thermo<so3_ele, distype>::VisNames(std::map<std::string, int>& names)
{
  so3_ele::VisNames(names);

  return;
}  // VisNames()

/*----------------------------------------------------------------------*
 | return visualization data (public)                        dano 04/13 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
bool DRT::ELEMENTS::So3_Thermo<so3_ele, distype>::VisData(
    const std::string& name, std::vector<double>& data)
{
  return so3_ele::VisData(name, data);

}  // VisData()

FOUR_C_NAMESPACE_CLOSE

/*----------------------------------------------------------------------*/
// include the file at the end of so3_thermo.cpp
#include "baci_so3_thermo_fwd.hpp"
