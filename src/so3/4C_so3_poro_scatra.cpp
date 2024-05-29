/*----------------------------------------------------------------------*/
/*! \file

 \brief implementation of the 3D solid-poro element including scatra functionality

 \level 2

*----------------------------------------------------------------------*/

#include "4C_so3_poro_scatra.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_so3_poro_scatra_eletypes.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
DRT::ELEMENTS::So3PoroScatra<so3_ele, distype>::So3PoroScatra(int id, int owner)
    : So3Poro<so3_ele, distype>(id, owner), impltype_(INPAR::SCATRA::impltype_undefined)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                    schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
DRT::ELEMENTS::So3PoroScatra<so3_ele, distype>::So3PoroScatra(
    const DRT::ELEMENTS::So3PoroScatra<so3_ele, distype>& old)
    : So3Poro<so3_ele, distype>(old), impltype_(old.impltype_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance and return pointer to it (public)           |
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
CORE::Elements::Element* DRT::ELEMENTS::So3PoroScatra<so3_ele, distype>::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::So3PoroScatra<so3_ele, distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data (public)                                    schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3PoroScatra<so3_ele, distype>::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  so3_ele::AddtoPack(data, type);
  // pack scalar transport impltype
  so3_ele::AddtoPack(data, impltype_);

  // add base class Element
  my::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data (public)                                  schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3PoroScatra<so3_ele, distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract scalar transport impltype_
  impltype_ = static_cast<INPAR::SCATRA::ImplType>(so3_ele::ExtractInt(position, data));

  // extract base class Element
  std::vector<char> basedata(0);
  my::ExtractfromPack(position, data, basedata);
  my::Unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                           schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3PoroScatra<so3_ele, distype>::Print(std::ostream& os) const
{
  os << "So3_Poro_Scatra ";
  os << CORE::FE::CellTypeToString(distype).c_str() << " ";
  CORE::Elements::Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  read this element (public)                             schmidt 09/17|
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
bool DRT::ELEMENTS::So3PoroScatra<so3_ele, distype>::ReadElement(
    const std::string& eletype, const std::string& eledistype, INPUT::LineDefinition* linedef)
{
  // read base element
  my::ReadElement(eletype, eledistype, linedef);

  // read scalar transport implementation type
  std::string impltype;
  linedef->ExtractString("TYPE", impltype);

  if (impltype == "Undefined")
    impltype_ = INPAR::SCATRA::impltype_undefined;
  else if (impltype == "AdvReac")
    impltype_ = INPAR::SCATRA::impltype_advreac;
  else if (impltype == "CardMono")
    impltype_ = INPAR::SCATRA::impltype_cardiac_monodomain;
  else if (impltype == "Chemo")
    impltype_ = INPAR::SCATRA::impltype_chemo;
  else if (impltype == "ChemoReac")
    impltype_ = INPAR::SCATRA::impltype_chemoreac;
  else if (impltype == "Loma")
    impltype_ = INPAR::SCATRA::impltype_loma;
  else if (impltype == "Poro")
    impltype_ = INPAR::SCATRA::impltype_poro;
  else if (impltype == "PoroReac")
    impltype_ = INPAR::SCATRA::impltype_pororeac;
  else if (impltype == "PoroReacECM")
    impltype_ = INPAR::SCATRA::impltype_pororeacECM;
  else if (impltype == "PoroMultiReac")
    impltype_ = INPAR::SCATRA::impltype_multipororeac;
  else if (impltype == "RefConcReac")
    impltype_ = INPAR::SCATRA::impltype_refconcreac;
  else if (impltype == "Std")
    impltype_ = INPAR::SCATRA::impltype_std;
  else
    FOUR_C_THROW("Invalid implementation type for So3_Poro_Scatra elements!");

  return true;
}

/*----------------------------------------------------------------------*
 | get the unique ParObject Id (public)                    schmidt 09/17|
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
int DRT::ELEMENTS::So3PoroScatra<so3_ele, distype>::UniqueParObjectId() const
{
  int parobjectid(-1);
  switch (distype)
  {
    case CORE::FE::CellType::tet4:
    {
      parobjectid = SoTet4PoroScatraType::Instance().UniqueParObjectId();
      break;
    }
    case CORE::FE::CellType::tet10:
    {
      parobjectid = SoTet10PoroScatraType::Instance().UniqueParObjectId();
      break;
    }
    case CORE::FE::CellType::hex8:
    {
      parobjectid = SoHex8PoroScatraType::Instance().UniqueParObjectId();
      break;
    }
    case CORE::FE::CellType::hex27:
    {
      parobjectid = SoHex27PoroScatraType::Instance().UniqueParObjectId();
      break;
    }
    case CORE::FE::CellType::nurbs27:
    {
      parobjectid = SoNurbs27PoroScatraType::Instance().UniqueParObjectId();
      break;
    }
    default:
    {
      FOUR_C_THROW("unknown element type!");
      break;
    }
  }
  return parobjectid;
}

/*----------------------------------------------------------------------*
 | get the element type (public)                           schmidt 09/17|
 *----------------------------------------------------------------------*/
template <class so3_ele, CORE::FE::CellType distype>
CORE::Elements::ElementType& DRT::ELEMENTS::So3PoroScatra<so3_ele, distype>::ElementType() const
{
  switch (distype)
  {
    case CORE::FE::CellType::tet4:
      return SoTet4PoroScatraType::Instance();
    case CORE::FE::CellType::tet10:
      return SoTet10PoroScatraType::Instance();
    case CORE::FE::CellType::hex8:
      return SoHex8PoroScatraType::Instance();
    case CORE::FE::CellType::hex27:
      return SoHex27PoroScatraType::Instance();
    case CORE::FE::CellType::nurbs27:
      return SoNurbs27PoroScatraType::Instance();
    default:
      FOUR_C_THROW("unknown element type!");
      break;
  }
  return SoHex8PoroScatraType::Instance();
};

/*----------------------------------------------------------------------*
 |                                                         schmidt 09/17|
 *----------------------------------------------------------------------*/
template class DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoTet4, CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoTet10, CORE::FE::CellType::tet10>;
template class DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoHex8, CORE::FE::CellType::hex8>;
template class DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::SoHex27, CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::So3PoroScatra<DRT::ELEMENTS::NURBS::SoNurbs27,
    CORE::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
