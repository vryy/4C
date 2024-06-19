/*----------------------------------------------------------------------*/
/*! \file

 \brief implementation of the 3D solid-poro element (p1, mixed approach) including scatra
 functionality

 \level 2

 *----------------------------------------------------------------------*/

#include "4C_so3_poro_p1_scatra.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_so3_poro_p1_scatra_eletypes.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <class so3_ele, Core::FE::CellType distype>
Discret::ELEMENTS::So3PoroP1Scatra<so3_ele, distype>::So3PoroP1Scatra(int id, int owner)
    : So3PoroP1<so3_ele, distype>(id, owner), impltype_(Inpar::ScaTra::impltype_undefined)
{
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                    schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <class so3_ele, Core::FE::CellType distype>
Discret::ELEMENTS::So3PoroP1Scatra<so3_ele, distype>::So3PoroP1Scatra(
    const Discret::ELEMENTS::So3PoroP1Scatra<so3_ele, distype>& old)
    : So3PoroP1<so3_ele, distype>(old), impltype_(old.impltype_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance and return pointer to it (public)           |
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <class so3_ele, Core::FE::CellType distype>
Core::Elements::Element* Discret::ELEMENTS::So3PoroP1Scatra<so3_ele, distype>::Clone() const
{
  auto* newelement = new Discret::ELEMENTS::So3PoroP1Scatra<so3_ele, distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data (public)                                    schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1Scatra<so3_ele, distype>::pack(
    Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  so3_ele::add_to_pack(data, type);

  // pack scalar transport impltype
  so3_ele::add_to_pack(data, impltype_);

  // add base class Element
  my::pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data (public)                                  schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1Scatra<so3_ele, distype>::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract scalar transport impltype_
  impltype_ = static_cast<Inpar::ScaTra::ImplType>(so3_ele::extract_int(position, data));

  // extract base class Element
  std::vector<char> basedata(0);
  my::extract_from_pack(position, data, basedata);
  my::unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                           schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1Scatra<so3_ele, distype>::Print(std::ostream& os) const
{
  os << "So3_Poro_P1_Scatra ";
  os << Core::FE::CellTypeToString(distype).c_str() << " ";
  Core::Elements::Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 | get the unique ParObject Id (public)                    schmidt 09/17|
 *----------------------------------------------------------------------*/
template <class so3_ele, Core::FE::CellType distype>
int Discret::ELEMENTS::So3PoroP1Scatra<so3_ele, distype>::UniqueParObjectId() const
{
  int parobjectid(-1);
  switch (distype)
  {
    case Core::FE::CellType::hex8:
    {
      parobjectid = SoHex8PoroP1ScatraType::Instance().UniqueParObjectId();
      break;
    }
    case Core::FE::CellType::tet4:
    {
      parobjectid = SoTet4PoroP1ScatraType::Instance().UniqueParObjectId();
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
template <class so3_ele, Core::FE::CellType distype>
Core::Elements::ElementType& Discret::ELEMENTS::So3PoroP1Scatra<so3_ele, distype>::ElementType()
    const
{
  switch (distype)
  {
    case Core::FE::CellType::tet4:
      return SoTet4PoroP1ScatraType::Instance();
    case Core::FE::CellType::hex8:
      return SoHex8PoroP1ScatraType::Instance();
    default:
      FOUR_C_THROW("unknown element type!");
      break;
  }
  return SoHex8PoroP1ScatraType::Instance();
}

/*----------------------------------------------------------------------*
 |  read this element (public)                             schmidt 09/17|
 *----------------------------------------------------------------------*/
template <class so3_ele, Core::FE::CellType distype>
bool Discret::ELEMENTS::So3PoroP1Scatra<so3_ele, distype>::ReadElement(
    const std::string& eletype, const std::string& eledistype, Input::LineDefinition* linedef)
{
  // read base element
  my::ReadElement(eletype, eledistype, linedef);

  // read scalar transport implementation type
  std::string impltype;
  linedef->extract_string("TYPE", impltype);

  if (impltype == "Undefined")
    impltype_ = Inpar::ScaTra::impltype_undefined;
  else if (impltype == "AdvReac")
    impltype_ = Inpar::ScaTra::impltype_advreac;
  else if (impltype == "CardMono")
    impltype_ = Inpar::ScaTra::impltype_cardiac_monodomain;
  else if (impltype == "Chemo")
    impltype_ = Inpar::ScaTra::impltype_chemo;
  else if (impltype == "ChemoReac")
    impltype_ = Inpar::ScaTra::impltype_chemoreac;
  else if (impltype == "Loma")
    impltype_ = Inpar::ScaTra::impltype_loma;
  else if (impltype == "Poro")
    impltype_ = Inpar::ScaTra::impltype_poro;
  else if (impltype == "PoroReac")
    impltype_ = Inpar::ScaTra::impltype_pororeac;
  else if (impltype == "PoroReacECM")
    impltype_ = Inpar::ScaTra::impltype_pororeacECM;
  else if (impltype == "PoroMultiReac")
    impltype_ = Inpar::ScaTra::impltype_multipororeac;
  else if (impltype == "RefConcReac")
    impltype_ = Inpar::ScaTra::impltype_refconcreac;
  else if (impltype == "Std")
    impltype_ = Inpar::ScaTra::impltype_std;
  else
    FOUR_C_THROW("Invalid implementation type for So3_Poro_P1_Scatra elements!");

  return true;
}

/*----------------------------------------------------------------------*
 |                                                         schmidt 09/17|
 *----------------------------------------------------------------------*/
template class Discret::ELEMENTS::So3PoroP1Scatra<Discret::ELEMENTS::SoTet4,
    Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::So3PoroP1Scatra<Discret::ELEMENTS::SoHex8,
    Core::FE::CellType::hex8>;

FOUR_C_NAMESPACE_CLOSE
