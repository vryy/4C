/*----------------------------------------------------------------------------*/
/*! \file
\brief 2D wall element for structure part of porous medium including scatra functionality.

\level 2


*/
/*---------------------------------------------------------------------------*/

#include "4C_w1_poro_scatra.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_w1_poro_scatra_eletypes.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::Wall1PoroScatra<distype>::Wall1PoroScatra(int id, int owner)
    : DRT::ELEMENTS::Wall1Poro<distype>(id, owner), impltype_(INPAR::SCATRA::impltype_undefined)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                    schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::Wall1PoroScatra<distype>::Wall1PoroScatra(
    const DRT::ELEMENTS::Wall1PoroScatra<distype>& old)
    : DRT::ELEMENTS::Wall1Poro<distype>(old), impltype_(old.impltype_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance and return pointer to it (public)           |
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::Element* DRT::ELEMENTS::Wall1PoroScatra<distype>::Clone() const
{
  DRT::ELEMENTS::Wall1PoroScatra<distype>* newelement =
      new DRT::ELEMENTS::Wall1PoroScatra<distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data (public)                                    schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Wall1PoroScatra<distype>::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  my::AddtoPack(data, type);
  // pack scalar transport impltype
  my::AddtoPack(data, impltype_);

  // add base class Element
  my::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data (public)                                  schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Wall1PoroScatra<distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract scalar transport impltype_
  impltype_ = static_cast<INPAR::SCATRA::ImplType>(my::ExtractInt(position, data));

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
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Wall1PoroScatra<distype>::Print(std::ostream& os) const
{
  os << "Wall1_Poro_Scatra ";
  Element::Print(os);
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 |  read this element (public)                             schmidt 09/17|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
bool DRT::ELEMENTS::Wall1PoroScatra<distype>::ReadElement(
    const std::string& eletype, const std::string& eledistype, INPUT::LineDefinition* linedef)
{
  // read base element
  my::ReadElement(eletype, eledistype, linedef);

  // read implementation type
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
    FOUR_C_THROW("Invalid implementation type for Wall1_Poro_Scatra elements!");

  return true;
}


/*----------------------------------------------------------------------*
 |                                                         schmidt 09/17|
 *----------------------------------------------------------------------*/
template class DRT::ELEMENTS::Wall1PoroScatra<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::Wall1PoroScatra<CORE::FE::CellType::quad4>;
template class DRT::ELEMENTS::Wall1PoroScatra<CORE::FE::CellType::quad9>;
template class DRT::ELEMENTS::Wall1PoroScatra<CORE::FE::CellType::nurbs4>;
template class DRT::ELEMENTS::Wall1PoroScatra<CORE::FE::CellType::nurbs9>;

FOUR_C_NAMESPACE_CLOSE
