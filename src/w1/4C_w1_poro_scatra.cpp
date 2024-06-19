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
template <Core::FE::CellType distype>
Discret::ELEMENTS::Wall1PoroScatra<distype>::Wall1PoroScatra(int id, int owner)
    : Discret::ELEMENTS::Wall1Poro<distype>(id, owner), impltype_(Inpar::ScaTra::impltype_undefined)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                    schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::Wall1PoroScatra<distype>::Wall1PoroScatra(
    const Discret::ELEMENTS::Wall1PoroScatra<distype>& old)
    : Discret::ELEMENTS::Wall1Poro<distype>(old), impltype_(old.impltype_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance and return pointer to it (public)           |
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Core::Elements::Element* Discret::ELEMENTS::Wall1PoroScatra<distype>::Clone() const
{
  Discret::ELEMENTS::Wall1PoroScatra<distype>* newelement =
      new Discret::ELEMENTS::Wall1PoroScatra<distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data (public)                                    schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::Wall1PoroScatra<distype>::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  my::add_to_pack(data, type);
  // pack scalar transport impltype
  my::add_to_pack(data, impltype_);

  // add base class Element
  my::pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data (public)                                  schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::Wall1PoroScatra<distype>::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract scalar transport impltype_
  impltype_ = static_cast<Inpar::ScaTra::ImplType>(my::extract_int(position, data));

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
template <Core::FE::CellType distype>
void Discret::ELEMENTS::Wall1PoroScatra<distype>::Print(std::ostream& os) const
{
  os << "Wall1_Poro_Scatra ";
  Core::Elements::Element::Print(os);
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 |  read this element (public)                             schmidt 09/17|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
bool Discret::ELEMENTS::Wall1PoroScatra<distype>::ReadElement(
    const std::string& eletype, const std::string& eledistype, Input::LineDefinition* linedef)
{
  // read base element
  my::ReadElement(eletype, eledistype, linedef);

  // read implementation type
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
    FOUR_C_THROW("Invalid implementation type for Wall1_Poro_Scatra elements!");

  return true;
}


/*----------------------------------------------------------------------*
 |                                                         schmidt 09/17|
 *----------------------------------------------------------------------*/
template class Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::quad4>;
template class Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::quad9>;
template class Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::nurbs4>;
template class Discret::ELEMENTS::Wall1PoroScatra<Core::FE::CellType::nurbs9>;

FOUR_C_NAMESPACE_CLOSE
