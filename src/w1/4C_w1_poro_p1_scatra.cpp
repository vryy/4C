/*----------------------------------------------------------------------------*/
/*! \file
\brief A 2D wall element for solid-part of porous medium using p1 (mixed) approach including scatra
 functionality.

\level 2


*/
/*---------------------------------------------------------------------------*/

#include "4C_w1_poro_p1_scatra.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_w1_poro_p1_scatra_eletypes.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::Wall1PoroP1Scatra<distype>::Wall1PoroP1Scatra(int id, int owner)
    : Discret::ELEMENTS::Wall1PoroP1<distype>(id, owner),
      impltype_(Inpar::ScaTra::impltype_undefined)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                    schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::Wall1PoroP1Scatra<distype>::Wall1PoroP1Scatra(
    const Discret::ELEMENTS::Wall1PoroP1Scatra<distype>& old)
    : Discret::ELEMENTS::Wall1PoroP1<distype>(old), impltype_(old.impltype_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance and return pointer to it (public)           |
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Core::Elements::Element* Discret::ELEMENTS::Wall1PoroP1Scatra<distype>::Clone() const
{
  Discret::ELEMENTS::Wall1PoroP1Scatra<distype>* newelement =
      new Discret::ELEMENTS::Wall1PoroP1Scatra<distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data (public)                                    schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::Wall1PoroP1Scatra<distype>::Pack(
    Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  my::add_to_pack(data, type);
  // pack scalar transport impltype
  my::add_to_pack(data, impltype_);

  // add base class Element
  my::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data (public)                                  schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::Wall1PoroP1Scatra<distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract scalar transport impltype
  impltype_ = static_cast<Inpar::ScaTra::ImplType>(my::extract_int(position, data));

  // extract base class Element
  std::vector<char> basedata(0);
  my::extract_from_pack(position, data, basedata);
  my::Unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}

/*----------------------------------------------------------------------*
 |  read this element (public)                             schmidt 09/17|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
bool Discret::ELEMENTS::Wall1PoroP1Scatra<distype>::ReadElement(
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
    FOUR_C_THROW("Invalid implementation type for Wall1_PoroP1Scatra elements!");

  return true;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                           schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::Wall1PoroP1Scatra<distype>::Print(std::ostream& os) const
{
  os << "Wall1_PoroP1Scatra ";
  Core::Elements::Element::Print(os);
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                           schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::Wall1PoroP1Scatra<distype>::UniqueParObjectId() const
{
  int parobjectid(-1);
  switch (distype)
  {
    case Core::FE::CellType::tri3:
    {
      parobjectid = Discret::ELEMENTS::WallTri3PoroP1ScatraType::Instance().UniqueParObjectId();
      break;
    }
    case Core::FE::CellType::quad4:
    {
      parobjectid = Discret::ELEMENTS::WallQuad4PoroP1ScatraType::Instance().UniqueParObjectId();
      break;
    }
    case Core::FE::CellType::quad9:
    {
      parobjectid = Discret::ELEMENTS::WallQuad9PoroP1ScatraType::Instance().UniqueParObjectId();
      break;
    }
    default:
    {
      FOUR_C_THROW("unknown element type");
      break;
    }
  }
  return parobjectid;
}

/*----------------------------------------------------------------------*
 | get the element type (public)                           schmidt 09/17|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Core::Elements::ElementType& Discret::ELEMENTS::Wall1PoroP1Scatra<distype>::ElementType() const
{
  switch (distype)
  {
    case Core::FE::CellType::tri3:
      return Discret::ELEMENTS::WallTri3PoroP1ScatraType::Instance();
      break;
    case Core::FE::CellType::quad4:
      return Discret::ELEMENTS::WallQuad4PoroP1ScatraType::Instance();
      break;
    case Core::FE::CellType::quad9:
      return Discret::ELEMENTS::WallQuad9PoroP1ScatraType::Instance();
      break;
    default:
      FOUR_C_THROW("unknown element type");
      break;
  }
  return Discret::ELEMENTS::WallQuad4PoroP1ScatraType::Instance();
}

/*----------------------------------------------------------------------*
 *                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
template class Discret::ELEMENTS::Wall1PoroP1Scatra<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::Wall1PoroP1Scatra<Core::FE::CellType::quad4>;
template class Discret::ELEMENTS::Wall1PoroP1Scatra<Core::FE::CellType::quad9>;

FOUR_C_NAMESPACE_CLOSE
