/*----------------------------------------------------------------------*/
/*! \file

\level 3


\brief Nonlinear Membrane Finite Element with ScaTra coupling

*----------------------------------------------------------------------*/

#include "4C_membrane_scatra.hpp"

#include "4C_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  constructor (public)                                   sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::MembraneScatra<distype>::MembraneScatra(int id, int owner)
    : Membrane<distype>(id, owner), impltype_(Inpar::ScaTra::impltype_undefined)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-constructor (public)                              sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::MembraneScatra<distype>::MembraneScatra(
    const Discret::ELEMENTS::MembraneScatra<distype>& old)
    : Membrane<distype>(old), impltype_(old.impltype_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of MembraneScatra              sfuchs 05/18 |
 |  and return pointer to it (public)                                   |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Core::Elements::Element* Discret::ELEMENTS::MembraneScatra<distype>::Clone() const
{
  Discret::ELEMENTS::MembraneScatra<distype>* newelement =
      new Discret::ELEMENTS::MembraneScatra<distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::MembraneScatra<distype>::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  Core::Communication::ParObject::add_to_pack(data, type);

  // pack scalar transport impltype_
  Core::Communication::ParObject::add_to_pack(data, impltype_);

  // add base class Element
  Membrane<distype>::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::MembraneScatra<distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract scalar transport impltype
  impltype_ = static_cast<Inpar::ScaTra::ImplType>(
      Core::Communication::ParObject::ExtractInt(position, data));

  // extract base class Element
  std::vector<char> basedata(0);
  Core::Communication::ParObject::extract_from_pack(position, data, basedata);
  Membrane<distype>::Unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                            sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::MembraneScatra<distype>::Print(std::ostream& os) const
{
  os << "MembraneScatra ";
  Membrane<distype>::Print(os);

  return;
}

/*----------------------------------------------------------------------*
 |  read this element (public)                             sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
bool Discret::ELEMENTS::MembraneScatra<distype>::ReadElement(
    const std::string& eletype, const std::string& eledistype, Input::LineDefinition* linedef)
{
  // read base element
  Membrane<distype>::ReadElement(eletype, eledistype, linedef);

  // read scalar transport implementation type
  std::string impltype;
  linedef->ExtractString("TYPE", impltype);

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
  else if (impltype == "RefConcReac")
    impltype_ = Inpar::ScaTra::impltype_refconcreac;
  else if (impltype == "Std")
    impltype_ = Inpar::ScaTra::impltype_std;
  else
    FOUR_C_THROW("Invalid implementation type for Wall1_Scatra elements!");

  return true;
}

/*----------------------------------------------------------------------*
 |  Get vector of ptrs to nodes (private)                  sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
inline Core::Nodes::Node** Discret::ELEMENTS::MembraneScatra<distype>::Nodes()
{
  return Membrane<distype>::Nodes();
}

/*----------------------------------------------------------------------*
 |  Get shape type of element (private)                    sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Core::FE::CellType Discret::ELEMENTS::MembraneScatra<distype>::Shape() const
{
  return Membrane<distype>::Shape();
}

template class Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::quad4>;
template class Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::quad9>;

FOUR_C_NAMESPACE_CLOSE
