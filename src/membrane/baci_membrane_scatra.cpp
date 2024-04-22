/*----------------------------------------------------------------------*/
/*! \file

\level 3


\brief Nonlinear Membrane Finite Element with ScaTra coupling

*----------------------------------------------------------------------*/

#include "baci_membrane_scatra.hpp"

#include "baci_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  constructor (public)                                   sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::MembraneScatra<distype>::MembraneScatra(int id, int owner)
    : Membrane<distype>(id, owner), impltype_(INPAR::SCATRA::impltype_undefined)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-constructor (public)                              sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::MembraneScatra<distype>::MembraneScatra(
    const DRT::ELEMENTS::MembraneScatra<distype>& old)
    : Membrane<distype>(old), impltype_(old.impltype_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of MembraneScatra              sfuchs 05/18 |
 |  and return pointer to it (public)                                   |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::Element* DRT::ELEMENTS::MembraneScatra<distype>::Clone() const
{
  DRT::ELEMENTS::MembraneScatra<distype>* newelement =
      new DRT::ELEMENTS::MembraneScatra<distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::MembraneScatra<distype>::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  CORE::COMM::ParObject::AddtoPack(data, type);

  // pack scalar transport impltype_
  CORE::COMM::ParObject::AddtoPack(data, impltype_);

  // add base class Element
  Membrane<distype>::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::MembraneScatra<distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract scalar transport impltype
  impltype_ =
      static_cast<INPAR::SCATRA::ImplType>(CORE::COMM::ParObject::ExtractInt(position, data));

  // extract base class Element
  std::vector<char> basedata(0);
  CORE::COMM::ParObject::ExtractfromPack(position, data, basedata);
  Membrane<distype>::Unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                            sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::MembraneScatra<distype>::Print(std::ostream& os) const
{
  os << "MembraneScatra ";
  Membrane<distype>::Print(os);

  return;
}

/*----------------------------------------------------------------------*
 |  read this element (public)                             sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
bool DRT::ELEMENTS::MembraneScatra<distype>::ReadElement(
    const std::string& eletype, const std::string& eledistype, INPUT::LineDefinition* linedef)
{
  // read base element
  Membrane<distype>::ReadElement(eletype, eledistype, linedef);

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
  else if (impltype == "RefConcReac")
    impltype_ = INPAR::SCATRA::impltype_refconcreac;
  else if (impltype == "Std")
    impltype_ = INPAR::SCATRA::impltype_std;
  else
    FOUR_C_THROW("Invalid implementation type for Wall1_Scatra elements!");

  return true;
}

/*----------------------------------------------------------------------*
 |  Get vector of ptrs to nodes (private)                  sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
inline DRT::Node** DRT::ELEMENTS::MembraneScatra<distype>::Nodes()
{
  return Membrane<distype>::Nodes();
}

/*----------------------------------------------------------------------*
 |  Get shape type of element (private)                    sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
CORE::FE::CellType DRT::ELEMENTS::MembraneScatra<distype>::Shape() const
{
  return Membrane<distype>::Shape();
}

template class DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::tri6>;
template class DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::quad4>;
template class DRT::ELEMENTS::MembraneScatra<CORE::FE::CellType::quad9>;

FOUR_C_NAMESPACE_CLOSE
