/*!----------------------------------------------------------------------

\level 3

\maintainer Sebastian Fuchs
            fuchs@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262

\brief Nonlinear Membrane Finite Element with ScaTra coupling

*----------------------------------------------------------------------*/

#include "membrane_scatra.H"

#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |  constructor (public)                                   sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::MembraneScatra<distype>::MembraneScatra(int id, int owner)
    : Membrane<distype>(id, owner), impltype_(INPAR::SCATRA::impltype_undefined)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-constructor (public)                              sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
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
template <DRT::Element::DiscretizationType distype>
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
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MembraneScatra<distype>::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  DRT::ParObject::AddtoPack(data, type);

  // pack scalar transport impltype_
  DRT::ParObject::AddtoPack(data, impltype_);

  // add base class Element
  Membrane<distype>::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MembraneScatra<distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  DRT::ParObject::ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // extract scalar transport impltype
  impltype_ = static_cast<INPAR::SCATRA::ImplType>(DRT::ParObject::ExtractInt(position, data));

  // extract base class Element
  std::vector<char> basedata(0);
  DRT::ParObject::ExtractfromPack(position, data, basedata);
  Membrane<distype>::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                            sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MembraneScatra<distype>::Print(std::ostream& os) const
{
  os << "MembraneScatra ";
  Membrane<distype>::Print(os);

  return;
}

/*----------------------------------------------------------------------*
 |  read this element (public)                             sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::MembraneScatra<distype>::ReadElement(
    const std::string& eletype, const std::string& eledistype, DRT::INPUT::LineDefinition* linedef)
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
  else if (impltype == "BondReac")
    impltype_ = INPAR::SCATRA::impltype_bondreac;
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
    dserror("Invalid implementation type for Wall1_Scatra elements!");

  return true;
}

/*----------------------------------------------------------------------*
 |  Get vector of ptrs to nodes (private)                  sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
inline DRT::Node** DRT::ELEMENTS::MembraneScatra<distype>::Nodes()
{
  return Membrane<distype>::Nodes();
}

/*----------------------------------------------------------------------*
 |  Get shape type of element (private)                    sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::Element::DiscretizationType DRT::ELEMENTS::MembraneScatra<distype>::Shape() const
{
  return Membrane<distype>::Shape();
}

template class DRT::ELEMENTS::MembraneScatra<DRT::Element::tri3>;
template class DRT::ELEMENTS::MembraneScatra<DRT::Element::tri6>;
template class DRT::ELEMENTS::MembraneScatra<DRT::Element::quad4>;
template class DRT::ELEMENTS::MembraneScatra<DRT::Element::quad9>;
