/*----------------------------------------------------------------------*/
/*!
 \file wall1_poro_p1_scatra.cpp

 \brief a 2D wall element for solid-part of porous medium using p1 (mixed) approach including scatra
 functionality

 \level 2

 \maintainer Christoph Schmidt
             schmidt@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 15251
 *----------------------------------------------------------------------*/

#include "wall1_poro_p1_scatra.H"
#include "wall1_poro_p1_scatra_eletypes.H"

#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Wall1_PoroP1Scatra<distype>::Wall1_PoroP1Scatra(int id, int owner)
    : DRT::ELEMENTS::Wall1_PoroP1<distype>(id, owner), impltype_(INPAR::SCATRA::impltype_undefined)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                    schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Wall1_PoroP1Scatra<distype>::Wall1_PoroP1Scatra(
    const DRT::ELEMENTS::Wall1_PoroP1Scatra<distype>& old)
    : DRT::ELEMENTS::Wall1_PoroP1<distype>(old), impltype_(old.impltype_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance and return pointer to it (public)           |
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::Element* DRT::ELEMENTS::Wall1_PoroP1Scatra<distype>::Clone() const
{
  DRT::ELEMENTS::Wall1_PoroP1Scatra<distype>* newelement =
      new DRT::ELEMENTS::Wall1_PoroP1Scatra<distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data (public)                                    schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1Scatra<distype>::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
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
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1Scatra<distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  my::ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract scalar transport impltype
  impltype_ = static_cast<INPAR::SCATRA::ImplType>(my::ExtractInt(position, data));

  // extract base class Element
  std::vector<char> basedata(0);
  my::ExtractfromPack(position, data, basedata);
  my::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}

/*----------------------------------------------------------------------*
 |  read this element (public)                             schmidt 09/17|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::Wall1_PoroP1Scatra<distype>::ReadElement(
    const std::string& eletype, const std::string& eledistype, DRT::INPUT::LineDefinition* linedef)
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
    dserror("Invalid implementation type for Wall1_PoroP1Scatra elements!");

  return true;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                           schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1Scatra<distype>::Print(std::ostream& os) const
{
  os << "Wall1_PoroP1Scatra ";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << my::data_;
  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                           schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Wall1_PoroP1Scatra<distype>::UniqueParObjectId() const
{
  int parobjectid(-1);
  switch (distype)
  {
    case DRT::Element::tri3:
    {
      parobjectid = DRT::ELEMENTS::WallTri3PoroP1ScatraType::Instance().UniqueParObjectId();
      break;
    }
    case DRT::Element::quad4:
    {
      parobjectid = DRT::ELEMENTS::WallQuad4PoroP1ScatraType::Instance().UniqueParObjectId();
      break;
    }
    case DRT::Element::quad9:
    {
      parobjectid = DRT::ELEMENTS::WallQuad9PoroP1ScatraType::Instance().UniqueParObjectId();
      break;
    }
    default:
    {
      dserror("unknown element type");
      break;
    }
  }
  return parobjectid;
}

/*----------------------------------------------------------------------*
 | get the element type (public)                           schmidt 09/17|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ElementType& DRT::ELEMENTS::Wall1_PoroP1Scatra<distype>::ElementType() const
{
  switch (distype)
  {
    case DRT::Element::tri3:
      return DRT::ELEMENTS::WallTri3PoroP1ScatraType::Instance();
      break;
    case DRT::Element::quad4:
      return DRT::ELEMENTS::WallQuad4PoroP1ScatraType::Instance();
      break;
    case DRT::Element::quad9:
      return DRT::ELEMENTS::WallQuad9PoroP1ScatraType::Instance();
      break;
    default:
      dserror("unknown element type");
      break;
  }
  return DRT::ELEMENTS::WallQuad4PoroP1ScatraType::Instance();
}

/*----------------------------------------------------------------------*
 *                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
template class DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::tri3>;
template class DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::quad4>;
template class DRT::ELEMENTS::Wall1_PoroP1Scatra<DRT::Element::quad9>;
