/*----------------------------------------------------------------------*/
/*!
 \file so3_poro_p1_scatra.cpp

 \brief implementation of the 3D solid-poro element (p1, mixed approach) including scatra
 functionality

 \level 2

 \maintainer Christoph Schmidt
             schmidt@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 15251
 *----------------------------------------------------------------------*/

#include "so3_poro_p1_scatra.H"
#include "so3_poro_p1_scatra_eletypes.H"

#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                         schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Poro_P1_Scatra<so3_ele, distype>::So3_Poro_P1_Scatra(int id, int owner)
    : So3_Poro_P1<so3_ele, distype>(id, owner), impltype_(INPAR::SCATRA::impltype_undefined)
{
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                    schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Poro_P1_Scatra<so3_ele, distype>::So3_Poro_P1_Scatra(
    const DRT::ELEMENTS::So3_Poro_P1_Scatra<so3_ele, distype>& old)
    : So3_Poro_P1<so3_ele, distype>(old), impltype_(old.impltype_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance and return pointer to it (public)           |
 |                                                        schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
DRT::Element* DRT::ELEMENTS::So3_Poro_P1_Scatra<so3_ele, distype>::Clone() const
{
  DRT::ELEMENTS::So3_Poro_P1_Scatra<so3_ele, distype>* newelement =
      new DRT::ELEMENTS::So3_Poro_P1_Scatra<so3_ele, distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data (public)                                    schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1_Scatra<so3_ele, distype>::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
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
template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1_Scatra<so3_ele, distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  so3_ele::ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // extract scalar transport impltype_
  impltype_ = static_cast<INPAR::SCATRA::ImplType>(so3_ele::ExtractInt(position, data));

  // extract base class Element
  std::vector<char> basedata(0);
  my::ExtractfromPack(position, data, basedata);
  my::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                           schmidt 09/17 |
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1_Scatra<so3_ele, distype>::Print(std::ostream& os) const
{
  os << "So3_Poro_P1_Scatra ";
  os << DRT::DistypeToString(distype).c_str() << " ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 | get the unique ParObject Id (public)                    schmidt 09/17|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Poro_P1_Scatra<so3_ele, distype>::UniqueParObjectId() const
{
  int parobjectid(-1);
  switch (distype)
  {
    case DRT::Element::hex8:
    {
      parobjectid = So_hex8PoroP1ScatraType::Instance().UniqueParObjectId();
      break;
    }
    case DRT::Element::tet4:
    {
      parobjectid = So_tet4PoroP1ScatraType::Instance().UniqueParObjectId();
      break;
    }
    default:
    {
      dserror("unknown element type!");
      break;
    }
  }
  return parobjectid;
}

/*----------------------------------------------------------------------*
 | get the element type (public)                           schmidt 09/17|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ElementType& DRT::ELEMENTS::So3_Poro_P1_Scatra<so3_ele, distype>::ElementType() const
{
  switch (distype)
  {
    case DRT::Element::tet4:
      return So_tet4PoroP1ScatraType::Instance();
    case DRT::Element::hex8:
      return So_hex8PoroP1ScatraType::Instance();
    default:
      dserror("unknown element type!");
      break;
  }
  return So_hex8PoroP1ScatraType::Instance();
}

/*----------------------------------------------------------------------*
 |  read this element (public)                             schmidt 09/17|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::So3_Poro_P1_Scatra<so3_ele, distype>::ReadElement(
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
    dserror("Invalid implementation type for So3_Poro_P1_Scatra elements!");

  return true;
}

/*----------------------------------------------------------------------*
 |                                                         schmidt 09/17|
 *----------------------------------------------------------------------*/
template class DRT::ELEMENTS::So3_Poro_P1_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>;
template class DRT::ELEMENTS::So3_Poro_P1_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>;
