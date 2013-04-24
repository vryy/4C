/*----------------------------------------------------------------------*/
/*!
 \file so3_poro_p1.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/

#include "so3_poro_p1.H"
#include "so3_poro_p1_eletypes.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            vuong 03/12|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Poro_P1<so3_ele,distype>::So3_Poro_P1(int id, int owner):
So3_Poro<so3_ele,distype>(id,owner)
{
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       vuong 03/12|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Poro_P1<so3_ele,distype>::So3_Poro_P1(const DRT::ELEMENTS::So3_Poro_P1<so3_ele,distype>& old):
So3_Poro<so3_ele,distype>(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            vuong 03/12|
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
DRT::Element* DRT::ELEMENTS::So3_Poro_P1<so3_ele,distype>::Clone() const
{
  DRT::ELEMENTS::So3_Poro_P1< so3_ele, distype>* newelement =
      new DRT::ELEMENTS::So3_Poro_P1< so3_ele, distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele,distype>::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  so3_ele::AddtoPack(data,type);

  // add base class Element
  my::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele,distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  so3_ele::ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  my::ExtractfromPack(position,data,basedata);
  my::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              vuong 03/12|
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele,distype>::Print(ostream& os) const
{
  os << "So3_Poro_P1 ";
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Poro_P1<so3_ele,distype>::UniqueParObjectId() const
{
  switch(distype)
  {
  case DRT::Element::hex8:
    return So_hex8PoroP1Type::Instance().UniqueParObjectId();
    break;
  default: dserror("unknown element type!");
    break;
  }
  return -1;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ElementType & DRT::ELEMENTS::So3_Poro_P1<so3_ele,distype>::ElementType() const
{
  switch(distype)
  {
  //case DRT::Element::tet4:
  //  return So_tet4PoroType::Instance();
  //case DRT::Element::tet10:
  //  return So_tet10PoroType::Instance();
  case DRT::Element::hex8:
    return So_hex8PoroP1Type::Instance();
  //case DRT::Element::hex27:
  //  return So_hex27PoroType::Instance();
  default: dserror("unknown element type!");
    break;
  }
  return So_hex8PoroP1Type::Instance();
};

#include "so3_poro_p1_fwd.hpp"
