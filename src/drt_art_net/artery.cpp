
/*!----------------------------------------------------------------------
\file artery.cpp
\brief

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ARTNET
#ifdef CCADISCRET

#include "artery.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 |  ctor (public)                                           ismail 01/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Artery::Artery(int id, int owner) :
DRT::Element(id,element_artery,owner),
is_ale_(false),
data_()
{
  gaussrule_ = intrule1D_undefined;

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      ismail 01/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Artery::Artery(const DRT::ELEMENTS::Artery& old) :
DRT::Element(old),
gaussrule_(old.gaussrule_),
is_ale_(old.is_ale_),
data_(old.data_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Artery and return pointer to it (public) |
 |                                                         ismail 01/09 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Artery::Clone() const
{
  DRT::ELEMENTS::Artery* newelement = new DRT::ELEMENTS::Artery(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                         ismail 01/09 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Artery::Shape() const
{
  switch (NumNode())
  {
  case  2: return line2;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         ismail 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Artery::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // Gaussrule
  AddtoPack(data,gaussrule_); //implicit conversion from enum to integer
  // is_ale_
  AddtoPack(data,is_ale_);


  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         ismail 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Artery::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);

  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class Element
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);
  // Gaussrule
  int gausrule_integer;
  ExtractfromPack(position,data,gausrule_integer);
  gaussrule_ = GaussRule1D(gausrule_integer); //explicit conversion from integer to enum
  // is_ale_
  ExtractfromPack(position,data,is_ale_);

  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                           ismail 01/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Artery::~Artery()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                             ismail 01/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Artery::Print(ostream& os) const
{
  os << "Artery ";
  Element::Print(os);

  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Artery2Register (public)            ismail 01/09|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Artery::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::ArteryRegister(Type()));
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           ismail 01/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ArteryRegister::ArteryRegister(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      ismail 01/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ArteryRegister::ArteryRegister(
                               const DRT::ELEMENTS::ArteryRegister& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                         ismail 01/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ArteryRegister* DRT::ELEMENTS::ArteryRegister::Clone() const
{
  return new DRT::ELEMENTS::ArteryRegister(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         ismail 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ArteryRegister::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class ElementRegister
  vector<char> basedata(0);
  ElementRegister::Pack(basedata);
  AddtoPack(data,basedata);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         ismail 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ArteryRegister::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // base class ElementRegister
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  ElementRegister::Unpack(basedata);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                          ismail 01/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ArteryRegister::~ArteryRegister()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                         ismail 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ArteryRegister::Print(ostream& os) const
{
  os << "ArteryRegister ";
  ElementRegister::Print(os);
  return;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_ARTNET
