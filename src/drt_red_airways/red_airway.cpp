
/*!----------------------------------------------------------------------
\file red_airway.cpp
\brief

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/
#ifdef D_RED_AIRWAY
#ifdef CCADISCRET

#include "red_airway.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 |  ctor (public)                                           ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAirway::RedAirway(int id, int owner) :
DRT::Element(id,element_red_airway,owner),
is_ale_(false),
data_()
{
  gaussrule_ = intrule1D_undefined;

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAirway::RedAirway(const DRT::ELEMENTS::RedAirway& old) :
DRT::Element(old),
gaussrule_(old.gaussrule_),
is_ale_(old.is_ale_),
data_(old.data_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of RedAirway and return pointer             |
 |  to it                                                      (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::RedAirway::Clone() const
{
  DRT::ELEMENTS::RedAirway* newelement = new DRT::ELEMENTS::RedAirway(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::RedAirway::Shape() const
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
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAirway::Pack(vector<char>& data) const
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
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAirway::Unpack(const vector<char>& data)
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
 |  dtor (public)                                           ismail 01/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAirway::~RedAirway()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                             ismail 01/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAirway::Print(ostream& os) const
{
  os << "RedAirway ";
  Element::Print(os);

  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return RedAirway2Register (public)            ismail 01/10|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::RedAirway::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::RedAirwayRegister(Type()));
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           ismail 01/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAirwayRegister::RedAirwayRegister(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      ismail 01/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAirwayRegister::RedAirwayRegister(
                               const DRT::ELEMENTS::RedAirwayRegister& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAirwayRegister* DRT::ELEMENTS::RedAirwayRegister::Clone() const
{
  return new DRT::ELEMENTS::RedAirwayRegister(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAirwayRegister::Pack(vector<char>& data) const
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
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAirwayRegister::Unpack(const vector<char>& data)
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
 |  dtor (public)                                          ismail 01/10 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAirwayRegister::~RedAirwayRegister()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAirwayRegister::Print(ostream& os) const
{
  os << "RedAirwayRegister ";
  ElementRegister::Print(os);
  return;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_RED_AIRWAYS
