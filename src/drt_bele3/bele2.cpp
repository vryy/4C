/*!----------------------------------------------------------------------
\file bele2.cpp

\brief  2D boundary elment

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "bele2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele2::Bele2(int id, int owner) :
DRT::Element(id,element_bele2,owner)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele2::Bele2(const DRT::ELEMENTS::Bele2& old) :
DRT::Element(old)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Bele2::Clone() const
{
  DRT::ELEMENTS::Bele2* newelement = new DRT::ELEMENTS::Bele2(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Bele2::Shape() const
{
  switch (NumNode())
  {
  case  2: return line2;
  case  3: return line3;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele2::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele2::Unpack(const vector<char>& data)
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

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele2::~Bele2()
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele2::Print(ostream& os) const
{
  os << "Bele2 " << DRT::DistypeToString(Shape());
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Bele2::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Bele2Register(Type()));
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gjb 05/08|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Bele2::Lines()
{
  vector<RCP<DRT::Element> > lines(1);
  lines[0]=rcp(this,false);
  return lines;
}


/*----------------------------------------------------------------------*
 |  get vector of Surfaces (length 1) (public)               gammi 04/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Bele2::Surfaces()
{
  vector<RCP<DRT::Element> > surfaces(0);
  return surfaces;
}



//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                            u.may 12/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele2Register::Bele2Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       u.may 12/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele2Register::Bele2Register(
                               const DRT::ELEMENTS::Bele2Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                          u.may 12/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele2Register* DRT::ELEMENTS::Bele2Register::Clone() const
{
  return new DRT::ELEMENTS::Bele2Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          u.may 02/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele2Register::Pack(vector<char>& data) const
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
 |                                                          u.may 02/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele2Register::Unpack(const vector<char>& data)
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
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            u.may 12/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele2Register::~Bele2Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           u.may 12/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele2Register::Print(ostream& os) const
{
  os << "Bele2Register ";
  ElementRegister::Print(os);
  return;
}





#endif  // #ifdef CCADISCRET
