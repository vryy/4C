/*----------------------------------------------------------------------*/
/*!
\file so_sh8p8.cpp
\brief

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* defintions */
#ifdef D_SOLID3
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "so_sh8p8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            bborn 03/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh8p8::So_sh8p8(int id, int owner) :
DRT::ELEMENTS::So_sh8(id,owner)
{
  SetType(element_sosh8p8);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       bborn 03/09|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh8p8::So_sh8p8(const DRT::ELEMENTS::So_sh8p8& old) :
DRT::ELEMENTS::So_sh8(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                          bborn 03/09 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_sh8p8::Clone() const
{
  DRT::ELEMENTS::So_sh8p8* newelement = new DRT::ELEMENTS::So_sh8p8(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          bborn 03/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Pack(std::vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class So_sh8 Element
  std::vector<char> basedata(0);
  DRT::ELEMENTS::So_sh8::Pack(basedata);
  AddtoPack(data,basedata);
  // techniques
  AddtoPack(data,stab_);
  AddtoPack(data,ans_);
  AddtoPack(data,lin_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          bborn 03/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Unpack(const std::vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class So_sh8 Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  DRT::ELEMENTS::So_sh8::Unpack(basedata);
  // techniques
  ExtractfromPack(position,data,stab_);
  ExtractfromPack(position,data,ans_);
  ExtractfromPack(position,data,lin_);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            bborn 03/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh8p8::~So_sh8p8()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              bborn 03/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Print(std::ostream& os) const
{
  os << "So_sh8p8 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Sosh8p8Register (public)             bborn 03/09|
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ElementRegister> DRT::ELEMENTS::So_sh8p8::ElementRegister() const
{
  return Teuchos::rcp(new DRT::ELEMENTS::Sosh8p8Register(Type()));
}

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                            bborn 03/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sosh8p8Register::Sosh8p8Register(
  DRT::Element::ElementType etype
  )
  : DRT::ELEMENTS::Sosh8Register::Sosh8Register(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       bborn 03/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sosh8p8Register::Sosh8p8Register(
  const DRT::ELEMENTS::Sosh8p8Register& old
  )
  : DRT::ELEMENTS::Sosh8Register::Sosh8Register(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                          bborn 03/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sosh8p8Register* DRT::ELEMENTS::Sosh8p8Register::Clone() const
{
  return new DRT::ELEMENTS::Sosh8p8Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          bborn 03/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Sosh8p8Register::Pack(std::vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class So_sh8 Element
  std::vector<char> basedata(0);
  DRT::ELEMENTS::Sosh8Register::Pack(basedata);
  AddtoPack(data,basedata);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          bborn 03/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Sosh8p8Register::Unpack(const std::vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // base class ElementRegister
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  DRT::ELEMENTS::Sosh8Register::Unpack(basedata);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           bborn 03/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Sosh8p8Register::Print(std::ostream& os) const
{
  os << "Sosh8p8Register ";
  ElementRegister::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            bborn 03/09|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sosh8p8Register::~Sosh8p8Register()
{
  return;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
