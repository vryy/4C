#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_hex8p1j1.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "so_hex8.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                               lw 12/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_Hex8P1J1::So_Hex8P1J1(int id, int owner) :
DRT::ELEMENTS::So_hex8(id,owner)
{
  SetType(element_so_hex8p1j1);

  K_pu_.PutScalar(0.0);
  K_tu_.PutScalar(0.0);

  R_t_.PutScalar(0.0);
  R_p_.PutScalar(0.0);

  K_tt_ = 0.0;
  K_pt_ = 0.0;

  p_.PutScalar(0.0);
  p_o_.PutScalar(0.0);

  t_.PutScalar(1.0);
  t_o_.PutScalar(1.0);

  m_.PutScalar(0.0);
  for (int i=0; i<3; ++i)
  {
    m_(i,0)=1.0;
  }

  Identity6_.PutScalar(0.0);
  for (int i=0; i<6; ++i)
  {
    Identity6_(i,i) = 1.0;
  }

  I_d_ = Identity6_;
  I_d_.MultiplyNT(-1.0/3.0, m_, m_, 1.0);

  I_0_.PutScalar(0.0);

  for (int i=0; i<3; ++i)
  {
    I_0_(i,i) = 1.0;
  }
  for (int i=3; i<6; ++i)
  {
    I_0_(i,i) = 0.5;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                          lw 12/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_Hex8P1J1::So_Hex8P1J1(const DRT::ELEMENTS::So_Hex8P1J1& old) :
DRT::ELEMENTS::So_hex8(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                              lw 12/08|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_Hex8P1J1::Clone() const
{
  DRT::ELEMENTS::So_Hex8P1J1* newelement = new DRT::ELEMENTS::So_Hex8P1J1(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                              lw 12/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_Hex8P1J1::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class So_hex8 Element
  vector<char> basedata(0);
  DRT::ELEMENTS::So_hex8::Pack(basedata);
  AddtoPack(data,basedata);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                              lw 12/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_Hex8P1J1::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class So_hex8 Element
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  DRT::ELEMENTS::So_hex8::Unpack(basedata);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                               lw 12/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_Hex8P1J1::~So_Hex8P1J1()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                 lw 12/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_Hex8P1J1::Print(ostream& os) const
{
  os << "So_Hex8P1J1 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return SoHex8P1J1Register (public)             lw 12/08|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::So_Hex8P1J1::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::SoHex8P1J1Register(Type()));
}

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                               lw 12/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoHex8P1J1Register::SoHex8P1J1Register(DRT::Element::ElementType etype) :
DRT::ELEMENTS::Soh8Register::Soh8Register(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                          lw 12/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoHex8P1J1Register::SoHex8P1J1Register(const DRT::ELEMENTS::SoHex8P1J1Register& old) :
DRT::ELEMENTS::Soh8Register::Soh8Register(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                              lw 12/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoHex8P1J1Register* DRT::ELEMENTS::SoHex8P1J1Register::Clone() const
{
  return new DRT::ELEMENTS::SoHex8P1J1Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                              lw 12/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8P1J1Register::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class So_hex8 Element
  vector<char> basedata(0);
  DRT::ELEMENTS::Soh8Register::Pack(basedata);
  AddtoPack(data,basedata);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                              lw 12/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8P1J1Register::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // base class ElementRegister
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  DRT::ELEMENTS::Soh8Register::Unpack(basedata);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                              lw 12/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoHex8P1J1Register::Print(ostream& os) const
{
  os << "SoHex8P1J1Register ";
  ElementRegister::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoHex8P1J1Register::~SoHex8P1J1Register()
{
  return;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3

