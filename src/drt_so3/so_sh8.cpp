/*!----------------------------------------------------------------------
\file so_sh8.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOH8
#ifdef CCADISCRET

#include "so_sh8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "so_hex8.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh8::So_sh8(int id, int owner) :
DRT::ELEMENTS::So_hex8(id,owner)
{
  SetType(element_sosh8);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh8::So_sh8(const DRT::ELEMENTS::So_sh8& old) :
DRT::ELEMENTS::So_hex8(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_sh8::Clone() const
{
  DRT::ELEMENTS::So_sh8* newelement = new DRT::ELEMENTS::So_sh8(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class So_hex8 Element
  vector<char> basedata(0);
  DRT::ELEMENTS::So_hex8::Pack(basedata);
  AddtoPack(data,basedata);
  // thickdir
  AddtoPack(data,thickdir_);
  //AddtoPack(data,thickvec_);
  AddtoPack(data,nodes_rearranged_);
  // original (input) nodeids
  AddtoPack(data,inp_nodeIds_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8::Unpack(const vector<char>& data)
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
  // thickdir
  ExtractfromPack(position,data,thickdir_);
  //ExtractfromPack(position,data,thickvec_);
  ExtractfromPack(position,data,nodes_rearranged_);
  ExtractfromPack(position,data,inp_nodeIds_);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh8::~So_sh8()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8::Print(ostream& os) const
{
  os << "So_sh8 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Sosh8Register (public)                 maf 04/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::So_sh8::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Sosh8Register(Type()));
}

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sosh8Register::Sosh8Register(DRT::Element::ElementType etype) :
DRT::ELEMENTS::Soh8Register::Soh8Register(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sosh8Register::Sosh8Register(
                               const DRT::ELEMENTS::Sosh8Register& old) :
DRT::ELEMENTS::Soh8Register::Soh8Register(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sosh8Register* DRT::ELEMENTS::Sosh8Register::Clone() const
{
  return new DRT::ELEMENTS::Sosh8Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Sosh8Register::Pack(vector<char>& data) const
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
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Sosh8Register::Unpack(const vector<char>& data)
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
 |  print (public)                                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Sosh8Register::Print(ostream& os) const
{
  os << "Sosh8Register ";
  ElementRegister::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sosh8Register::~Sosh8Register()
{
  return;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
