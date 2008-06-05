/*!----------------------------------------------------------------------
\file shell8.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#ifdef CCADISCRET

#include "shell8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"




/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Shell8::Shell8(int id, int owner) :
DRT::Element(id,element_shell8,owner),
forcetype_(s8_none),
thickness_(0.0),
ngptri_(0),
nhyb_(0),
ans_(0),
sdc_(1.0),
material_(0),
data_()
{
  ngp_[0] = ngp_[1] = ngp_[2] = 0;
  eas_[0] = eas_[1] = eas_[2] = eas_[3] = eas_[4] = 0;
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Shell8::Shell8(const DRT::ELEMENTS::Shell8& old) :
DRT::Element(old),
forcetype_(old.forcetype_),
thickness_(old.thickness_),
ngptri_(old.ngptri_),
nhyb_(old.nhyb_),
ans_(old.ans_),
sdc_(old.sdc_),
material_(old.material_),
data_(old.data_)
{
  for (int i=0; i<3; ++i) ngp_[i] = old.ngp_[i];
  for (int i=0; i<5; ++i) eas_[i] = old.eas_[i];
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Shell8 and return pointer to it (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Shell8::Clone() const
{
  DRT::ELEMENTS::Shell8* newelement = new DRT::ELEMENTS::Shell8(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Shell8::Shape() const
{
  switch (NumNode())
  {
  case 4: return quad4;
  case 8: return quad8;
  case 9: return quad9;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Shell8::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // forcetype_
  AddtoPack(data,forcetype_);
  // thickness_
  AddtoPack(data,thickness_);
  // ngp_
  AddtoPack(data,ngp_,3*sizeof(int));
  // ngptri_
  AddtoPack(data,ngptri_);
  // nhyb_
  AddtoPack(data,nhyb_);
  // eas_
  AddtoPack(data,eas_,5*sizeof(int));
  // ans_
  AddtoPack(data,ans_);
  // sdc_
  AddtoPack(data,sdc_);
  // material_
  AddtoPack(data,material_);
  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Shell8::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);
  // forcetype_
  ExtractfromPack(position,data,forcetype_);
  // thickness_
  ExtractfromPack(position,data,thickness_);
  // ngp_
  ExtractfromPack(position,data,ngp_,3*sizeof(int));
  // ngptri_
  ExtractfromPack(position,data,ngptri_);
  // nhyb_
  ExtractfromPack(position,data,nhyb_);
  // eas_
  ExtractfromPack(position,data,eas_,5*sizeof(int));
  // ans_
  ExtractfromPack(position,data,ans_);
  // sdc_
  ExtractfromPack(position,data,sdc_);
  // material_
  ExtractfromPack(position,data,material_);
  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Shell8::~Shell8()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Shell8::Print(ostream& os) const
{
  os << "Shell8 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Shell8Register (public)              mwgee 12/06|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Shell8::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Shell8Register(Type()));
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             mwgee 01/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Shell8::Lines()
{
  // do NOT store line or surface elements inside the parent element 
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization, 
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Shell8Line,Shell8>(DRT::UTILS::buildLines,this);
}

/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          mwgee 01/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Shell8::Surfaces()
{
  vector<RCP<Element> > surfaces(1);
  surfaces[0]= rcp(this, false);
  return surfaces;
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Shell8Register::Shell8Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Shell8Register::Shell8Register(
                               const DRT::ELEMENTS::Shell8Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Shell8Register* DRT::ELEMENTS::Shell8Register::Clone() const
{
  return new DRT::ELEMENTS::Shell8Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Shell8Register::Pack(vector<char>& data) const
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
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Shell8Register::Unpack(const vector<char>& data)
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
 |  dtor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Shell8Register::~Shell8Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Shell8Register::Print(ostream& os) const
{
  os << "Shell8Register ";
  ElementRegister::Print(os);
  return;
}





#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SHELL8
