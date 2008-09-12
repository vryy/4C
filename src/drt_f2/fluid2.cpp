/*!----------------------------------------------------------------------
\file fluid2.cpp
\brief

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET

#include "fluid2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"

using namespace DRT::UTILS;

map<string,DRT::ELEMENTS::Fluid2::StabilisationAction> DRT::ELEMENTS::Fluid2::stabstrtoact_;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2::Fluid2(int id, int owner) :
DRT::Element(id,element_fluid2,owner),
is_ale_(false),
data_()
{
  gaussrule_ = intrule2D_undefined;
  
  sub_acc_old_.resize(0,0);
  sub_vel_.resize(0,0);
  sub_vel_old_.resize(0,0);
  sub_pre_.resize(0);
  sub_pre_old_.resize(0);
  
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2::Fluid2(const DRT::ELEMENTS::Fluid2& old) :
DRT::Element(old),
gaussrule_(old.gaussrule_),
is_ale_(old.is_ale_),
data_(old.data_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Fluid2 and return pointer to it (public) |
 |                                                          gammi 11/06 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Fluid2::Clone() const
{
  DRT::ELEMENTS::Fluid2* newelement = new DRT::ELEMENTS::Fluid2(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Fluid2::Shape() const
{
  switch (NumNode())
  {
  case  3: return tri3;
  case  4: return quad4;
  case  6: return tri6;
  case  8: return quad8;
  case  9: return quad9;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          gammi 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2::Pack(vector<char>& data) const
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

  // history variables
  AddtoPack(data,sub_acc_old_.extent(blitz::firstDim));
  AddtoPack(data,sub_acc_old_.extent(blitz::secondDim));

  int size = sub_acc_old_.extent(blitz::firstDim)
             *sub_acc_old_.extent(blitz::secondDim)
             *sizeof(double);
  AddtoPack(data,sub_acc_old_.data(),size);
  AddtoPack(data,sub_vel_.data()    ,size);
  AddtoPack(data,sub_vel_old_.data(),size);

  size = sub_acc_old_.extent(blitz::secondDim)*sizeof(double);
  AddtoPack(data,sub_pre_.data()    ,size);
  AddtoPack(data,sub_pre_old_.data(),size);


  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          gammi 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2::Unpack(const vector<char>& data)
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
  gaussrule_ = GaussRule2D(gausrule_integer); //explicit conversion from integer to enum
  // is_ale_
  ExtractfromPack(position,data,is_ale_);

  // history variables (subscale velocities, accelerations and pressure)
  {
    int firstdim;
    int secondim;

    ExtractfromPack(position,data,firstdim);
    ExtractfromPack(position,data,secondim);

    sub_acc_old_.resize(firstdim,secondim);
    sub_vel_    .resize(firstdim,secondim);
    sub_vel_old_.resize(firstdim,secondim);

    int size = firstdim*secondim*sizeof(double);

    ExtractfromPack(position,data,&(sub_acc_old_.data()[0]),size);
    ExtractfromPack(position,data,&(sub_vel_.data()[0])    ,size);
    ExtractfromPack(position,data,&(sub_vel_old_.data()[0]),size);

    sub_pre_    .resize(secondim);
    sub_pre_old_.resize(secondim);

    ExtractfromPack(position,data,&(sub_pre_.data()[0])    ,secondim*sizeof(double));
    ExtractfromPack(position,data,&(sub_pre_old_.data()[0]),secondim*sizeof(double));
  }

  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 11/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2::~Fluid2()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 11/06|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2::Print(ostream& os) const
{
  os << "Fluid2 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Fluid2Register (public)              gammi 04/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Fluid2::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Fluid2Register(Type()));
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gjb 05/08|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> >  DRT::ELEMENTS::Fluid2::Lines()
{
  // do NOT store line or surface elements inside the parent element 
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization, 
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Fluid2Line,Fluid2>(DRT::UTILS::buildLines,this);
}


/*----------------------------------------------------------------------*
 |  get vector of Surfaces (length 1) (public)               gammi 04/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> >  DRT::ELEMENTS::Fluid2::Surfaces()
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
DRT::ELEMENTS::Fluid2Register::Fluid2Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2Register::Fluid2Register(
                               const DRT::ELEMENTS::Fluid2Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2Register* DRT::ELEMENTS::Fluid2Register::Clone() const
{
  return new DRT::ELEMENTS::Fluid2Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2Register::Pack(vector<char>& data) const
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
void DRT::ELEMENTS::Fluid2Register::Unpack(const vector<char>& data)
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
DRT::ELEMENTS::Fluid2Register::~Fluid2Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2Register::Print(ostream& os) const
{
  os << "Fluid2Register ";
  ElementRegister::Print(os);
  return;
}





#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2
