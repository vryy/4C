/*!----------------------------------------------------------------------
\file fluid3.cpp
\brief

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "fluid3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*/
// map to convert strings tao actions (stabilisation)
/*----------------------------------------------------------------------*/
map<string,DRT::ELEMENTS::Fluid3::StabilisationAction> DRT::ELEMENTS::Fluid3::stabstrtoact_;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 02/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3::Fluid3(int id, int owner) :
DRT::Element(id,element_fluid3,owner),
is_ale_(false),
data_()
{
    gaussrule_ = intrule3D_undefined;
    surfaces_.resize(0);
    surfaceptrs_.resize(0);
    lines_.resize(0);
    lineptrs_.resize(0);

    Cs_delta_sq_=0;

    sub_acc_old_.resize(0,0);
    sub_vel_.resize(0,0);
    sub_vel_old_.resize(0,0);
    sub_pre_.resize(0);
    sub_pre_old_.resize(0);

    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 02/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3::Fluid3(const DRT::ELEMENTS::Fluid3& old) :
DRT::Element(old),
gaussrule_(old.gaussrule_),
is_ale_(old.is_ale_),
data_(old.data_),
surfaces_(old.surfaces_),
surfaceptrs_(old.surfaceptrs_),
lines_(old.lines_),
lineptrs_(old.lineptrs_),
Cs_delta_sq_(old.Cs_delta_sq_)
{
    return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Fluid3 and return pointer to it (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Fluid3::Clone() const
{
  DRT::ELEMENTS::Fluid3* newelement = new DRT::ELEMENTS::Fluid3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Fluid3::Shape() const
{
  switch (NumNode())
  {
  case  4: return tet4;
  case  5: return pyramid5;
  case  6: return wedge6;
  case  8: return hex8;
  case 10: return tet10;
  case 15: return wedge15;
  case 20: return hex20;
  case 27: return hex27;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3::Pack(vector<char>& data) const
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
  // Cs_delta_sq_, the Smagorinsky constant for the dynamic Smagorinsky model
  AddtoPack(data,Cs_delta_sq_);

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
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3::Unpack(const vector<char>& data)
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
  gaussrule_ = GaussRule3D(gausrule_integer); //explicit conversion from integer to enum
  // is_ale_
  ExtractfromPack(position,data,is_ale_);
  // extract Cs_delta_sq_, the Smagorinsky constant for the dynamic
  // Smagorinsky model
  ExtractfromPack(position,data,Cs_delta_sq_);

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

  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 02/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3::~Fluid3()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 02/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3::Print(ostream& os) const
{
  os << "Fluid3 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}


/*----------------------------------------------------------------------*
 |  allocate and return Fluid3Register (public)              mwgee 02/08|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Fluid3::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Fluid3Register(Type()));
}


/*----------------------------------------------------------------------*
 |  get vector of lines              (public)                  gjb 03/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::Fluid3::Lines()
{
  // once constructed do not reconstruct again
  // make sure they exist
  if ((int)lines_.size()    == NumLine() &&
      (int)lineptrs_.size() == NumLine() &&
      dynamic_cast<DRT::ELEMENTS::Fluid3Line*>(lineptrs_[0]) )
    return (DRT::Element**)(&(lineptrs_[0]));
  
  // so we have to allocate new line elements
  DRT::UTILS::ElementBoundaryFactory<Fluid3Line,Fluid3>(false,lines_,lineptrs_,this);

  return (DRT::Element**)(&(lineptrs_[0]));
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                            gjb 05/08|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::Fluid3::Surfaces()
{
  // once constructed do not reconstruct again
  // make sure they exist
  if ((int)surfaces_.size()    == NumSurface() &&
      (int)surfaceptrs_.size() == NumSurface() &&
      dynamic_cast<DRT::ELEMENTS::Fluid3Surface*>(surfaceptrs_[0]) )
    return (DRT::Element**)(&(surfaceptrs_[0]));

  // so we have to allocate new surface elements
  DRT::UTILS::ElementBoundaryFactory<Fluid3Surface,Fluid3>(false,surfaces_,surfaceptrs_,this);
  
  return (DRT::Element**)(&(surfaceptrs_[0]));
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                g.bau 03/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::Fluid3::Volumes()
{
  volume_.resize(1);
  volume_[0] = this; //points to Fluid3 element itself
  return &volume_[0];
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3Register::Fluid3Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3Register::Fluid3Register(
                               const DRT::ELEMENTS::Fluid3Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3Register* DRT::ELEMENTS::Fluid3Register::Clone() const
{
  return new DRT::ELEMENTS::Fluid3Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Register::Pack(vector<char>& data) const
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
void DRT::ELEMENTS::Fluid3Register::Unpack(const vector<char>& data)
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
DRT::ELEMENTS::Fluid3Register::~Fluid3Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Register::Print(ostream& os) const
{
  os << "Fluid3Register ";
  ElementRegister::Print(os);
  return;
}





#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
