/*!----------------------------------------------------------------------
\file beam2.cpp
\brief two dimensional nonlinear corotational Timoshenko beam element

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BEAM2
#ifdef CCADISCRET

#include "beam2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_elementregister.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2::Beam2(int id, int owner) :
DRT::Element(id,element_beam2,owner),
data_(),
material_(0),
crosssec_(0),
crosssecshear_(0),
mominer_(0),
lumpedflag_(0),
thermalenergy_(0),
//since lines_ is a vector it calls its constructor automatically -> 
//no initialization of lines_ necessary here

//note: for corotational approach integration for Neumann conditions only
//hence enough to integrate 3rd order polynomials exactly
gaussrule_(DRT::UTILS::intrule_line_2point)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2::Beam2(const DRT::ELEMENTS::Beam2& old) :
DRT::Element(old),
data_(old.data_),
material_(old.material_),
crosssec_(old.crosssec_),
crosssecshear_(old.crosssecshear_),
mominer_(old.mominer_),
lumpedflag_(old.lumpedflag_),
thermalenergy_(old.thermalenergy_),
lines_(old.lines_),
gaussrule_(old.gaussrule_)
{
  return;
}
/*----------------------------------------------------------------------*
 |  Deep copy this instance of Beam2 and return pointer to it (public) |
 |                                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Beam2::Clone() const
{
  DRT::ELEMENTS::Beam2* newelement = new DRT::ELEMENTS::Beam2(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2::~Beam2()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 01/08
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2::Print(ostream& os) const
{
  os << "Beam2 ";
  Element::Print(os);
  os << " gaussrule_: " << gaussrule_ << " ";
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Beam2Register (public)               cyron 01/08|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Beam2::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Beam2Register(Type()));
}



/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          cyron 01/08 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Beam2::Shape() const
{
  return line2;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 01/08/
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  //material type
  AddtoPack(data,material_);
  //cross section
  AddtoPack(data,crosssec_);
   //cross section with shear correction
  AddtoPack(data,crosssecshear_);
  //moment of inertia of area
  AddtoPack(data,mominer_);
  //flag determining if consistent or lumped mass matrix
  AddtoPack(data,lumpedflag_);
  //thermal energy responsible for statistical forces
  AddtoPack(data,thermalenergy_);
  // gaussrule_
  AddtoPack(data,gaussrule_); //implicit conversion from enum to integer
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2::Unpack(const vector<char>& data)
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
  //material type
  ExtractfromPack(position,data,material_);
  //cross section
  ExtractfromPack(position,data,crosssec_);
  //cross section with shear correction
  ExtractfromPack(position,data,crosssecshear_);
  //moment of inertia of area
  ExtractfromPack(position,data,mominer_);
  //flag determining if consistent or lumped mass matrix
  ExtractfromPack(position,data,lumpedflag_);
  //thermal energy responsible for statistical forces
  ExtractfromPack(position,data,thermalenergy_);
  // gaussrule_
  int gausrule_integer;
  ExtractfromPack(position,data,gausrule_integer);
  gaussrule_ = DRT::UTILS::GaussRule1D(gausrule_integer); //explicit conversion from integer to enum
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                          cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::Beam2::Lines()
{
  lines_.resize(1);
  lines_[0] = this;
  return &lines_[0];
}




//------------- class Beam2Register: -------------------------------------


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2Register::Beam2Register(DRT::Element::ElementType etype):
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2Register::Beam2Register(
                               const DRT::ELEMENTS::Beam2Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2Register* DRT::ELEMENTS::Beam2Register::Clone() const
{
  return new DRT::ELEMENTS::Beam2Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2Register::Pack(vector<char>& data) const
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


/*-----------------------------------------------------------------------*
 |  Unpack data (public)                                      cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2Register::Unpack(const vector<char>& data)
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
 |  dtor (public)                                            cyron 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Beam2Register::~Beam2Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           cyron 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2Register::Print(ostream& os) const
{
  os << "Beam2Register ";
  ElementRegister::Print(os);
  return;
}


int DRT::ELEMENTS::Beam2Register::Initialize(DRT::Discretization& dis)
{	
	 return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM2
