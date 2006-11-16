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
#ifdef TRILINOS_PACKAGE

#include "shell8.H"
#include "dserror.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Shell8::Shell8(int id) :
CCADISCRETIZATION::Element(id,element_shell8),
forcetype_(s8_none),
thickness_(0.0),
ngptri_(0),
nhyb_(0),
ans_(0),
sdc_(1.0),
material_(0)
{
  ngp_[0] = ngp_[1] = ngp_[2] = 0;
  eas_[0] = eas_[1] = eas_[2] = eas_[3] = eas_[4] = 0;
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Shell8::Shell8(const CCADISCRETIZATION::Shell8& old) :
CCADISCRETIZATION::Element(old),
forcetype_(old.forcetype_),
thickness_(old.thickness_),
ngptri_(old.ngptri_),
nhyb_(old.nhyb_),
ans_(old.ans_),
sdc_(old.sdc_),
material_(old.material_)
{
  for (int i=0; i<3; ++i) ngp_[i] = old.ngp_[i];
  for (int i=0; i<5; ++i) eas_[i] = old.eas_[i];
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Shell8 and return pointer to it (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Element* CCADISCRETIZATION::Shell8::Clone() const
{
  CCADISCRETIZATION::Shell8* newelement = new CCADISCRETIZATION::Shell8(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data from this element into vector of length size     (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
const char* CCADISCRETIZATION::Shell8::Pack(int& size) const
{
  const int sizeint    = sizeof(int);
  const int sizeforcetype = sizeof(enum ForceType);
  const int sizedouble = sizeof(double);
  //const int sizechar   = sizeof(char);

  // pack base class
  int         basesize = 0;
  const char* basedata = Element::Pack(basesize);
   
  // calculate size of vector data
  size = 
  sizeint +        // size itself
  basesize +       // base class data
  sizeforcetype +  // forcetype_
  sizedouble +     // thickness_
  3*sizeint +      // ngp_
  sizeint +        // ngptri_  
  sizeint +        // nhyb_
  5*sizeint +      // eas_
  sizeint +        // ans_
  sizedouble +     // sdc_
  sizeint +        // material_
  0;            // continue to add data here...
  
  char* data = new char[size];
  
  // pack stuff into vector
  int position = 0;
  
  // add size
  AddtoPack(position,data,size);
  // add base class
  AddtoPack(position,data,basedata,basesize);
  delete [] basedata;
  // forcetype_
  AddtoPack(position,data,&forcetype_,sizeof(enum ForceType));
  // thickness_
  AddtoPack(position,data,thickness_);
  // ngp_
  AddtoPack(position,data,ngp_,3*sizeint);
  // ngptri_
  AddtoPack(position,data,ngptri_);
  // nhyb_
  AddtoPack(position,data,nhyb_);
  // eas_
  AddtoPack(position,data,eas_,5*sizeint);
  // ans_
  AddtoPack(position,data,ans_);
  // sdc_
  AddtoPack(position,data,sdc_);
  // material_
  AddtoPack(position,data,material_);
  // continue to add stuff here  
    
  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);

  return data;
}

/*----------------------------------------------------------------------*
 |  Unpack data into this element                              (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool CCADISCRETIZATION::Shell8::Unpack(const char* data)
{
  const int sizeint    = sizeof(int);
  //const int sizeforcetype = sizeof(enum ForceType);
  //const int sizedouble = sizeof(double);
  //const int sizechar   = sizeof(char);

  int position = 0;

  // extract size
  int size = 0;
  ExtractfromPack(position,data,size);
  
  // extract base class
  int basesize = SizePack(&data[position]);
  Element::Unpack(&data[position]);
  position += basesize;
  
  // forcetype_
  ExtractfromPack(position,data,forcetype_);
  // thickness_
  ExtractfromPack(position,data,thickness_);
  // ngp_
  ExtractfromPack(position,data,ngp_,3*sizeint);
  // ngptri_
  ExtractfromPack(position,data,ngptri_);
  // nhyb_
  ExtractfromPack(position,data,nhyb_);
  // eas_
  ExtractfromPack(position,data,eas_,5*sizeint);
  // ans_
  ExtractfromPack(position,data,ans_);
  // sdc_
  ExtractfromPack(position,data,sdc_);
  // material_
  ExtractfromPack(position,data,material_);
  // extract more stuff here

  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);

  return true;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Shell8::~Shell8()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Shell8::Print(ostream& os) const
{
  os << "Shell8 ";
  Element::Print(os);
  return;
}
















#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SHELL8
