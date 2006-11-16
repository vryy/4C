/*!----------------------------------------------------------------------
\file designelement.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "designelement.H"
#include "designdiscretization.H"
#include "dserror.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::DesignElement::DesignElement(int id, enum ElementType type) :
Element(id,type)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::DesignElement::DesignElement(const CCADISCRETIZATION::DesignElement& old) :
Element(old),
lentityid_(old.lentityid_),
lorientation_(old.lorientation_),
lentity_(old.lentity_),
hentityid_(old.hentityid_),
horientation_(old.horientation_),
hentity_(old.hentity_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::DesignElement* CCADISCRETIZATION::DesignElement::Clone() const
{
  CCADISCRETIZATION::DesignElement* newelement = new CCADISCRETIZATION::DesignElement(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data from this element into vector of length size     (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
const char* CCADISCRETIZATION::DesignElement::Pack(int& size) const
{
  const int sizeint    = sizeof(int);
  //const int sizedouble = sizeof(double);
  //const int sizechar   = sizeof(char);

  // pack base class
  int         basesize = 0;
  const char* basedata = Element::Pack(basesize);
   
  // calculate size of vector data
  size = 
  sizeint +                     // size itself
  basesize +                    // base class data
  SizeVector(lentityid_) +      // lentityid_
  SizeVector(lorientation_) +   // lorientation_
  SizeVector(hentityid_) +      // hentityid_
  SizeVector(horientation_) +   // horientation_
  0;                            // continue to add data here...
  
  char* data = new char[size];
  
  // pack stuff into vector
  int position = 0;
  
  // add size
  AddtoPack(position,data,size);
  // add base class
  AddtoPack(position,data,basedata,basesize);
  delete [] basedata;
  // lentityid_
  AddVectortoPack(position,data,lentityid_);
  // lorientation_
  AddVectortoPack(position,data,lorientation_);
  // hentityid_
  AddVectortoPack(position,data,hentityid_);
  // horientation_
  AddVectortoPack(position,data,horientation_);
  // continue to add stuff here  
    
  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);

  return data;
}

/*----------------------------------------------------------------------*
 |  Unpack data into this element                              (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool CCADISCRETIZATION::DesignElement::Unpack(const char* data)
{
  //const int sizeint    = sizeof(int);
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
  // extract lentityid_
  ExtractVectorfromPack(position,data,lentityid_);
  // extract lorientation_
  ExtractVectorfromPack(position,data,lorientation_);
  // extract hentityid_
  ExtractVectorfromPack(position,data,hentityid_);
  // extract horientation_
  ExtractVectorfromPack(position,data,horientation_);
  // extract more stuff here

  // ptrs to lower and higher entities are not valid anymore
  lentity_.resize(0);
  hentity_.resize(0);

  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);

  return true;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::DesignElement::~DesignElement()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::DesignElement::Print(ostream& os) const
{
  os << "DesignElement ";
  Element::Print(os);
  if (NumHigherEntityIds())
  {
    os << " HigherEntities ";
    for (int i=0; i<NumHigherEntityIds(); ++i) 
      os << HigherEntityIds()[i] << " ";
  }
  if (NumLowerEntityIds())
  {
    os << " LowerEntities ";
    for (int i=0; i<NumLowerEntityIds(); ++i)
      os << LowerEntityIds()[i] << " ";
  }
  return;
}


/*----------------------------------------------------------------------*
 | set lower entity ids (public)                             mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::DesignElement::SetLowerEntities(const int nele, const int* ids, const int* orientation)
{
  lentityid_.resize(nele);
  lorientation_.resize(nele);
  for (int i=0; i<nele; ++i)
  {
    lentityid_[i]      = ids[i];
    lorientation_[i] = orientation[i];
  }
  return;
}

/*----------------------------------------------------------------------*
 | set higher entity ids (public)                            mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::DesignElement::SetHigherEntities(const int nele, const int* ids, const int* orientation)
{
  hentityid_.resize(nele);
  horientation_.resize(nele);
  for (int i=0; i<nele; ++i)
  {
    hentityid_[i]      = ids[i];
    horientation_[i] = orientation[i];
  }
  return;
}


/*----------------------------------------------------------------------*
 | get pointers to lower entities (public)                   mwgee 11/06|
 *----------------------------------------------------------------------*/
bool CCADISCRETIZATION::DesignElement::BuildLowerElementPointers(const Discretization& lower)
{
  const int nlowerid = NumLowerEntityIds();
  if (!nlowerid) return true;
  lentity_.resize(nlowerid);
  const int* lowerid = LowerEntityIds();
  for (int i=0; i<nlowerid; ++i)
  {
    Element* ele = lower.gElement(lowerid[i]);
    if (!ele) return false;
    lentity_[i] = ele;
  }
  return true;
}












#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
