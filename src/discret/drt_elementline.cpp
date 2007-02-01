/*!----------------------------------------------------------------------
\file drt_elementsurface.cpp
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

#include "drt_elementline.H"
#include "drt_dserror.H"




/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 01/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ElementLine::ElementLine(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes) :
DRT::Element(id,element_elementline,owner)
{
  SetNodeIds(nnode,nodeids);
  node_.resize(nnode);
  for (int i=0; i<nnode; ++i) node_[i] = nodes[i];
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ElementLine::ElementLine(const DRT::ElementLine& old) :
DRT::Element(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ElementLine::Clone() const
{
  DRT::ElementLine* newelement = new DRT::ElementLine(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data from this element into vector of length size     (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
const char* DRT::ElementLine::Pack(int& size) const
{
  const int sizeint    = sizeof(int);
  //const int sizedouble = sizeof(double);
  //const int sizechar   = sizeof(char);

  // pack base class
  int         basesize = 0;
  const char* basedata = Element::Pack(basesize);
  
  // calculate size of vector data
  size = 
  sizeint +                    // size itself
  sizeint +                    // type of this instance of ParObject, see top of ParObject.H
  basesize +                   // base class data
  0;            // continue to add data here...
  
  char* data = new char[size];
  
  // pack stuff into vector
  int position = 0;
  
  // add size
  AddtoPack(position,data,size);
  // ParObject type
  int type = UniqueParObjectId();
  AddtoPack(position,data,type);
  // add base class
  AddtoPack(position,data,basedata,basesize);
  delete [] basedata;
  // continue to add stuff here  
    
  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);

  return data;
}

/*----------------------------------------------------------------------*
 |  Unpack data into this element                              (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
bool DRT::ElementLine::Unpack(const char* data)
{
  //const int sizeint    = sizeof(int);
  //const int sizeforcetype = sizeof(enum ForceType);
  //const int sizedouble = sizeof(double);
  //const int sizechar   = sizeof(char);

  int position = 0;

  // extract size
  int size = 0;
  ExtractfromPack(position,data,size);
  // ParObject instance type
  int type=0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("Wrong instance type in data");
  // extract base class
  int basesize = SizePack(&data[position]);
  Element::Unpack(&data[position]);
  position += basesize;
  // extract more stuff here

  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);

  return true;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ElementLine::~ElementLine()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::ElementLine::Print(ostream& os) const
{
  os << "ElementLine ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return ElementSurfaceRegister (public)      mwgee 01/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ElementLine::ElementRegister() const
{
  return null;
}




#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
