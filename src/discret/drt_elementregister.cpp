/*!----------------------------------------------------------------------
\file drt_elementregister.cpp
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

#include "drt_elementregister.H"
#include "drt_dserror.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::ElementRegister::ElementRegister(DRT::Element::ElementType etype) :
ParObject(),
etype_(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::ElementRegister::ElementRegister(const DRT::ElementRegister& old) :
ParObject(old),
etype_(old.etype_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::ElementRegister::~ElementRegister()
{
  return;
}


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const DRT::ElementRegister& eler)
{
  eler.Print(os); 
  return os;
}


/*----------------------------------------------------------------------*
 |  print element (public)                                   mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::ElementRegister::Print(ostream& os) const
{
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ElementRegister::Pack(vector<char>& data) const
{
  data.resize(0);
  
  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add type of element
  AddtoPack(data,etype_);
  
  return;
}

#if 0
/*----------------------------------------------------------------------*
 |  Pack data from this element into vector of length size     (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
const char* DRT::ElementRegister::Pack(int& size) const
{
  const int sizeint    = sizeof(int);
  const int sizetype   = sizeof(enum DRT::Element::ElementType);

  // calculate size of vector
  size = 
  sizeint                +   // holds size itself
  sizeint                +   // type of this instance of ParObject, see top of ParObject.H
  sizetype               +   // holds type of element
  0;                         // continue to add data here...
  
  char* data = new char[size];
  int position = 0;

  // add size
  AddtoPack(position,data,size);
  // ParObject type
  int type = ParObject_ElementRegister;
  AddtoPack(position,data,type);
  // add type of element
  AddtoPack(position,data,etype_);

  return data;
}
#endif

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ElementRegister::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  ExtractfromPack(position,data,etype_);
  
  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
} 

#if 0
/*----------------------------------------------------------------------*
 |  Unpack data into this element                              (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool DRT::ElementRegister::Unpack(const char* data)
{
  int position = 0;
  
  // extract size
  int size = 0;
  ExtractfromPack(position,data,size);
  // ParObject instance type
  int type=0;
  ExtractfromPack(position,data,type);
  if (type != ParObject_ElementRegister) dserror("Wrong instance type in data");
  // extract type
  ExtractfromPack(position,data,etype_);

  return true;
}
#endif









#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
