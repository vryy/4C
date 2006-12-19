/*!----------------------------------------------------------------------
\file drt_condition.cpp
\brief A condition of any kind

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "drt_condition.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Condition::Condition(ConditionType type) :
Container(),
type_(type)
{
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Condition::Condition() :
Container(),
type_(condition_none)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Condition::Condition(const DRT::Condition& old) :
Container(old),
type_(old.type_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Condition::~Condition()
{
  return;
}


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const DRT::Condition& cond)
{
  cond.Print(os); 
  return os;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Condition::Print(ostream& os) const
{
  if (Type()==condition_Dirichlet)           os << "Dirichlet boundary condition: ";
  else if (Type()==condition_Neumann) os << "Neumann boundary condition: ";
  Container::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data from this element into vector of length size     (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
const char* DRT::Condition::Pack(int& size) const
{
  const int sizeint    = sizeof(int);
  const int sizetype   = sizeof(enum ConditionType);

  // pack base class Container
  int basesize=0;
  const char* basedata = Container::Pack(basesize);
  
  size = 
  sizeint +      // size itself 
  sizeint +      // type of this instance of ParObject, see top of ParObject.H
  basesize +     // base class Container
  sizetype +     // type_
  0;             // continue to add data here...
  
  
  char* data = new char[size];
  // pack stuff into vector
  int position = 0;

  // size
  AddtoPack(position,data,size);
  // ParObject type
  int type = ParObject_Condition;
  AddtoPack(position,data,type);
  // add base class
  AddtoPack(position,data,basedata,basesize);
  delete [] basedata;
  // add type_
  AddtoPack(position,data,type_);

  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);
  return data;
}


/*----------------------------------------------------------------------*
 |  Unpack data into this element                              (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool DRT::Condition::Unpack(const char* data)
{
  int position=0;
  
  // extract size
  int size = 0;
  ExtractfromPack(position,data,size);
  // ParObject instance type
  int type=0;
  ExtractfromPack(position,data,type);
  if (type != ParObject_Condition) dserror("Wrong instance type in data");
  // extract base class
  int basesize = SizePack(&data[position]);
  Container::Unpack(&data[position]);
  position += basesize;
  // extract type_
  ExtractfromPack(position,data,type_);
    
  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);
  return true;
}









#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
