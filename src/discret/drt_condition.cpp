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
#include "drt_element.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Condition::Condition(const int id, const ConditionType type) :
Container(),
id_(id),
type_(type),
comm_(null)
{ 
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Condition::Condition() :
Container(),
id_(-1),
type_(none),
comm_(null)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Condition::Condition(const DRT::Condition& old) :
Container(old),
id_(old.id_),
type_(old.type_),
comm_(old.comm_)
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
  os << "Condition " << Id() << " ";
  if      (Type()==PointDirichlet)   os << "Point Dirichlet boundary condition: ";
  else if (Type()==LineDirichlet)    os << "Line Dirichlet boundary condition: ";
  else if (Type()==SurfaceDirichlet) os << "Surface Dirichlet boundary condition: ";
  else if (Type()==VolumeDirichlet)  os << "Volume Dirichlet boundary condition: ";
  else if (Type()==PointNeumann)     os << "Point Neumann boundary condition: ";
  else if (Type()==LineNeumann)      os << "Line Neumann boundary condition: ";
  else if (Type()==SurfaceNeumann)   os << "Surface Neumann boundary condition: ";
  else if (Type()==VolumeNeumann)    os << "Volume Neumann boundary condition: ";
  Container::Print(os);
  if ((int)geometry_.size())
  {
    cout << endl;
    cout << "Elements of this condition:\n";
    map<int,RefCountPtr<DRT::Element> >::const_iterator curr;
    for (curr=geometry_.begin(); curr!=geometry_.end(); ++curr)
      cout << "      " << *(curr->second) << endl;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Condition::Pack(vector<char>& data) const
{
  data.resize(0);
  
  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class container
  vector<char> basedata;
  Container::Pack(basedata);
  AddVectortoPack(data,basedata);
  // id_
  AddtoPack(data,id_);
  // type_
  AddtoPack(data,type_);
  
  return;
}

#if 0
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
  sizeint  +     // id_
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
  // id_
  AddtoPack(position,data,id_);
  // add type_
  AddtoPack(position,data,type_);

  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);
  return data;
}
#endif

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Condition::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Container
  vector<char> basedata(0);
  ExtractVectorfromPack(position,data,basedata);
  Container::Unpack(basedata);
  // id_
  ExtractfromPack(position,data,id_);
  // type_
  ExtractfromPack(position,data,type_);
  
  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
} 

#if 0
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
  // id_
  ExtractfromPack(position,data,id_);
  // extract type_
  ExtractfromPack(position,data,type_);
    
  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);
  return true;
}
#endif








#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
