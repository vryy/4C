/*!----------------------------------------------------------------------
\file drt_dofset.cpp
\brief A set of degrees of freedom

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "drt_dofset.H"
#include "iostream"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::DofSet::DofSet() :
ParObject()
{
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::DofSet::DofSet(const DRT::DofSet& old) :
ParObject(old),
dofs_(old.dofs_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::DofSet::~DofSet()
{
  return;
}


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const DRT::DofSet& dofset)
{
  dofset.Print(os); 
  return os;
}


/*----------------------------------------------------------------------*
 |  print this  (public)                                     mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::DofSet::Print(ostream& os) const
{
  os << "NumDof " << NumDof() << " Dof(s) ";
  for (int i=0; i<NumDof(); ++i)
    os << Dofs()[i] << " ";
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data from this element into vector of length size     (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
const char* DRT::DofSet::Pack(int& size) const
{
  const int sizeint    = sizeof(int);

  size = 
  sizeint +          // size itself 
  sizeint +          // type of this instance of ParObject, see top of ParObject.H
  SizeVector(dofs_)+ // no. of degrees of freedom and degrees of freedom
  0;                 // continue to add data here...
  
  
  char* data = new char[size];
  // pack stuff into vector
  int position = 0;

  // size
  AddtoPack(position,data,size);
  // ParObject type
  int type = ParObject_DofSet;    // see top of ParObject.H
  AddtoPack(position,data,type);
  // dofs_
  AddVectortoPack(position,data,dofs_);

  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);
  return data;
}


/*----------------------------------------------------------------------*
 |  Unpack data into this element                              (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool DRT::DofSet::Unpack(const char* data)
{
  int position=0;
  
  // extract size
  int size = 0;
  ExtractfromPack(position,data,size);
  // ParObject instance type
  int type=0;
  ExtractfromPack(position,data,type);
  if (type != ParObject_DofSet) dserror("Wrong instance type in data");
  // dofs_
  ExtractVectorfromPack(position,data,dofs_);
    
  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);
  return true;
}


/*----------------------------------------------------------------------*
 |  set no. of degrees of freedom                              (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
void DRT::DofSet::SetNumDof(const int size)
{
  int oldsize = NumDof();
  dofs_.resize(size);
  for (int i=oldsize; i<size; ++i) dofs_[i] = -1;
  return;
}

/*----------------------------------------------------------------------*
 |  set no. of degrees of freedom                              (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
void DRT::DofSet::SetDof(const int* dofs, const int size)
{
  dofs_.resize(size);
  for (int i=0; i<size; ++i) dofs_[i] = dofs[i];
  return;
}






#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
