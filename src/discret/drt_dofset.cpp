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
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::DofSet::Pack(vector<char>& data) const
{
  data.resize(0);
  // type of parobject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // dofs_
  AddtoPack(data,dofs_);
  
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::DofSet::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // dofs_
  ExtractfromPack(position,data,dofs_);
  
  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
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
