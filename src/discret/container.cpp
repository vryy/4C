/*!----------------------------------------------------------------------
\file container.cpp
\brief A data storage container

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "container.H"
#include "dserror.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Container::Container() :
ParObject()
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Container::Container(const CCADISCRETIZATION::Container& old) :
ParObject(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Container::~Container()
{
  return;
}


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CCADISCRETIZATION::Container& cont)
{
  cont.Print(os); 
  return os;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Container::Print(ostream& os) const
{
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data from this element into vector of length size     (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
const char* CCADISCRETIZATION::Container::Pack(int& size) const
{
  return NULL;
}


/*----------------------------------------------------------------------*
 |  Unpack data into this element                              (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool CCADISCRETIZATION::Container::Unpack(const char* data)
{
  return true;
}

/*----------------------------------------------------------------------*
 |  Add stuff to the container                                 (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Container::Add(const string& name, const int* data, const int num)
{
  // check whether this name was used before
  map<string,RefCountPtr<vector<int> > >::iterator curr = intdata_.find(name);
  if (curr != intdata_.end()) dserror("Record with name %s already exists",&name[0]);
  
  // get data in a vector
  RefCountPtr<vector<int> > storage = rcp(new vector<int>(num));
  vector<int>& access = *storage;
  for (int i=0; i<num; ++i) access[i] = data[i];
  
  // store the vector
  intdata_[name] = storage;
  return;
}

/*----------------------------------------------------------------------*
 |  Add stuff to the container                                 (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Container::Add(const string& name, const double* data, const int num)
{
  // check whether this name was used before
  map<string,RefCountPtr<vector<double> > >::iterator curr = doubledata_.find(name);
  if (curr != doubledata_.end()) dserror("Record with name %s already exists",&name[0]);

  // get data in a vector
  RefCountPtr<vector<double> > storage = rcp(new vector<double>(num));
  vector<double>& access = *storage;
  for (int i=0; i<num; ++i) access[i] = data[i];

  // store the vector
  doubledata_[name] = storage;
  return;
}










#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
