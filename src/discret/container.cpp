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
ParObject(old),
intdata_(old.intdata_),
doubledata_(old.doubledata_)
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
 |  Pack data from this element into vector of length size     (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
const char* CCADISCRETIZATION::Container::Pack(int& size) const
{
  const int sizeint    = sizeof(int);

  size = 
  sizeint +      // size itself 
  sizeint +      // number of entries in intdata_
  sizeint +      // number of entries in doubledata_
  0;             // continue to add data here...
  
  
  // iterate through data and record sizes
  map<string,RefCountPtr<vector<int> > >::const_iterator curr;
  for (curr = intdata_.begin(); curr != intdata_.end(); ++curr)
  {
    size += SizeString(curr->first);     // size of the key
    size += SizeVector(*(curr->second)); // size of data
  }
  map<string,RefCountPtr<vector<double> > >::const_iterator dcurr;
  for (dcurr = doubledata_.begin(); dcurr != doubledata_.end(); ++dcurr)
  {
    size += SizeString(dcurr->first);      // size of the key
    size += SizeVector(*(dcurr->second));  // size of data
  }


  char* data = new char[size];
  // pack stuff into vector
  int position = 0;

  // size
  AddtoPack(position,data,size);
  // intdata_.size()
  int tmp = (int)intdata_.size();
  AddtoPack(position,data,tmp);
  // doubledata_.size()
  tmp = (int)doubledata_.size();
  AddtoPack(position,data,tmp);
  
  // continue to pack data here...
  
  for (curr = intdata_.begin(); curr != intdata_.end(); ++curr)
  {
    AddStringtoPack(position,data,curr->first);
    AddVectortoPack(position,data,*(curr->second));
  }
  for (dcurr = doubledata_.begin(); dcurr != doubledata_.end(); ++dcurr)
  {
    AddStringtoPack(position,data,dcurr->first);
    AddVectortoPack(position,data,*(dcurr->second));
  }
    
  
  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);
  return data;
}


/*----------------------------------------------------------------------*
 |  Unpack data into this element                              (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool CCADISCRETIZATION::Container::Unpack(const char* data)
{
  int position=0;
  
  // extract size
  int size = 0;
  ExtractfromPack(position,data,size);
  // extract number of entries in intdata_
  int numintdata=0;
  ExtractfromPack(position,data,numintdata);
  // extract number of entries in doubledata_
  int numdoubledata=0;
  ExtractfromPack(position,data,numdoubledata);
  
  // iterate and extract
  for (int i=0; i<numintdata; ++i)
  {
    string tmp;
    ExtractStringfromPack(position,data,tmp);
    vector<int> tmpvector;
    ExtractVectorfromPack(position,data,tmpvector);
    // add to the hits container
    Add(tmp,tmpvector);
  }  
  for (int i=0; i<numdoubledata; ++i)
  {
    string tmp;
    ExtractStringfromPack(position,data,tmp);
    vector<double> tmpvector;
    ExtractVectorfromPack(position,data,tmpvector);
    Add(tmp,tmpvector);
  }  
    
  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);
  return true;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Container::Print(ostream& os) const
{
  map<string,RefCountPtr<vector<int> > >::const_iterator curr;
  for (curr = intdata_.begin(); curr != intdata_.end(); ++curr)
  {
    vector<int>& data = *(curr->second);
    os << curr->first << " : ";
    for (int i=0; i<(int)data.size(); ++i) os << data[i] << " ";
    os << endl;
  }
  map<string,RefCountPtr<vector<double> > >::const_iterator dcurr;
  for (dcurr = doubledata_.begin(); dcurr != doubledata_.end(); ++dcurr)
  {
    vector<double>& data = *(dcurr->second);
    os << dcurr->first << " : ";
    for (int i=0; i<(int)data.size(); ++i) os << data[i] << " ";
    os << endl;
  }
  return;
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

/*----------------------------------------------------------------------*
 |  Delete stuff from the container                            (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
void CCADISCRETIZATION::Container::Delete(const string& name)
{
  map<string,RefCountPtr<vector<int> > >::iterator icurr = intdata_.find(name);
  if (icurr != intdata_.end()) 
  {
    intdata_.erase(name);
    return;
  }

  map<string,RefCountPtr<vector<double> > >::iterator dcurr = doubledata_.find(name);
  if (dcurr != doubledata_.end()) 
  {
    doubledata_.erase(name);
    return;
  }
  return;
}









#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
