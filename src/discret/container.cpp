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
DRT::Container::Container() :
ParObject()
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Container::Container(const DRT::Container& old) :
ParObject(old),
intdata_(old.intdata_),
doubledata_(old.doubledata_),
stringdata_(old.stringdata_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Container::~Container()
{
  return;
}


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const DRT::Container& cont)
{
  cont.Print(os); 
  return os;
}


/*----------------------------------------------------------------------*
 |  Pack data from this element into vector of length size     (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
const char* DRT::Container::Pack(int& size) const
{
  const int sizeint    = sizeof(int);

  size = 
  sizeint +      // size itself 
  sizeint +      // type of this instance of ParObject, see top of ParObject.H
  sizeint +      // number of entries in intdata_
  sizeint +      // number of entries in doubledata_
  sizeint +      // number of entries in stringdata_
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
  map<string,string>::const_iterator scurr;
  for (scurr = stringdata_.begin(); scurr != stringdata_.end(); ++scurr)
  {
    size += SizeString(scurr->first);      // size of the key
    size += SizeString(scurr->second);     // size of data
  }


  char* data = new char[size];
  // pack stuff into vector
  int position = 0;

  // size
  AddtoPack(position,data,size);
  // ParObject type
  int type = ParObject_Container;
  AddtoPack(position,data,type);
  // intdata_.size()
  int tmp = (int)intdata_.size();
  AddtoPack(position,data,tmp);
  // doubledata_.size()
  tmp = (int)doubledata_.size();
  AddtoPack(position,data,tmp);
  // stringdata_.size()
  tmp = (int)stringdata_.size();
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
  for (scurr = stringdata_.begin(); scurr != stringdata_.end(); ++scurr)
  {
    AddStringtoPack(position,data,scurr->first);
    AddStringtoPack(position,data,scurr->second);
  }
    
  
  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);

  return data;
}


/*----------------------------------------------------------------------*
 |  Unpack data into this element                              (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool DRT::Container::Unpack(const char* data)
{
  int position=0;
  
  // extract size
  int size = 0;
  ExtractfromPack(position,data,size);
  // ParObject instance type
  int type=0;
  ExtractfromPack(position,data,type);
  if (type != ParObject_Container) dserror("Wrong instance type in data");
  // extract number of entries in intdata_
  int numintdata=0;
  ExtractfromPack(position,data,numintdata);
  // extract number of entries in doubledata_
  int numdoubledata=0;
  ExtractfromPack(position,data,numdoubledata);
  // extract number of entries in stringdata_
  int numstringdata=0;
  ExtractfromPack(position,data,numstringdata);
  
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
  for (int i=0; i<numstringdata; ++i)
  {
    string tmp;
    ExtractStringfromPack(position,data,tmp);
    string tmp2;
    ExtractStringfromPack(position,data,tmp2);
    Add(tmp,tmp2);
  }  
    
  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);
  return true;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Container::Print(ostream& os) const
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
  map<string,string>::const_iterator scurr;
  for (scurr = stringdata_.begin(); scurr != stringdata_.end(); ++scurr)
    os << scurr->first << " : " << scurr->second << endl;
  return;
}

/*----------------------------------------------------------------------*
 |  Add stuff to the container                                 (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
void DRT::Container::Add(const string& name, const int* data, const int num)
{
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
void DRT::Container::Add(const string& name, const double* data, const int num)
{
  // get data in a vector
  RefCountPtr<vector<double> > storage = rcp(new vector<double>(num));
  vector<double>& access = *storage;
  for (int i=0; i<num; ++i) access[i] = data[i];

  // store the vector
  doubledata_[name] = storage;
  return;
}

/*----------------------------------------------------------------------*
 |  Add stuff to the container                                 (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
void DRT::Container::Add(const string& name, const string& data)
{
  // store the data
  stringdata_[name] = data;
  return;
}

/*----------------------------------------------------------------------*
 |  Delete stuff from the container                            (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
void DRT::Container::Delete(const string& name)
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

  map<string,string>::iterator scurr = stringdata_.find(name);
  if (scurr != stringdata_.end()) 
  {
    stringdata_.erase(name);
    return;
  }
  return;
}


/*----------------------------------------------------------------------*
 |  Get a string                                               (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
const string* DRT::Container::GetString(const string& name) const
{
  map<string,string>::const_iterator curr = stringdata_.find(name);
  if (curr != stringdata_.end())
    return &(curr->second);
  else return NULL;
  return NULL;
}







#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
