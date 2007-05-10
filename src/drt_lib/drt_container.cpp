/*!----------------------------------------------------------------------
\file drt_container.cpp

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

#include "drt_container.H"
#include "drt_dserror.H"



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
stringdata_(old.stringdata_)
{
  // we want a true deep copy of the data, not a reference
  map<string,RefCountPtr<vector<int> > >::const_iterator ifool;
  for (ifool=old.intdata_.begin(); ifool!=old.intdata_.end(); ++ifool)
    intdata_[ifool->first] = rcp(new vector<int>(*(ifool->second)));

  map<string,RefCountPtr<vector<double> > >::const_iterator dfool;
  for (dfool=old.doubledata_.begin(); dfool!=old.doubledata_.end(); ++dfool)
    doubledata_[dfool->first] = rcp(new vector<double>(*(dfool->second)));

  map<string,RefCountPtr<Epetra_SerialDenseMatrix> >::const_iterator curr;
  for (curr=old.matdata_.begin();curr!=old.matdata_.end(); ++curr)
    matdata_[curr->first] = rcp(new Epetra_SerialDenseMatrix(*(curr->second)));

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
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Container::Pack(vector<char>& data) const
{
  data.resize(0);
  // no. of objects in maps
  const int indatasize = intdata_.size();
  const int doubledatasize = doubledata_.size();
  const int stringdatasize = stringdata_.size();
  const int matdatasize = matdata_.size();
  
  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // indatasize
  AddtoPack(data,indatasize);
  // doubledatasize
  AddtoPack(data,doubledatasize);
  // stringdatasize
  AddtoPack(data,stringdatasize);
  // matdatasize
  AddtoPack(data,matdatasize);
  // iterate through intdata_ and add to pack
  map<string,RefCountPtr<vector<int> > >::const_iterator icurr;
  for (icurr = intdata_.begin(); icurr != intdata_.end(); ++icurr)
  {
    AddtoPack(data,icurr->first);
    AddtoPack(data,*(icurr->second));
  }
  // iterate though doubledata_ and add to pack
  map<string,RefCountPtr<vector<double> > >::const_iterator dcurr;
  for (dcurr = doubledata_.begin(); dcurr != doubledata_.end(); ++dcurr)
  {
    AddtoPack(data,dcurr->first);
    AddtoPack(data,*(dcurr->second));
  }
  // iterate through stringdata_ and add to pack
  map<string,string>::const_iterator scurr;
  for (scurr = stringdata_.begin(); scurr != stringdata_.end(); ++scurr)
  {
    AddtoPack(data,scurr->first);
    AddtoPack(data,scurr->second);
  }
  // iterate though matdata_ and add to pack
  map<string,RefCountPtr<Epetra_SerialDenseMatrix> >::const_iterator mcurr;
  for (mcurr=matdata_.begin(); mcurr!=matdata_.end(); ++mcurr)
  {
    AddtoPack(data,mcurr->first);
    AddtoPack(data,*(mcurr->second));
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Container::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract no. objects in intdata_
  int intdatasize = 0;
  ExtractfromPack(position,data,intdatasize);
  // extract no. objects in doubledata_
  int doubledatasize = 0;
  ExtractfromPack(position,data,doubledatasize);
  // extract no. objects in stringdata_
  int stringdatasize = 0;
  ExtractfromPack(position,data,stringdatasize);
  // extract no. objects in matdata_
  int matdatasize = 0;
  ExtractfromPack(position,data,matdatasize);
  
  // iterate though records of intdata_ and extract
  for (int i=0; i<intdatasize; ++i)
  {
    string key;
    ExtractfromPack(position,data,key);
    vector<int> value(0);
    ExtractfromPack(position,data,value);
    Add(key,value);
  }
  
  // iterate though records of doubledata_ and extract
  for (int i=0; i<doubledatasize; ++i)
  {
    string key;
    ExtractfromPack(position,data,key);
    vector<double> value(0);
    ExtractfromPack(position,data,value);
    Add(key,value);
  }
  
  // iterate though records of stringdata_ and extract
  for (int i=0; i<stringdatasize; ++i)
  {
    string key;
    ExtractfromPack(position,data,key);
    string value;
    ExtractfromPack(position,data,value);
    Add(key,value);
  }

  // iterate though records of matdata_ and extract
  for (int i=0; i<matdatasize; ++i)
  {
    string key;
    ExtractfromPack(position,data,key);
    Epetra_SerialDenseMatrix value(0,0);
    ExtractfromPack(position,data,value);
    Add(key,value);
  }
  
  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
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
    //os << endl;
  }

  map<string,RefCountPtr<vector<double> > >::const_iterator dcurr;
  for (dcurr = doubledata_.begin(); dcurr != doubledata_.end(); ++dcurr)
  {
    vector<double>& data = *(dcurr->second);
    os << dcurr->first << " : ";
    for (int i=0; i<(int)data.size(); ++i) os << data[i] << " ";
    //os << endl;
  }

  map<string,string>::const_iterator scurr;
  for (scurr = stringdata_.begin(); scurr != stringdata_.end(); ++scurr)
    os << scurr->first << " : " << scurr->second << " ";
    
  map<string,RefCountPtr<Epetra_SerialDenseMatrix> >::const_iterator matcurr;
  for (matcurr=matdata_.begin(); matcurr!=matdata_.end(); ++matcurr)
    os << endl << matcurr->first << " :\n" << *(matcurr->second);

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
 |  Add stuff to the container                                 (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
void DRT::Container::Add(const string& name, const Epetra_SerialDenseMatrix& matrix)
{
  matdata_[name] = rcp(new Epetra_SerialDenseMatrix(matrix));
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
  
  map<string,RefCountPtr<Epetra_SerialDenseMatrix> >::iterator matcurr = matdata_.find(name);
  if (matcurr != matdata_.end()) 
  {
    matdata_.erase(name);
    return;
  }
  
  
  return;
}


// specializations have to be in the same namespace as the template
// compilers can be pedantic!
namespace DRT
{
/*----------------------------------------------------------------------*
 |  Get a vector<int> specialization                           (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
  template<> vector<int>* Container::Get(const string& name)
  {
    map<string,RefCountPtr<vector<int> > >::iterator icurr = intdata_.find(name);
    if (icurr != intdata_.end())
      return icurr->second.get();
    else return NULL;
  }
/*----------------------------------------------------------------------*
 |  Get a vector<double> specialization                        (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
  template<> vector<double>* Container::Get(const string& name)
  {
    map<string,RefCountPtr<vector<double> > >::iterator dcurr = doubledata_.find(name);
    if (dcurr != doubledata_.end())
      return dcurr->second.get();
    else return NULL;
  }
/*----------------------------------------------------------------------*
 |  Get a string specialization                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
  template<> string* Container::Get(const string& name)
  {
    map<string,string>::iterator scurr = stringdata_.find(name);
    if (scurr != stringdata_.end())
      return &(scurr->second);
    else return NULL;
  }
/*----------------------------------------------------------------------*
 |  Get a Epetra_SerialDensMatrix specialization               (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
  template<> Epetra_SerialDenseMatrix* Container::Get(const string& name)
  {
    map<string,RefCountPtr<Epetra_SerialDenseMatrix> >::iterator mcurr = matdata_.find(name);
    if (mcurr != matdata_.end())
      return mcurr->second.get();
    else return NULL;
  }
} // end of namespace DRT



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
