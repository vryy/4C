/*!----------------------------------------------------------------------
\file drt_container.cpp

\brief A data storage container

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_container.H"
#include "drt_dserror.H"
#include "Epetra_Vector.h"


DRT::ContainerType DRT::ContainerType::instance_;


DRT::ParObject* DRT::ContainerType::Create( const std::vector<char> & data )
{
  DRT::Container* object = new DRT::Container();
  object->Unpack(data);
  return object;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Container::Container() :
ParObject(), data_( NULL )
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Container::Container(const DRT::Container& old) :
ParObject(old), data_( NULL )
{
  if ( old.data_ != NULL )
  {
    data_ = new Teuchos::ParameterList();

    // we want a true deep copy of the data, not a reference

    for ( Teuchos::ParameterList::ConstIterator i=old.data_->begin();
          i!=old.data_->end();
          ++i )
    {
      const std::string & n = old.data_->name( i );
      const ParameterEntry & e = old.data_->entry( i );

      if ( e.isType<std::string>() )
      {
        std::string & v = e.getValue<std::string>( NULL );
        data_->set( n, v );
      }
      else if ( e.isType<Teuchos::RCP<std::vector<int> > >() )
      {
        Teuchos::RCP<std::vector<int> > & v = e.getValue<Teuchos::RCP<std::vector<int> > >( NULL );
        data_->set( n, Teuchos::rcp( new std::vector<int>( *v ) ) );
      }
      else if ( e.isType<Teuchos::RCP<std::vector<double> > >() )
      {
        Teuchos::RCP<std::vector<double> > & v = e.getValue<Teuchos::RCP<std::vector<double> > >( NULL );
        data_->set( n, Teuchos::rcp( new std::vector<double>( *v ) ) );
      }
      else if ( e.isType<Teuchos::RCP<Epetra_SerialDenseMatrix> >() )
      {
        Teuchos::RCP<Epetra_SerialDenseMatrix> & v = e.getValue<Teuchos::RCP<Epetra_SerialDenseMatrix> >( NULL );
        data_->set( n, Teuchos::rcp( new Epetra_SerialDenseMatrix( *v ) ) );
      }
      else if ( e.isType<Teuchos::RCP<Epetra_MultiVector> >() )
      {
        Teuchos::RCP<Epetra_MultiVector> & v = e.getValue<Teuchos::RCP<Epetra_MultiVector> >( NULL );
        data_->set( n, Teuchos::rcp( new Epetra_MultiVector( *v ) ) );
      }
      else
      {
        dserror( "unsuppored entity type" );
      }
    }
  }
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Container::~Container()
{
  Clear();
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
void DRT::Container::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  map<string,RCP<vector<int> > >             intdata;
  map<string,RCP<vector<double> > >          doubledata;
  map<string,string>                         stringdata;
  map<string,RCP<Epetra_SerialDenseMatrix> > matdata;
  map<string,RCP<Epetra_MultiVector> >       evecdata;

  if ( data_ != NULL )
  {
    // we want a true deep copy of the data, not a reference

    for ( Teuchos::ParameterList::ConstIterator i=data_->begin();
          i!=data_->end();
          ++i )
    {
      const std::string & n = data_->name( i );
      const ParameterEntry & e = data_->entry( i );

      if ( e.isType<std::string>() )
      {
        std::string & v = e.getValue<std::string>( NULL );
        stringdata[n] = v;
      }
      else if ( e.isType<Teuchos::RCP<std::vector<int> > >() )
      {
        Teuchos::RCP<std::vector<int> > & v = e.getValue<Teuchos::RCP<std::vector<int> > >( NULL );
        intdata[n] = v;
      }
      else if ( e.isType<Teuchos::RCP<std::vector<double> > >() )
      {
        Teuchos::RCP<std::vector<double> > & v = e.getValue<Teuchos::RCP<std::vector<double> > >( NULL );
        doubledata[n] = v;
      }
      else if ( e.isType<Teuchos::RCP<Epetra_SerialDenseMatrix> >() )
      {
        Teuchos::RCP<Epetra_SerialDenseMatrix> & v = e.getValue<Teuchos::RCP<Epetra_SerialDenseMatrix> >( NULL );
        matdata[n] = v;
      }
      else if ( e.isType<Teuchos::RCP<Epetra_MultiVector> >() )
      {
        Teuchos::RCP<Epetra_MultiVector> & v = e.getValue<Teuchos::RCP<Epetra_MultiVector> >( NULL );
        evecdata[n] = v;
      }
      else
      {
        dserror( "unsuppored entity type" );
      }
    }
  }

  // no. of objects in maps
  const int indatasize = intdata.size();
  const int doubledatasize = doubledata.size();
  const int stringdatasize = stringdata.size();
  const int matdatasize = matdata.size();

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
  map<string,RCP<vector<int> > >::const_iterator icurr;
  for (icurr = intdata.begin(); icurr != intdata.end(); ++icurr)
  {
    AddtoPack(data,icurr->first);
    AddtoPack(data,*(icurr->second));
  }
  // iterate though doubledata_ and add to pack
  map<string,RCP<vector<double> > >::const_iterator dcurr;
  for (dcurr = doubledata.begin(); dcurr != doubledata.end(); ++dcurr)
  {
    AddtoPack(data,dcurr->first);
    AddtoPack(data,*(dcurr->second));
  }
  // iterate through stringdata_ and add to pack
  map<string,string>::const_iterator scurr;
  for (scurr = stringdata.begin(); scurr != stringdata.end(); ++scurr)
  {
    AddtoPack(data,scurr->first);
    AddtoPack(data,scurr->second);
  }
  // iterate though matdata_ and add to pack
  map<string,RCP<Epetra_SerialDenseMatrix> >::const_iterator mcurr;
  for (mcurr=matdata.begin(); mcurr!=matdata.end(); ++mcurr)
  {
    AddtoPack(data,mcurr->first);
    AddtoPack(data,*(mcurr->second));
  }

  // on purpose the map<string,RCP<Epetra_MultiVector> > evecdata_
  // is NOT included in the Pack/Unpack stuff

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Container::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
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

  // on purpose the map<string,RCP<Epetra_MultiVector> > evecdata_
  // is NOT included in the Pack/Unpack stuff

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Container::Print(ostream& os) const
{
  map<string,RCP<vector<int> > >             intdata;
  map<string,RCP<vector<double> > >          doubledata;
  map<string,string>                         stringdata;
  map<string,RCP<Epetra_SerialDenseMatrix> > matdata;
  map<string,RCP<Epetra_MultiVector> >       evecdata;

  if ( data_ != NULL )
  {
    // we want a true deep copy of the data, not a reference

    for ( Teuchos::ParameterList::ConstIterator i=data_->begin();
          i!=data_->end();
          ++i )
    {
      const std::string & n = data_->name( i );
      const ParameterEntry & e = data_->entry( i );

      if ( e.isType<std::string>() )
      {
        std::string & v = e.getValue<std::string>( NULL );
        stringdata[n] = v;
      }
      else if ( e.isType<Teuchos::RCP<std::vector<int> > >() )
      {
        Teuchos::RCP<std::vector<int> > & v = e.getValue<Teuchos::RCP<std::vector<int> > >( NULL );
        intdata[n] = v;
      }
      else if ( e.isType<Teuchos::RCP<std::vector<double> > >() )
      {
        Teuchos::RCP<std::vector<double> > & v = e.getValue<Teuchos::RCP<std::vector<double> > >( NULL );
        doubledata[n] = v;
      }
      else if ( e.isType<Teuchos::RCP<Epetra_SerialDenseMatrix> >() )
      {
        Teuchos::RCP<Epetra_SerialDenseMatrix> & v = e.getValue<Teuchos::RCP<Epetra_SerialDenseMatrix> >( NULL );
        matdata[n] = v;
      }
      else if ( e.isType<Teuchos::RCP<Epetra_MultiVector> >() )
      {
        Teuchos::RCP<Epetra_MultiVector> & v = e.getValue<Teuchos::RCP<Epetra_MultiVector> >( NULL );
        evecdata[n] = v;
      }
      else
      {
        dserror( "unsuppored entity type" );
      }
    }
  }

  map<string,RCP<vector<int> > >::const_iterator curr;
  for (curr = intdata.begin(); curr != intdata.end(); ++curr)
  {
    vector<int>& data = *(curr->second);
    os << curr->first << " : ";
    for (int i=0; i<(int)data.size(); ++i) os << data[i] << " ";
    //os << endl;
  }

  map<string,RCP<vector<double> > >::const_iterator dcurr;
  for (dcurr = doubledata.begin(); dcurr != doubledata.end(); ++dcurr)
  {
    vector<double>& data = *(dcurr->second);
    os << dcurr->first << " : ";
    for (int i=0; i<(int)data.size(); ++i) os << data[i] << " ";
    //os << endl;
  }

  map<string,string>::const_iterator scurr;
  for (scurr = stringdata.begin(); scurr != stringdata.end(); ++scurr)
    os << scurr->first << " : " << scurr->second << " ";

  map<string,RCP<Epetra_SerialDenseMatrix> >::const_iterator matcurr;
  for (matcurr=matdata.begin(); matcurr!=matdata.end(); ++matcurr)
    os << endl << matcurr->first << " :\n" << *(matcurr->second);

  map<string,RCP<Epetra_MultiVector> >::const_iterator eveccurr;
  for (eveccurr=evecdata.begin(); eveccurr!=evecdata.end(); ++eveccurr)
    os << endl << eveccurr->first << " (type Epetra_Vector or Epetra_MultiVector) \n";

  return;
}

/*----------------------------------------------------------------------*
 |  Add stuff to the container                                 (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
void DRT::Container::Add(const string& name, const int* data, const int num)
{
  // get data in a vector
  RCP<vector<int> > storage = rcp(new vector<int>(data, data+num));
  Add( name, storage );
}

/*----------------------------------------------------------------------*
 |  Add stuff to the container                                 (public) |
 |                                                             lw 04/08 |
 *----------------------------------------------------------------------*/
void DRT::Container::Add(const string& name, RCP<vector<int> > data)
{
  if ( data_==NULL )
    data_ = new Teuchos::ParameterList();
  data_->set( name, data );
}

/*----------------------------------------------------------------------*
 |  Add stuff to the container                                 (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
void DRT::Container::Add(const string& name, const double* data, const int num)
{
  // get data in a vector
  RCP<vector<double> > storage = rcp(new vector<double>(data, data+num));
  Add( name, storage );
}

/*----------------------------------------------------------------------*
 |  Add stuff to the container                                 (public) |
 |                                                             lw 04/08 |
 *----------------------------------------------------------------------*/
void DRT::Container::Add(const string& name, RCP<vector<double> > data)
{
  if ( data_==NULL )
    data_ = new Teuchos::ParameterList();
  data_->set( name, data );
}

/*----------------------------------------------------------------------*
 |  Add stuff to the container                                 (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
void DRT::Container::Add(const string& name, const string& data)
{
  if ( data_==NULL )
    data_ = new Teuchos::ParameterList();
  data_->set( name, data );
}

/*----------------------------------------------------------------------*
 |  Add stuff to the container                                 (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
void DRT::Container::Add(const string& name, const Epetra_SerialDenseMatrix& matrix)
{
  if ( data_==NULL )
    data_ = new Teuchos::ParameterList();
  data_->set( name, rcp(new Epetra_SerialDenseMatrix(matrix)) );
}

/*----------------------------------------------------------------------*
 |  Add stuff to the container                                 (public) |
 |                                                             lw 04/08 |
 *----------------------------------------------------------------------*/
void DRT::Container::Add(const string& name, RCP<Epetra_SerialDenseMatrix> matrix)
{
  if ( data_==NULL )
    data_ = new Teuchos::ParameterList();
  data_->set( name, matrix );
}

/*----------------------------------------------------------------------*
 |  Add stuff to the container                                 (public) |
 |                                                            gee 03/10 |
 *----------------------------------------------------------------------*/
void DRT::Container::Add(const string& name, Epetra_MultiVector& data)
{
  if ( data_==NULL )
    data_ = new Teuchos::ParameterList();
  data_->set( name, rcp(new Epetra_MultiVector(data)) );
}

/*----------------------------------------------------------------------*
 |  Add stuff to the container                                 (public) |
 |                                                            gee 03/10 |
 *----------------------------------------------------------------------*/
void DRT::Container::Add(const string& name, Epetra_Vector& data)
{
  if ( data_==NULL )
    data_ = new Teuchos::ParameterList();
  data_->set( name, rcp(new Epetra_Vector(data)) );
}

/*----------------------------------------------------------------------*
 |  Delete stuff from the container                            (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
void DRT::Container::Delete(const string& name)
{
  if ( data_!=NULL )
  {
    data_->remove( name );
  }
}


// note that there exist no general implementation of this template,
// there are specializations only! This means that nothing else works other
// than explicitly implemented here
namespace DRT
{
/*----------------------------------------------------------------------*
 |  Get a vector<int> specialization                           (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
  template<> const vector<int>* Container::Get(const string& name) const
  {
    if ( data_!=NULL )
    {
      const ParameterEntry * e = data_->getEntryPtr( name );
      if ( e != NULL and e->isType<Teuchos::RCP<std::vector<int> > >() )
      {
        return &*e->getValue<Teuchos::RCP<std::vector<int> > >( NULL );
      }
    }
    return NULL;
  }
/*----------------------------------------------------------------------*
 |  Get a vector<double> specialization                        (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
  template<> const vector<double>* Container::Get(const string& name) const
  {
    if ( data_!=NULL )
    {
      const ParameterEntry * e = data_->getEntryPtr( name );
      if ( e != NULL and e->isType<Teuchos::RCP<std::vector<double> > >() )
      {
        return &*e->getValue<Teuchos::RCP<std::vector<double> > >( NULL );
      }
    }
    return NULL;
  }
/*----------------------------------------------------------------------*
 |  Get a string specialization                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
  template<> const string* Container::Get(const string& name) const
  {
    if ( data_!=NULL )
    {
      const ParameterEntry * e = data_->getEntryPtr( name );
      if ( e != NULL and e->isType<std::string>() )
      {
        return & e->getValue<std::string>( NULL );
      }
    }
    return NULL;
  }
/*----------------------------------------------------------------------*
 |  Get a Epetra_SerialDensMatrix specialization               (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
  template<> const Epetra_SerialDenseMatrix* Container::Get(const string& name) const
  {
    if ( data_!=NULL )
    {
      const ParameterEntry * e = data_->getEntryPtr( name );
      if ( e != NULL and e->isType<Teuchos::RCP<Epetra_SerialDenseMatrix> >() )
      {
        return &*e->getValue<Teuchos::RCP<Epetra_SerialDenseMatrix> >( NULL );
      }
    }
    return NULL;
  }
/*----------------------------------------------------------------------*
 |  Get a Epetra_MultiVector specialization                    (public) |
 |                                                            gee 03/10 |
 *----------------------------------------------------------------------*/
  template<> const Epetra_MultiVector* Container::Get(const string& name) const
  {
    if ( data_!=NULL )
    {
      const ParameterEntry * e = data_->getEntryPtr( name );
      if ( e != NULL and e->isType<Teuchos::RCP<Epetra_MultiVector> >() )
      {
        return &*e->getValue<Teuchos::RCP<Epetra_MultiVector> >( NULL );
      }
    }
    return NULL;
  }
/*----------------------------------------------------------------------*
 |  Get a Epetra_MultiVector specialization                    (public) |
 |                                                            gee 03/10 |
 *----------------------------------------------------------------------*/
  template<> const Epetra_Vector* Container::Get(const string& name) const
  {
    if ( data_!=NULL )
    {
      const ParameterEntry * e = data_->getEntryPtr( name );
      if ( e != NULL and e->isType<Teuchos::RCP<Epetra_MultiVector> >() )
      {
        Epetra_MultiVector * mv = &*e->getValue<Teuchos::RCP<Epetra_MultiVector> >( NULL );
        Epetra_Vector * ev = dynamic_cast<Epetra_Vector*>( mv );
        if ( ev==NULL )
        {
          dserror("Object in container is NOT Epetra_Vector");
        }
        return ev;
      }
    }
    return NULL;
  }
/*----------------------------------------------------------------------*
 |  Get a vector<int> specialization                           (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
  template<> vector<int>* Container::GetMutable(const string& name)
  {
    if ( data_!=NULL )
    {
      ParameterEntry * e = data_->getEntryPtr( name );
      if ( e != NULL and e->isType<Teuchos::RCP<std::vector<int> > >() )
      {
        return &*e->getValue<Teuchos::RCP<std::vector<int> > >( NULL );
      }
    }
    return NULL;
  }
/*----------------------------------------------------------------------*
 |  Get a vector<double> specialization                        (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
  template<> vector<double>* Container::GetMutable(const string& name)
  {
    if ( data_!=NULL )
    {
      ParameterEntry * e = data_->getEntryPtr( name );
      if ( e != NULL and e->isType<Teuchos::RCP<std::vector<double> > >() )
      {
        return &*e->getValue<Teuchos::RCP<std::vector<double> > >( NULL );
      }
    }
    return NULL;
  }
/*----------------------------------------------------------------------*
 |  Get a string specialization                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
  template<> string* Container::GetMutable(const string& name)
  {
    if ( data_!=NULL )
    {
      ParameterEntry * e = data_->getEntryPtr( name );
      if ( e != NULL and e->isType<std::string>() )
      {
        return & e->getValue<std::string>( NULL );
      }
    }
    return NULL;
  }
/*----------------------------------------------------------------------*
 |  Get a Epetra_SerialDensMatrix specialization               (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
  template<> Epetra_SerialDenseMatrix* Container::GetMutable(const string& name)
  {
    if ( data_!=NULL )
    {
      ParameterEntry * e = data_->getEntryPtr( name );
      if ( e != NULL and e->isType<Teuchos::RCP<Epetra_SerialDenseMatrix> >() )
      {
        return &*e->getValue<Teuchos::RCP<Epetra_SerialDenseMatrix> >( NULL );
      }
    }
    return NULL;
  }
/*----------------------------------------------------------------------*
 |  Get a Epetra_MultiVector specialization                    (public) |
 |                                                            gee 03/10 |
 *----------------------------------------------------------------------*/
  template<> Epetra_MultiVector* Container::GetMutable(const string& name)
  {
    if ( data_!=NULL )
    {
      ParameterEntry * e = data_->getEntryPtr( name );
      if ( e != NULL and e->isType<Teuchos::RCP<Epetra_MultiVector> >() )
      {
        return &*e->getValue<Teuchos::RCP<Epetra_MultiVector> >( NULL );
      }
    }
    return NULL;
  }
/*----------------------------------------------------------------------*
 |  Get a Epetra_MultiVector specialization                    (public) |
 |                                                            gee 03/10 |
 *----------------------------------------------------------------------*/
  template<> Epetra_Vector* Container::GetMutable(const string& name)
  {
    if ( data_!=NULL )
    {
      ParameterEntry * e = data_->getEntryPtr( name );
      if ( e != NULL and e->isType<Teuchos::RCP<Epetra_MultiVector> >() )
      {
        Epetra_MultiVector * mv = &*e->getValue<Teuchos::RCP<Epetra_MultiVector> >( NULL );
        Epetra_Vector * ev = dynamic_cast<Epetra_Vector*>( mv );
        if ( ev==NULL )
        {
          dserror("Object in container is NOT Epetra_Vector");
        }
        return ev;
      }
    }
    return NULL;
  }
} // end of namespace DRT

/*----------------------------------------------------------------------*
 |  just get an int back                                       (public) |
 |                                                          chfoe 11/07 |
 *----------------------------------------------------------------------*/
int DRT::Container::GetInt(const string& name) const
{
  const vector<int>* vecptr = Get<vector<int> >(name);
  if(vecptr==NULL) dserror("Integer %s cannot be read from the container.",name.c_str());
  if( vecptr->size()!=1 ) dserror("Trying to read integer from vector of wrong length.");
  return (*vecptr)[0];
}

/*----------------------------------------------------------------------*
 |  just get a double back                                     (public) |
 |                                                             lw 12/07 |
 *----------------------------------------------------------------------*/
double DRT::Container::GetDouble(const string& name) const
{
  const vector<double>* vecptr = Get<vector<double> >(name);
  if(vecptr==NULL) dserror("Double %s cannot be read from the container.",name.c_str());
  if( vecptr->size()!=1 ) dserror("Trying to read double from vector of wrong length.");
  return (*vecptr)[0];
}




#endif  // #ifdef CCADISCRET
