/*!----------------------------------------------------------------------
\file drt_parobject.cpp
\brief

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

#include "drt_parobject.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::ParObject::ParObject()
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::ParObject::ParObject(const DRT::ParObject& old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::ParObject::~ParObject()
{
  return;
}

void DRT::ParObject::AddtoPack(PackBuffer& data, const ParObject& obj)
{
  obj.Pack( data );
}

void DRT::ParObject::AddtoPack(PackBuffer& data, const ParObject* obj)
{
  obj->Pack( data );
}

/*----------------------------------------------------------------------*
 |      a vector<int> specialization                           (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::AddtoPack(PackBuffer& data, const std::vector<int>& stuff)
{
  int numele = stuff.size();
  AddtoPack(data,numele);
  AddtoPack(data,&stuff[0],numele*sizeof(int));
  return;
}
/*----------------------------------------------------------------------*
 |      a vector<double> specialization                        (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::AddtoPack(PackBuffer& data, const std::vector<double>& stuff)
{
  int numele = stuff.size();
  AddtoPack(data,numele);
  AddtoPack(data,&stuff[0],numele*sizeof(double));
  return;
}
/*----------------------------------------------------------------------*
 |        a vector<char> specialization                        (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::AddtoPack(PackBuffer& data, const std::vector<char>& stuff)
{
  int numele = stuff.size();
  AddtoPack(data,numele);
  AddtoPack(data,&stuff[0],numele*sizeof(char));
  return;
}

/*----------------------------------------------------------------------*
 |        a map <int,double> specialization                    (public) |
 |                                                           mgit 09/10 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::AddtoPack(PackBuffer& data, const std::map<int,double> & stuff)
{

  int numentries = (int) stuff.size();
  AddtoPack(data,numentries);

  // iterator
  std::map<int,double>::const_iterator colcurr;

  int i=0;
  for(colcurr=stuff.begin();colcurr!=stuff.end();++colcurr)
  {
    AddtoPack(data,colcurr->first);
    AddtoPack(data,colcurr->second);
    ++i;
  }

  if(i!=numentries)
    dserror("Something wrong with number of elements");

  return;
}

/*----------------------------------------------------------------------*
 |        a map <int,int>    specialization                    (public) |
 |                                                         schott 03/12 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::AddtoPack(PackBuffer& data, const std::map<int,int> & stuff)
{

  int numentries = (int) stuff.size();
  AddtoPack(data,numentries);

  // iterator
  std::map<int,int>::const_iterator colcurr;

  int i=0;
  for(colcurr=stuff.begin();colcurr!=stuff.end();++colcurr)
  {
    AddtoPack(data,colcurr->first);
    AddtoPack(data,colcurr->second);
    ++i;
  }

  if(i!=numentries)
    dserror("Something wrong with number of elements");

  return;
}

/*----------------------------------------------------------------------*
 |        a map <string,int>    specialization                 (public) |
 |                                                       nagler 07/2012 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::AddtoPack(PackBuffer& data, const std::map<std::string,int> & stuff)
{

  int numentries = (int) stuff.size();
  AddtoPack(data,numentries);

  // iterator
  std::map<std::string,int>::const_iterator colcurr;

  int i=0;
  for(colcurr=stuff.begin();colcurr!=stuff.end();++colcurr)
  {
    AddtoPack(data,colcurr->first);
    AddtoPack(data,colcurr->second);
    ++i;
  }

  if(i!=numentries)
    dserror("Something wrong with number of elements");

  return;
}

/*----------------------------------------------------------------------*
 |        a set<int> specialization                            (public) |
 |                                                           mgit 09/10 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::AddtoPack(PackBuffer& data, const std::set<int> & stuff)
{

  int numentries = (int) stuff.size();
  AddtoPack(data,numentries);

  // iterator
  std::set<int>::const_iterator colcurr;

  int i=0;
  for(colcurr=stuff.begin();colcurr!=stuff.end();++colcurr)
  {
    AddtoPack(data,*colcurr);
    ++i;
  }

  if(i!=numentries)
    dserror("Something wrong with number of elements");

   return;
}

/*----------------------------------------------------------------------*
 | a Epetra_SerialDenseMatrix specialization                   (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::AddtoPack(PackBuffer& data, const Epetra_SerialDenseMatrix& stuff)
{
  int m = stuff.M();
  int n = stuff.N();
  AddtoPack(data,m);
  AddtoPack(data,n);
  double* A = stuff.A();
  AddtoPack(data,A,n*m*sizeof(double));
  return;
}
/*----------------------------------------------------------------------*
 | a Epetra_SerialDenseVector specialization                   (public) |
 |                                                     TK & MAF  05/08  |
 *----------------------------------------------------------------------*/
void DRT::ParObject::AddtoPack(PackBuffer& data, const Epetra_SerialDenseVector& stuff)
{
  int m = stuff.Length();
  AddtoPack(data,m);
  double* A = stuff.Values();
  AddtoPack(data,A,m*sizeof(double));
  return;
}
/*----------------------------------------------------------------------*
 | a LINALG::SerialDenseMatrix specialization                  (public) |
 |                                                          henke 12/09 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::AddtoPack(PackBuffer& data, const LINALG::SerialDenseMatrix& stuff)
{
  int m = stuff.M();
  int n = stuff.N();
  AddtoPack(data,m);
  AddtoPack(data,n);
  double* A = stuff.A();
  AddtoPack(data,A,n*m*sizeof(double));
  return;
}
/*----------------------------------------------------------------------*
 | a string specialization                                     (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::AddtoPack(PackBuffer& data, const string& stuff)
{
  int numele = stuff.size();
  AddtoPack(data,numele);
  AddtoPack(data,&stuff[0],numele*sizeof(char));
  return;
}
/*----------------------------------------------------------------------*
 | a int vector specialization                                 (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::ExtractfromPack(std::vector<char>::size_type& position, const std::vector<char>& data, std::vector<int>& stuff)
{
  int dim = 0;
  ExtractfromPack(position,data,dim);
  stuff.resize(dim);
  int size = dim*sizeof(int);
  ExtractfromPack(position,data,&stuff[0],size);
  return;
}
/*----------------------------------------------------------------------*
 | a double vector specialization                                 (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::ExtractfromPack(std::vector<char>::size_type& position, const std::vector<char>& data, std::vector<double>& stuff)
{
  int dim = 0;
  ExtractfromPack(position,data,dim);
  stuff.resize(dim);
  int size = dim*sizeof(double);
  ExtractfromPack(position,data,&stuff[0],size);
  return;
}
/*----------------------------------------------------------------------*
 | a char vector specialization                                 (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::ExtractfromPack(std::vector<char>::size_type& position, const std::vector<char>& data, std::vector<char>& stuff)
{
  int dim = 0;
  ExtractfromPack(position,data,dim);
  stuff.resize(dim);
  int size = dim*sizeof(char);
  ExtractfromPack(position,data,&stuff[0],size);
  return;
}

/*----------------------------------------------------------------------*
 | a map specialization                                        (public) |
 |                                                         schott 03/12 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::ExtractfromPack(std::vector<char>::size_type& position, const std::vector<char>& data, std::map<int,int>& stuff)
{

  int numentries = 0;
  ExtractfromPack(position,data,numentries);

  stuff.clear();

  for(int i=0;i<numentries;i++)
  {
    int dof;
    int value;
    ExtractfromPack(position,data,dof);
    ExtractfromPack(position,data,value);

    //add to map
    stuff.insert (std::pair<int,int>(dof,value));

  }

  return;
}

/*----------------------------------------------------------------------*
 | a map specialization                                        (public) |
 |                                                           mgit 09/10 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::ExtractfromPack(std::vector<char>::size_type& position, const std::vector<char>& data, std::map<int,double>& stuff)
{

  int numentries = 0;
  ExtractfromPack(position,data,numentries);

  stuff.clear();

  for(int i=0;i<numentries;i++)
  {
    int dof;
    double value;
    ExtractfromPack(position,data,dof);
    ExtractfromPack(position,data,value);

    //add to map
    stuff.insert (std::pair<int,double>(dof,value));

  }

  return;
}

/*----------------------------------------------------------------------*
 | a map specialization                                        (public) |
 |                                                         nagler 07/12 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::ExtractfromPack(std::vector<char>::size_type& position, const std::vector<char>& data, std::map<std::string,int>& stuff)
{

  int numentries = 0;
  ExtractfromPack(position,data,numentries);

  stuff.clear();

  for(int i=0;i<numentries;i++)
  {
  	std::string keys;
    int values;
    ExtractfromPack(position,data,keys);
    ExtractfromPack(position,data,values);

    //add to map
    stuff.insert(std::pair<std::string,int>(keys,values));

  }

  return;
}

/*----------------------------------------------------------------------*
 | a set specialization                                        (public) |
 |                                                           mgit 09/10 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::ExtractfromPack(std::vector<char>::size_type& position, const std::vector<char>& data, std::set<int>& stuff)
{

  int numentries = 0;
  ExtractfromPack(position,data,numentries);

  stuff.clear();

  for(int i=0;i<numentries;i++)
  {
    int value;
    ExtractfromPack(position,data,value);

    //add to map
    stuff.insert(value);
  }

  return;
}

/*----------------------------------------------------------------------*
 | a Epetra_SerialDenseMatrix specialization                   (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::ExtractfromPack(std::vector<char>::size_type& position, const std::vector<char>& data,
                                     Epetra_SerialDenseMatrix& stuff)
{
  int m = 0;
  ExtractfromPack(position,data,m);
  int n = 0;
  ExtractfromPack(position,data,n);
  stuff.Reshape(m,n);
  double* a = stuff.A();
  if(m*n)
    ExtractfromPack(position,data,a,n*m*sizeof(double));
  return;
}
/*----------------------------------------------------------------------*
 | a Epetra_SerialDenseVector specialization                   (public) |
 |                                                     TK & MAF  05/08  |
 *----------------------------------------------------------------------*/
void DRT::ParObject::ExtractfromPack(std::vector<char>::size_type& position, const std::vector<char>& data,
                                     Epetra_SerialDenseVector& stuff)
{
  int m = 0;
  ExtractfromPack(position,data,m);
  stuff.Resize(m);
  double* a = stuff.Values();
  if(m)
    ExtractfromPack(position,data,a,m*sizeof(double));
  return;
}

/*----------------------------------------------------------------------*
 | a LINALG::SerialDenseMatrix specialization                  (public) |
 |                                                          henke 12/09 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::ExtractfromPack(std::vector<char>::size_type& position, const std::vector<char>& data,
                                     LINALG::SerialDenseMatrix& stuff)
{
  int m = 0;
  ExtractfromPack(position,data,m);
  int n = 0;
  ExtractfromPack(position,data,n);
  stuff.Reshape(m,n);
  double* a = stuff.A();
  if(m*n)
    ExtractfromPack(position,data,a,n*m*sizeof(double));
  return;
}

/*----------------------------------------------------------------------*
 | a string specialization                                     (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ParObject::ExtractfromPack(std::vector<char>::size_type& position, const std::vector<char>& data, string& stuff)
{
  int dim = 0;
  ExtractfromPack(position,data,dim);
  stuff.resize(dim);
  int size = dim*sizeof(char);
  ExtractfromPack(position,data,&stuff[0],size);
  return;
}



