/*!----------------------------------------------------------------------
\file contchainnetw.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "Epetra_SerialDenseSolver.h"
#include "contchainnetw.H"
#include "../drt_lib/linalg_serialdensevector.H"


extern struct _MATERIAL *mat;


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         06/08|
 *----------------------------------------------------------------------*/
MAT::ContChainNetw::ContChainNetw()
  : matdata_(NULL)
{
  isinit_=false;
  histstresscurr_=rcp(new vector<Epetra_SerialDenseVector>);
  artstresscurr_=rcp(new vector<Epetra_SerialDenseVector>);
  histstresslast_=rcp(new vector<Epetra_SerialDenseVector>);
  artstresslast_=rcp(new vector<Epetra_SerialDenseVector>);
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          06/08|
 *----------------------------------------------------------------------*/
MAT::ContChainNetw::ContChainNetw(MATERIAL* matdata)
  : matdata_(matdata)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)         06/08|
 *----------------------------------------------------------------------*/
void MAT::ContChainNetw::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matdata
  int matdata = matdata_ - mat;   // pointer difference to reach 0-entry
  AddtoPack(data,matdata);
  int histsize;
  if (!Initialized())
  {
    histsize=0;
  }
  else 
  {
    histsize = histstresslast_->size();
  }
  AddtoPack(data,2*histsize);  // lenght of history vector(s)
  for (int var = 0; var < histsize; ++var) 
  {
    AddtoPack(data,histstresslast_->at(var));
    AddtoPack(data,artstresslast_->at(var));
  }
 
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)         06/08|
 *----------------------------------------------------------------------*/
void MAT::ContChainNetw::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matdata
  int matdata;
  ExtractfromPack(position,data,matdata);
  matdata_ = &mat[matdata];     // unpack pointer to my specific matdata_

  int twicehistsize;
  ExtractfromPack(position,data,twicehistsize);
  
  if (twicehistsize == 0) isinit_=false;
  
  histstresscurr_=rcp(new vector<Epetra_SerialDenseVector>);
  artstresscurr_=rcp(new vector<Epetra_SerialDenseVector>);
  histstresslast_=rcp(new vector<Epetra_SerialDenseVector>);
  artstresslast_=rcp(new vector<Epetra_SerialDenseVector>);
  for (int var=0; var<twicehistsize; var+=2)
  {
    Epetra_SerialDenseVector tmp(6);
    histstresscurr_->push_back(tmp);
    artstresscurr_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    histstresslast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    artstresslast_->push_back(tmp);
  }
  
  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
}


/*----------------------------------------------------------------------*
 |  Return density                                (public)         06/08|
 *----------------------------------------------------------------------*/
double MAT::ContChainNetw::Density()
{
  return matdata_->m.contchainnetw->density;  // density, returned to evaluate mass matrix
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ContChainNetw::Initialize(const int numgp) 
{
  histstresscurr_=rcp(new vector<Epetra_SerialDenseVector>);
  artstresscurr_=rcp(new vector<Epetra_SerialDenseVector>);
  histstresslast_=rcp(new vector<Epetra_SerialDenseVector>);
  artstresslast_=rcp(new vector<Epetra_SerialDenseVector>);
  const Epetra_SerialDenseVector emptyvec(6);
  histstresscurr_->resize(numgp);
  histstresslast_->resize(numgp);
  artstresscurr_->resize(numgp);
  artstresslast_->resize(numgp);
  for (int j=0; j<numgp; ++j) 
  {
    histstresscurr_->at(j) = emptyvec;
    histstresslast_->at(j) = emptyvec;
    artstresscurr_->at(j) = emptyvec;
    artstresslast_->at(j) = emptyvec;
  }
  isinit_=true;
  return ;
  
}

void MAT::ContChainNetw::Update()
{
  histstresslast_=histstresscurr_;
  artstresslast_=artstresscurr_;
  const Epetra_SerialDenseVector emptyvec(6);//6 stresses for 3D
  histstresscurr_=rcp(new vector<Epetra_SerialDenseVector>);
  artstresscurr_=rcp(new vector<Epetra_SerialDenseVector>);
  const int numgp=histstresslast_->size();
  histstresscurr_->resize(numgp);
  artstresscurr_->resize(numgp);
  for (int j=0; j<numgp; ++j) 
  {
    histstresscurr_->at(j) = emptyvec;
    artstresscurr_->at(j) = emptyvec;
  }
}  

/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)         06/08|
 *----------------------------------------------------------------------*

*/

void MAT::ContChainNetw::Evaluate(const Epetra_SerialDenseVector* glstrain,
                                  const int gp,
                                  Teuchos::ParameterList& params,
                                  Epetra_SerialDenseMatrix* cmat,
                                  Epetra_SerialDenseVector* stress)

{
  return;
}


#endif
