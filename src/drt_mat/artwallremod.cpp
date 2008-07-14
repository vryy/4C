/*!----------------------------------------------------------------------
\file artwallremod.cpp
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
#include "artwallremod.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/linalg_utils.H"

extern struct _MATERIAL *mat;



/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         06/08|
 *----------------------------------------------------------------------*/
MAT::ArtWallRemod::ArtWallRemod()
  : matdata_(NULL)
{
  isinit_=false;
  gamma_ = rcp(new vector<double>);
  lambda_ = rcp(new vector<vector<double> >);
  phi_ = rcp(new vector<Epetra_SerialDenseMatrix>);
  stresses_ = rcp(new vector<Epetra_SerialDenseMatrix>);
  mytime_ = rcp(new vector<double>);
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          06/08|
 *----------------------------------------------------------------------*/
MAT::ArtWallRemod::ArtWallRemod(MATERIAL* matdata)
  : matdata_(matdata)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)         06/08|
 *----------------------------------------------------------------------*/
void MAT::ArtWallRemod::Pack(vector<char>& data) const
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
    histsize = gamma_->size();
  }
  AddtoPack(data,histsize);  // lenght of history vector(s)
  for (int var = 0; var < histsize; ++var) 
  {
    AddtoPack(data,gamma_->at(var));
    AddtoPack(data,phi_->at(var));
    AddtoPack(data,stresses_->at(var));
    AddtoPack(data,mytime_->at(var));
  }
  
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)         06/08|
 *----------------------------------------------------------------------*/
void MAT::ArtWallRemod::Unpack(const vector<char>& data)
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

  // history data
  isinit_ = true;
  int histsize;
  ExtractfromPack(position,data,histsize);
  
  if (histsize == 0) isinit_=false;
  gamma_ = rcp(new vector<double>);
  phi_ = rcp(new vector<Epetra_SerialDenseMatrix>);
  stresses_ = rcp(new vector<Epetra_SerialDenseMatrix>);
  mytime_ = rcp(new vector<double>);
  for (int var = 0; var < histsize; ++var) {
    double gamma;
    Epetra_SerialDenseMatrix tmp(3,3);
    ExtractfromPack(position,data,gamma);
    ExtractfromPack(position,data,tmp);
    gamma_->push_back(gamma);
    phi_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    stresses_->push_back(tmp);
    double mytime;
    ExtractfromPack(position,data,mytime);
    mytime_->push_back(mytime);
  }

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ArtWallRemod::Initialize(const int numgp, const int eleid) 
{
  srand ( time(NULL) + 5 + eleid*numgp );

  gamma_ = rcp(new vector<double>);
  lambda_ = rcp(new vector<vector<double> >);
  phi_ = rcp(new vector<Epetra_SerialDenseMatrix>);
  stresses_ = rcp(new vector<Epetra_SerialDenseMatrix>);
  mytime_ = rcp(new vector<double>);
  // initial basis is identity
  Epetra_SerialDenseMatrix id(3,3);
  for (int i=0; i<3; ++i) id(i,i) = 1.0;
  Epetra_SerialDenseMatrix initstress(3,3);
  
  vector<double> randominit(8);
  for (int i = 0; i < numgp; ++i) {
    randominit[i] = ( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
  }
  
  // initialize remodelling parameters
  for(int gp=0; gp<numgp; ++gp){
    lambda_->at(gp).resize(3);
    if (matdata_->m.artwallremod->initran == 1){
      // random init
      gamma_->push_back(randominit[gp]);
    } else if (matdata_->m.artwallremod->initran == 0){
      // pseudo-isotropic init
      gamma_->push_back(1.0);
    } else dserror("Unknown remodeling initialization");
    for (int i = 0; i < 3; ++i){
      lambda_->at(gp)[i] = 0.0;
    }
    phi_->push_back(id);
    stresses_->push_back(initstress);
    mytime_->push_back(matdata_->m.artwallremod->rembegt);
  }

  //mytime_ = 1.0;  // carefull!
  isinit_ = true;
  
  return ;
  
}


/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)         06/08|
 *----------------------------------------------------------------------*

*/

void MAT::ArtWallRemod::Evaluate(const Epetra_SerialDenseVector* glstrain,
                                  const int gp,
                                  Teuchos::ParameterList& params,
                                  Epetra_SerialDenseMatrix* cmat,
                                  Epetra_SerialDenseVector* stress,
                                  int eleId)

{
  return;
}




#endif

