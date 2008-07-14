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
  const double mue = matdata_->m.artwallremod->mue;
  const double kappa = matdata_->m.artwallremod->kappa;
  const double k1 = matdata_->m.artwallremod->k1;
  const double k2 = matdata_->m.artwallremod->k2;
  
  // right Cauchy-Green Tensor  C = 2 * E + I
  // build identity tensor I
  Epetra_SerialDenseVector Id(6);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  Epetra_SerialDenseVector C(*glstrain);
  C.Scale(2.0);
  C += Id;
  
  // invariants
  const double I1 = C(0) + C(1) + C(2);  // 1st invariant, trace
  const double I3 = C(0)*C(1)*C(2)
        + 0.25 * C(3)*C(4)*C(5)
        - 0.25 * C(1)*C(5)*C(5)
        - 0.25 * C(2)*C(3)*C(3)
        - 0.25 * C(0)*C(4)*C(4);    // 3rd invariant, determinant
  const double J = sqrt(I3);
  const double I3invcubroot = pow(I3,-1.0/3.0);
  
  // invert C
  Epetra_SerialDenseVector Cinv(6);

  Cinv(0) = C(1)*C(2) - 0.25*C(4)*C(4);
  Cinv(1) = C(0)*C(2) - 0.25*C(5)*C(5);
  Cinv(2) = C(0)*C(1) - 0.25*C(3)*C(3);
  Cinv(3) = 0.25*C(5)*C(4) - 0.5*C(3)*C(2);
  Cinv(4) = 0.25*C(3)*C(5) - 0.5*C(0)*C(4);
  Cinv(5) = 0.25*C(3)*C(4) - 0.5*C(5)*C(1);
  Cinv.Scale(1.0/I3);

  // isotropic part: NeoHooke  ************************************************
  // NeoHooke with penalty W = W^dev(C) + U(J) 
  // W = 1/2 mue (^I1-3) + 1/2 kappa (J-1)^2

  // Volumetric part of PK2 stress
  Epetra_SerialDenseVector SVol(Cinv);
  SVol.Scale(kappa*(J-1.0)*J);
  *stress+=SVol;

  // Deviatoric elastic part (2 d W^dev/d C)
  Epetra_SerialDenseVector SDev(Cinv);
  SDev.Scale(-1.0/3.0*I1);
  SDev+=Id;
  SDev.Scale(mue*I3invcubroot);  //mue*I3^(-1/3) (Id-1/3*I1*Cinv)
  *stress+=SDev;

  // elasticity matrix
  double scalar1 = 2.0*kappa*J*J - kappa*J;
  double scalar2 = -2.0*kappa*J*J + 2.0*kappa*J;
  double scalar3 = 2.0/3.0*mue*I3invcubroot*I1;
  double scalar4 = 2.0/3.0*mue*I3invcubroot;

  // add scalar2 Cinv o Cinv (see Holzapfel p. 254)
  AddtoCmatHolzapfelProduct((*cmat),Cinv,scalar2);
  for (int i=0; i<6; ++i)
  {
     for (int j=0; j<6; ++j)
     {
       // add volumetric elastic part 2
       (*cmat)(i,j) += scalar1 * Cinv(i) * Cinv(j) // add scalar Cinv x Cinv
       // add visco-elastic deviatoric part 2
           + (-scalar4)*Id(i)*Cinv(j)// add scalar Id x Cinv
           + (-scalar4)*Id(j)*Cinv(i)// add scalar Cinv x Id
           + (scalar3)*Cinv(i)*Cinv(j)/3.0;// add scalar Cinv x Cinv
     }
  }
  // anisotropic part: ***************************************************
  // W_aniso=(k1/(2.0*k2))*(exp(k2*pow((I_{5,6} - 1.0),2)-1.0)); fiber SEF
  
  // fiber directions
  vector<double> a1(3,0.0);
  vector<double> a2(3,0.0);
  
  // structural tensors in voigt-notation
  Epetra_SerialDenseVector A1(6);
  Epetra_SerialDenseVector A2(6);
  for (int i = 0; i < 3; ++i) {
    A1(i) = a1[i]*a1[i];
    A2(i) = a2[i]*a2[i];
  }
  A1(3) = a1[0]*a1[1]; A1(4) = a1[1]*a1[2]; A1(5) = a1[0]*a1[2];
  A2(3) = a2[0]*a2[1]; A2(4) = a2[1]*a2[2]; A2(5) = a2[0]*a2[2];
  
  // modified (fiber-) invariants
  double I4 = 0.0;  // I4 = A1:C^dev
  double I6 = 0.0;  // I6 = A2:C^dev
  for (int i=0; i<6; ++i){
    for (int j = 0; j < 6; ++j) {
      I4 +=  I3invcubroot*C(i) * A1(j);
      I6 +=  I3invcubroot*C(i) * A2(j);
    }
  }
  
  // scalars
  const double expo1 = k2 * (I4 - 1.0) * (I4 - 1.0);
  const double dWdI4 = I3invcubroot * k1 * exp(expo1);
  const double expo2 = k2 * (I6 - 1.0) * (I6 - 1.0);
  const double dWdI6 = I3invcubroot * k1 * exp(expo2);
  
  // add anisotropic PK2 stress
  Epetra_SerialDenseVector Saniso(A1);
  Saniso.Scale(dWdI4);
  *stress += Saniso;
  Saniso = A2;
  Saniso.Scale(dWdI6);
  *stress += Saniso;
  
  // add anistropic elasticity
  const double scalar = I3invcubroot*I3invcubroot * k1 *2*k2;
  for (int i=0; i<6; ++i){
     for (int j=0; j<6; ++j){
       // add volumetric elastic part 2
       (*cmat)(i,j) += scalar * exp(expo1) * Cinv(i) * A1(j)  // add d2WdI4 Cinv x A1
                    +  scalar * exp(expo1) * Cinv(j) * A1(i)  // add d2WdI4 A1 x Cinv
                    +  scalar * exp(expo2) * Cinv(i) * A2(j)  // add d2WdI6 Cinv x A2
                    +  scalar * exp(expo2) * Cinv(j) * A2(i); // add d2WdI6 A2 x Cinv
     }
  }

  return;
}




#endif

