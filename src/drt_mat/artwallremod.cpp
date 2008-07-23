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
#include "../drt_io/io_gmsh.H"

extern struct _MATERIAL *mat;

//#define PI        (asin(1.0)*2.0)


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

void MAT::ArtWallRemod::Setup(const int numgp, const int eleid)
{
  a1_ = rcp(new vector<vector<double> > (numgp));
  a2_ = rcp(new vector<vector<double> > (numgp));
  int initflag = matdata_->m.artwallremod->init;
  double gamma = matdata_->m.artwallremod->gamma;
  gamma = gamma * PI/180.0;  // convert to radians
  // switch how to setup/initialize fiber directions
  if (initflag==0){
  // fibers aligned in YZ-plane with gamma around Z in global cartesian cosy
    Epetra_SerialDenseMatrix id(3,3);
    // basis is identity
    for (int i=0; i<3; ++i) id(i,i) = 1.0;
    Epetra_SerialDenseMatrix initstress(3,3);
    for (int gp = 0; gp < numgp; ++gp) {
      a1_->at(gp).resize(3);
      a2_->at(gp).resize(3);
      EvaluateFiberVecs(gp,gamma,id);
    }
  } else if (initflag==1){
  // fibers aligned in local element cosy with gamma around circumferential direction
    vector<double> rad(3);
    vector<double> axi(3);
    vector<double> cir(3);
    int ierr=0;
    // read local (cylindrical) cosy-directions at current element
    frdouble_n("RAD",&rad[0],3,&ierr);
    frdouble_n("AXI",&axi[0],3,&ierr);
    frdouble_n("CIR",&cir[0],3,&ierr);
    if (ierr!=1) dserror("Reading of SO_HEX8 element local cosy failed");
    Epetra_SerialDenseMatrix locsys(3,3);
    // basis is local cosy with third vec e3 = circumferential dir and e2 = axial dir
    for (int i=0; i<3; ++i){
      locsys(i,0) = rad[i];
      locsys(i,1) = axi[i];
      locsys(i,2) = cir[i];
    }
    for (int gp = 0; gp < numgp; ++gp) {
      a1_->at(gp).resize(3);
      a2_->at(gp).resize(3);
      EvaluateFiberVecs(gp,gamma,locsys);
    }
  } else if (initflag==2){
    dserror("Random init not yet implemented for ARTWALLREMOD");
  } else dserror("Unknown init for ARTWALLREMOD");

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ArtWallRemod::Initialize(const int numgp, const int eleid) 
{
  srand ( time(NULL) + 5 + eleid*numgp );

  gamma_ = rcp(new vector<double>);
  a1_ = rcp(new vector<vector<double> > (numgp));
  a2_ = rcp(new vector<vector<double> > (numgp));
  lambda_ = rcp(new vector<vector<double> > (numgp));
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
    a1_->at(gp).resize(3);
    a2_->at(gp).resize(3);
    if (matdata_->m.artwallremod->init == 1){
      // random init
      gamma_->push_back(randominit[gp]);
    } else if (matdata_->m.artwallremod->init == 0){
      // pseudo-isotropic init
      gamma_->push_back(45.*PI/180);
    } else dserror("Unknown remodeling initialization");
    for (int i = 0; i < 3; ++i){
      lambda_->at(gp)[i] = 0.0;
    }
    phi_->push_back(id);
    stresses_->push_back(initstress);
    mytime_->push_back(matdata_->m.artwallremod->rembegt);
    //EvaluateFiberVecs(gp);
    
//    cout << "gamma: " << gamma_->at(gp) << endl;
//    cout << "Phi: " << phi_->at(gp) << endl;
//    cout << "a1: " << PrintVec(a1_->at(gp)) << endl;
//    cout << "a2: " << PrintVec(a2_->at(gp)) << endl;
  }

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
                                  Epetra_SerialDenseVector* stress)

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
  const double incJ = pow(I3,-1.0/3.0);  // J^{-2/3}
  
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

  // S = Svol + Siso
  // Svol = J*kappa*(J-1)
  // Isochoric (deviatoric) part via projection PP:Sbar, see Holzapfel p. 230
  // Siso = J^{-2/3}  Dev[Sbar] = J^{-2/3} [Sbar - 1/3 trace(Sbar C) Cinv
  // for this Wiso trace(C Sbar) = trace(mue I C) = mue I1
  const double third = 1./3.;
  const double p = kappa*(J-1);
  for (int i = 0; i < 6; ++i) {
    (*stress)(i) = J*p * Cinv(i)  + incJ* (mue*Id(i) - third*mue*I1*Cinv(i));
  }

  // Elasticity = Cvol + Ciso, via projection see Holzapfel p. 255
  
  // Cvol = J(p + J dp/dJ) Cinv x Cinv  -  2 J p Cinv o Cinv
  // Ciso = 0 + 2/3 J^{-2/3} Sbar:C Psl - 2/3 (Cinv x Siso + Siso x Cinv)
  
  AddtoCmatHolzapfelProduct((*cmat),Cinv,(-2*J*p));  // -2 J p Cinv o Cinv
  
  const double fac = 2*third*incJ*mue*I1;  // 2/3 J^{-2/3} Sbar:C
  // fac Psl = fac (Cinv o Cinv) - fac/3 (Cinv x Cinv)
  //AddtoCmatHolzapfelProduct((*cmat),Cinv,fac);  // fac Cinv o Cinv

  Epetra_SerialDenseMatrix Psl(6,6);        // Psl = Cinv o Cinv - 1/3 Cinv x Cinv
  AddtoCmatHolzapfelProduct(Psl,Cinv,1.0);  // Cinv o Cinv 
  
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      double Siso_i = incJ* (mue*Id(i) - third*mue*I1*Cinv(i));
      double Siso_j = incJ* (mue*Id(j) - third*mue*I1*Cinv(j));
      (*cmat)(i,j) += J*(p+J*kappa) * Cinv(i) * Cinv(j)  // J(p + J dp/dJ) Cinv x Cinv
             + fac * Psl(i,j)                            // fac Cinv o Cinv 
             - fac*third * Cinv(i) * Cinv(j)             // - fac/3 Cinv x Cinv
             - 2*third * Cinv(i) * Siso_j                // -2/3 Cinv x Siso
             - 2*third * Cinv(j) * Siso_i;               // -2/3 Siso x Cinv
      Psl(i,j) += (-third) * Cinv(i) * Cinv(j);    // on the fly complete Psl needed later
    }
  }
  
  // anisotropic part: ***************************************************
  // W_aniso=(k1/(2.0*k2))*(exp(k2*pow((Ibar_{4,6} - 1.0),2)-1.0)); fiber SEF
  
//  cout << "a1: " << PrintVec(a1_->at(gp)) << endl;
//  cout << "a2: " << PrintVec(a2_->at(gp)) << endl;
  
  // structural tensors in voigt notation
  Epetra_SerialDenseVector A1(6);
  Epetra_SerialDenseVector A2(6);
  for (int i = 0; i < 3; ++i) {
    A1(i) = a1_->at(gp)[i]*a1_->at(gp)[i];
    A2(i) = a2_->at(gp)[i]*a2_->at(gp)[i];
  }
  A1(3) = a1_->at(gp)[0]*a1_->at(gp)[1]; A1(4) = a1_->at(gp)[1]*a1_->at(gp)[2]; A1(5) = a1_->at(gp)[0]*a1_->at(gp)[2];
  A2(3) = a2_->at(gp)[0]*a2_->at(gp)[1]; A2(4) = a2_->at(gp)[1]*a2_->at(gp)[2]; A2(5) = a2_->at(gp)[0]*a2_->at(gp)[2];

  // modified (fiber-) invariants Ibar_{4,6} = J_{4,6} = J^{-2/3} I_{4,6}
  // Voigt: trace(AB) =  a11 b11 + 2 a12 b12 + 2 a13 b13 + a22 b22 + 2 a23 b23 + a33 b33
  const double J4 = incJ * ( A1(0)*C(0) + A1(1)*C(1) + A1(2)*C(2) 
                    + 2.*(A1(3)*C(3) + A1(4)*C(4) + A1(5)*C(5))); //J4 = trace(A1:C^dev)
  const double J6 = incJ * ( A2(0)*C(0) + A2(1)*C(1) + A2(2)*C(2) 
                    + 2.*(A2(3)*C(3) + A2(4)*C(4) + A2(5)*C(5))); //J6 = trace(A2:C^dev)
  const double exp1 = exp(k2*(J4-1.)*(J4-1.));
  const double exp2 = exp(k2*(J6-1.)*(J6-1.));
  
  // fibers take compression only
  double fib1_tension = 1.;
  double fib2_tension = 1.;
  if (J4 < 0.0) fib1_tension = 0.;
  if (J6 < 0.0) fib2_tension = 0.;
  
  // PK2 fiber part in splitted formulation, see Holzapfel p. 271 
  Epetra_SerialDenseVector Sfiso(A1); // first compute Sfbar = dWf/dJ4 A1 + dWf/dJ6 A2
  const double fib1 = fib1_tension* 2.*(k1*(J4-1.)*exp1);  // 2 dWf/dJ4
  const double fib2 = fib2_tension* 2.*(k1*(J6-1.)*exp2);  // 2 dWf/dJ6
  Sfiso.Scale(fib1);
  Epetra_SerialDenseVector Stemp(A2);
  Stemp.Scale(fib2);
  Sfiso += Stemp;
  
  const double traceCSfbar =  Sfiso(0)*C(0) + Sfiso(1)*C(1) + Sfiso(2)*C(2) 
                 + 2.*(Sfiso(3)*C(3) + Sfiso(4)*C(4) + Sfiso(5)*C(5)); // trace(Sfbar C)
  // compute Sfiso = J^{-2/3} * (Sfbar - 1/3 trace(Sfbar C) Cinv
  for (int i = 0; i < 6; ++i) {
    Sfiso(i) = incJ * (Sfiso(i) - third*traceCSfbar*Cinv(i));
  }
  (*stress) += Sfiso;
  
  // Elasticity fiber part in splitted fromulation, see Holzapfel p. 255 and 272
  const double delta7bar1 = fib1_tension* 4.*(k1*exp1 + 2.*k1*k2*(J4-1.)*(J4-1.)*exp1); // 4 d^2Wf/dJ4dJ4
  const double delta7bar2 = fib2_tension* 4.*(k1*exp2 + 2.*k1*k2*(J6-1.)*(J6-1.)*exp2); // 4 d^2Wf/dJ6dJ6
  
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      double A1iso_i = incJ*A1(i)-third*J4*Cinv(i);  // A1iso = J^{-2/3} A1 - 1/3 J4 Cinv 
      double A1iso_j = incJ*A1(j)-third*J4*Cinv(j);
      double A2iso_i = incJ*A2(i)-third*J6*Cinv(i);  // A2iso = J^{-2/3} A2 - 1/3 J6 Cinv
      double A2iso_j = incJ*A2(j)-third*J6*Cinv(j);
      (*cmat)(i,j) += delta7bar1 * A1iso_i * A1iso_j  // delta7bar1 A1iso x A1iso
                    + delta7bar2 * A2iso_i * A2iso_j  // delta7bar2 A2iso x A2iso
                    + 2.*third*incJ*traceCSfbar * Psl(i,j)  // 2/3 J^{-2/3} trace(Sfbar C) Psl
                    - 2.*third* (Cinv(i) * Sfiso(j) + Cinv(j) * Sfiso(i)); // -2/3 (Cinv x Sfiso + Sfiso x Cinv)
    }
  }

  return;
}

void MAT::ArtWallRemod::EvaluateFiberVecs(const int gp, const double gamma, const Epetra_SerialDenseMatrix& locsys)
{
  for (int i = 0; i < 3; ++i) {
    // a1 = cos gamma e1 + sin gamma e2 with e1 related to maximal princ stress, e2 2nd largest
    a1_->at(gp)[i] = cos(gamma)*locsys(i,2) + sin(gamma)*locsys(i,1);
    // a2 = cos gamma e1 - sin gamma e2 with e1 related to maximal princ stress, e2 2nd largest
    a2_->at(gp)[i] = cos(gamma)*locsys(i,2) - sin(gamma)*locsys(i,1);
  }
  
  return;
}

std::string MAT::ArtWallRemod::PrintVec(const vector<double> actvec)
{
  std::stringstream out;
  vector<double>::const_iterator i;
  for (i=actvec.begin(); i<actvec.end(); ++i) {
    out << *i << " ";
  }
  return out.str();
}



#endif

