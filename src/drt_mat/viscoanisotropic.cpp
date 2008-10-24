/*!----------------------------------------------------------------------
\file viscoanisotropic.cpp
\brief

<pre>
Maintainer: Moritz Frenzel & Thomas Kloeppel
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
#include "viscoanisotropic.H"
#include "../drt_lib/linalg_serialdensevector.H"


extern struct _MATERIAL *mat;  ///< C-style material struct


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         05/08|
 *----------------------------------------------------------------------*/
MAT::ViscoAnisotropic::ViscoAnisotropic()
  : matdata_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          05/08|
 *----------------------------------------------------------------------*/
MAT::ViscoAnisotropic::ViscoAnisotropic(MATERIAL* matdata)
  : matdata_(matdata)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoAnisotropic::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matdata
  int matdata = matdata_ - mat;   // pointer difference to reach 0-entry
  AddtoPack(data,matdata);

  int numgp;
  int numhist;
  if (!Initialized())
  {
    numgp = 0;
    numhist = 0;
  }
  else
  {
    numgp = a1_->size();   // size is number of gausspoints
    numhist = histstresscurr_->size();
  }
  AddtoPack(data,numgp);
  // Pack internal variables
  for (int gp = 0; gp < numgp; ++gp)
  {
    AddtoPack(data,a1_->at(gp));
    AddtoPack(data,a2_->at(gp));
  }
  //Pack history data
  if (numhist != 0) AddtoPack(data,numhist); // Length of history vector(s)
  for (int var = 0; var < numhist; ++var)
  {
    AddtoPack(data,histstresslast_->at(var));
    AddtoPack(data,artstresslast_->at(var));
  }
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoAnisotropic::Unpack(const vector<char>& data)
{
  isinit_=true;
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matdata
  int matdata;
  ExtractfromPack(position,data,matdata);
  matdata_ = &mat[matdata];     // unpack pointer to my specific matdata_
  int numgp, numhist;
  ExtractfromPack(position,data,numgp);
  if (numgp == 0){ // no history data to unpack
    isinit_=false;
    if (position != (int)data.size())
      dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
    return;
  }
  // unpack fiber internal variables
  a1_ = rcp(new vector<vector<double> >(numgp));
  a2_ = rcp(new vector<vector<double> >(numgp));
  for (int gp = 0; gp < numgp; ++gp) {
    vector<double> a;
    ExtractfromPack(position,data,a);
    a1_->at(gp) = a;
    ExtractfromPack(position,data,a);
    a2_->at(gp) = a;
  }

  // unpack history
  ExtractfromPack(position,data,numhist);
//  histstresscurr_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
//  artstresscurr_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
//  histstresslast_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
//  artstresslast_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
  histstresscurr_=rcp(new vector<Epetra_SerialDenseVector>);
  artstresscurr_=rcp(new vector<Epetra_SerialDenseVector>);
  histstresslast_=rcp(new vector<Epetra_SerialDenseVector>);
  artstresslast_=rcp(new vector<Epetra_SerialDenseVector>);
  for (int var=0; var<numhist; var++)
  {
    // current vectors have to be initialized
//    LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> tmp;
    Epetra_SerialDenseVector tmp(NUM_STRESS_3D);
    histstresscurr_->push_back(tmp);
    artstresscurr_->push_back(tmp);

    // last vectors are unpacked
    ExtractfromPack(position,data,tmp);
    histstresslast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    artstresslast_->push_back(tmp);
  }


  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ViscoAnisotropic::Setup(const int numgp)
{
  a1_ = rcp(new vector<vector<double> > (numgp));
  a2_ = rcp(new vector<vector<double> > (numgp));
  if ((matdata_->m.viscoanisotropic->gamma<0) || (matdata_->m.viscoanisotropic->gamma >90)) dserror("Fiber angle not in [0,90]");
  const double gamma = (matdata_->m.viscoanisotropic->gamma*PI)/180.; //convert

  /* fibers are always related to a local element cosy
     which has to be specified in the element line */
  int ierr=0;
  int error = 1;
  // read local (cylindrical) cosy-directions at current element
  vector<double> rad(3);
  vector<double> axi(3);
  vector<double> cir(3);
  frdouble_n("RAD",&rad[0],3,&ierr); error*=ierr;
  frdouble_n("AXI",&axi[0],3,&ierr); error*=ierr;
  frdouble_n("CIR",&cir[0],3,&ierr); error*=ierr;
  if (error == 0) dserror("Problems with reading element local cosy, check RAD,AXI,CIR!");

  LINALG::FixedSizeSerialDenseMatrix<3,3> locsys;
  // basis is local cosy with third vec e3 = circumferential dir and e2 = axial dir
  double radnorm=0.; double axinorm=0.; double cirnorm=0.;
  for (int i = 0; i < 3; ++i) {
    radnorm += rad[i]*rad[i]; axinorm += axi[i]*axi[i]; cirnorm += cir[i]*cir[i];
  }
  radnorm = sqrt(radnorm); axinorm = sqrt(axinorm); cirnorm = sqrt(cirnorm);
  for (int i=0; i<3; ++i){
    locsys(i,0) = rad[i]/radnorm;
    locsys(i,1) = axi[i]/axinorm;
    locsys(i,2) = cir[i]/cirnorm;
  }
  for (int gp = 0; gp < numgp; ++gp) {
    a1_->at(gp).resize(3);
    a2_->at(gp).resize(3);
    for (int i = 0; i < 3; ++i) {
      // a1 = cos gamma e3 + sin gamma e2
      a1_->at(gp)[i] = cos(gamma)*locsys(i,2) + sin(gamma)*locsys(i,1);
      // a2 = cos gamma e3 - sin gamma e2
      a2_->at(gp)[i] = cos(gamma)*locsys(i,2) - sin(gamma)*locsys(i,1);
    }
  }

  // some input checking of viscous parameters for convenience
  if ((matdata_->m.viscoanisotropic->beta[0]<0) || (matdata_->m.viscoanisotropic->beta[1]<0)
     || (matdata_->m.viscoanisotropic->relax[0]<=0) || (matdata_->m.viscoanisotropic->relax[1]<=0))
        dserror("Check visocus parameters! Found beta < 0 or relax <= 0!");

  // initialize hist variables
//  histstresscurr_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
//  artstresscurr_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
//  histstresslast_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
//  artstresslast_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
//  const LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> emptyvec;
  
  histstresscurr_=rcp(new vector<Epetra_SerialDenseVector>);
  artstresscurr_=rcp(new vector<Epetra_SerialDenseVector>);
  histstresslast_=rcp(new vector<Epetra_SerialDenseVector>);
  artstresslast_=rcp(new vector<Epetra_SerialDenseVector>);
  const Epetra_SerialDenseVector emptyvec(NUM_STRESS_3D);
  // how many stress types are used?
  const int numst = matdata_->m.viscoanisotropic->numstresstypes;
  histstresscurr_->resize(numst*numgp);
  histstresslast_->resize(numst*numgp);
  artstresscurr_->resize(numst*numgp);
  artstresslast_->resize(numst*numgp);
  for (int j=0; j<numst*numgp; ++j)
  {
    histstresscurr_->at(j) = emptyvec;
    histstresslast_->at(j) = emptyvec;
    artstresscurr_->at(j) = emptyvec;
    artstresslast_->at(j) = emptyvec;
  }

  isinit_ = true;
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ViscoAnisotropic::Update()
{
  histstresslast_=histstresscurr_;
  artstresslast_=artstresscurr_;
//  const LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> emptyvec;//6 stresses for 3D
//  histstresscurr_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
//  artstresscurr_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
  const Epetra_SerialDenseVector emptyvec(NUM_STRESS_3D);//6 stresses for 3D
  histstresscurr_=rcp(new vector<Epetra_SerialDenseVector>);
  artstresscurr_=rcp(new vector<Epetra_SerialDenseVector>);
  const int histsize=histstresslast_->size();
  histstresscurr_->resize(histsize);
  artstresscurr_->resize(histsize);
  for (int j=0; j<histsize; ++j)
  {
    histstresscurr_->at(j) = emptyvec;
    artstresscurr_->at(j) = emptyvec;
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ViscoAnisotropic::Evaluate
(
  const Epetra_SerialDenseVector* glstrain,
  const int gp,
  Teuchos::ParameterList& params,
  Epetra_SerialDenseMatrix* cmat,
  Epetra_SerialDenseVector* stress
)
{
  const double mue = matdata_->m.viscoanisotropic->mue;
  const double kappa = matdata_->m.viscoanisotropic->kappa;
  const double k1 = matdata_->m.viscoanisotropic->k1;
  const double k2 = matdata_->m.viscoanisotropic->k2;

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
  Epetra_SerialDenseVector SisoEla_nh(NUM_STRESS_3D);  // isochoric elastic S from NeoHooke
  const double third = 1./3.;
  const double p = kappa*(J-1);
  for (int i = 0; i < 6; ++i) {
    (*stress)(i) = J*p * Cinv(i);  // volumetric part, not affected by viscosity
    SisoEla_nh(i) = incJ* (mue*Id(i) - third*mue*I1*Cinv(i));
  }

  // Elasticity =  Cvol + Ciso, via projection see Holzapfel p. 255
  //             + viscous C

  // Cvol = J(p + J dp/dJ) Cinv x Cinv  -  2 J p Cinv o Cinv
  // Ciso = 0 + 2/3 J^{-2/3} Sbar:C Psl - 2/3 (Cinv x Siso + Siso x Cinv)
  // Cvol not affected by viscosity
  Epetra_SerialDenseMatrix CisoEla_nh(NUM_STRESS_3D,NUM_STRESS_3D); // isochoric elastic C from NeoHooke

  AddtoCmatHolzapfelProduct((*cmat),Cinv,(-2*J*p));  // -2 J p Cinv o Cinv

  const double fac = 2*third*incJ*mue*I1;  // 2/3 J^{-2/3} Sbar:C
  // fac Psl = fac (Cinv o Cinv) - fac/3 (Cinv x Cinv)

  Epetra_SerialDenseMatrix Psl(6,6);        // Psl = Cinv o Cinv - 1/3 Cinv x Cinv
  AddtoCmatHolzapfelProduct(Psl,Cinv,1.0);  // first part Psl = Cinv o Cinv

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      double Siso_i = incJ* (mue*Id(i) - third*mue*I1*Cinv(i));
      double Siso_j = incJ* (mue*Id(j) - third*mue*I1*Cinv(j));
      (*cmat)(i,j) += J*(p+J*kappa) * Cinv(i) * Cinv(j);  // J(p + J dp/dJ) Cinv x Cinv
      Psl(i,j) += (-third) * Cinv(i) * Cinv(j);    // on the fly complete Psl needed later
      CisoEla_nh(i,j) = fac * Psl(i,j)                            // fac Psl
                      - 2*third * Cinv(i) * Siso_j                // -2/3 Cinv x Siso
                      - 2*third * Cinv(j) * Siso_i;               // -2/3 Siso x Cinv
    }
  }

  // anisotropic part: ***************************************************
  // W_aniso=(k1/(2.0*k2))*(exp(k2*pow((Ibar_{4,6} - 1.0),2)-1.0)); fiber SEF

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
  // however factor 2 for shear terms is already in C
  const double J4 = incJ * ( A1(0)*C(0) + A1(1)*C(1) + A1(2)*C(2)
                    + 1.*(A1(3)*C(3) + A1(4)*C(4) + A1(5)*C(5))); //J4 = trace(A1:C^dev)
  const double J6 = incJ * ( A2(0)*C(0) + A2(1)*C(1) + A2(2)*C(2)
                    + 1.*(A2(3)*C(3) + A2(4)*C(4) + A2(5)*C(5))); //J6 = trace(A2:C^dev)
  const double exp1 = exp(k2*(J4-1.)*(J4-1.));
  const double exp2 = exp(k2*(J6-1.)*(J6-1.));

  // fibers take compression only
  double fib1_tension = 1.;
  double fib2_tension = 1.;
  if (J4 < 1.0) fib1_tension = 0.;
  if (J6 < 1.0) fib2_tension = 0.;

  // PK2 fiber part in splitted formulation, see Holzapfel p. 271
  Epetra_SerialDenseVector SisoEla_fib1(A1); // first compute Sfbar1 = dWf/dJ4 A1
  Epetra_SerialDenseVector SisoEla_fib2(A2); // first compute Sfbar2 = dWf/dJ6 A2
  const double fib1 = fib1_tension* 2.*(k1*(J4-1.)*exp1);  // 2 dWf/dJ4
  const double fib2 = fib2_tension* 2.*(k1*(J6-1.)*exp2);  // 2 dWf/dJ6
  SisoEla_fib1.Scale(fib1);
  SisoEla_fib2.Scale(fib2);

  const double traceCSfbar1 =  SisoEla_fib1(0)*C(0) + SisoEla_fib1(1)*C(1) + SisoEla_fib1(2)*C(2)
                 + 1.*(SisoEla_fib1(3)*C(3) + SisoEla_fib1(4)*C(4) + SisoEla_fib1(5)*C(5)); // trace(Sfbar1 C)
  const double traceCSfbar2 =  SisoEla_fib2(0)*C(0) + SisoEla_fib2(1)*C(1) + SisoEla_fib2(2)*C(2)
                 + 1.*(SisoEla_fib2(3)*C(3) + SisoEla_fib2(4)*C(4) + SisoEla_fib2(5)*C(5)); // trace(Sfbar2 C)
  // compute Sfiso_a = J^{-2/3} * (Sfbar_a - 1/3 trace(Sfbar_a C) Cinv
  for (int i = 0; i < 6; ++i) {
    SisoEla_fib1(i) = incJ * (SisoEla_fib1(i) - third*traceCSfbar1*Cinv(i));
    SisoEla_fib2(i) = incJ * (SisoEla_fib2(i) - third*traceCSfbar2*Cinv(i));
  }

  // Elasticity fiber part in splitted formulation, see Holzapfel p. 255 and 272
  const double delta7bar1 = fib1_tension* 4.*(k1*exp1 + 2.*k1*k2*(J4-1.)*(J4-1.)*exp1); // 4 d^2Wf/dJ4dJ4
  const double delta7bar2 = fib2_tension* 4.*(k1*exp2 + 2.*k1*k2*(J6-1.)*(J6-1.)*exp2); // 4 d^2Wf/dJ6dJ6

  Epetra_SerialDenseMatrix CisoEla_fib1(NUM_STRESS_3D,NUM_STRESS_3D); // isochoric elastic C from Fib1
  Epetra_SerialDenseMatrix CisoEla_fib2(NUM_STRESS_3D,NUM_STRESS_3D); // isochoric elastic C from Fib2

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      double A1iso_i = incJ*A1(i)-third*J4*Cinv(i);  // A1iso = J^{-2/3} A1 - 1/3 J4 Cinv
      double A1iso_j = incJ*A1(j)-third*J4*Cinv(j);
      double A2iso_i = incJ*A2(i)-third*J6*Cinv(i);  // A2iso = J^{-2/3} A2 - 1/3 J6 Cinv
      double A2iso_j = incJ*A2(j)-third*J6*Cinv(j);
      CisoEla_fib1(i,j) = delta7bar1 * A1iso_i * A1iso_j  // delta7bar1 A1iso x A1iso
                        + 2.*third*incJ*traceCSfbar1 * Psl(i,j)  // 2/3 J^{-2/3} trace(Sfbar C) Psl
                        - 2.*third* (Cinv(i) * SisoEla_fib1(j) + Cinv(j) * SisoEla_fib1(i)); // -2/3 (Cinv x Sfiso1 + Sfiso1 x Cinv)
      CisoEla_fib2(i,j) =  delta7bar2 * A2iso_i * A2iso_j  // delta7bar2 A2iso x A2iso
                        + 2.*third*incJ*traceCSfbar2 * Psl(i,j)  // 2/3 J^{-2/3} trace(Sfbar2 C) Psl
                        - 2.*third* (Cinv(i) * SisoEla_fib2(j) + Cinv(j) * SisoEla_fib2(i)); // -2/3 (Cinv x Sfiso2 + Sfiso2 x Cinv)
    }
  }


  /* -------[-------  Viscous Part -------[-------
   * based on paper "A viscoelastic model for fiber-reinforced composites..." Holzapfel&Gasser, CMAME 2001
   */

  /* Layout for history vectors:
   * [SisoEla_nh(gp1) SisoEla_fib1(gp1) SisoEla_fib2(gp1)  SisoEla_nh(gp2) SisoEla_fib1(gp2)... ]
   * and the same for Q's
   */

  // read history
  const int numst = matdata_->m.viscoanisotropic->numstresstypes;
  Epetra_SerialDenseVector SisoEla_nh_old (histstresslast_->at(numst*gp   + 0));
  Epetra_SerialDenseVector SisoEla_fib1_old (histstresslast_->at(numst*gp + 1));
  Epetra_SerialDenseVector SisoEla_fib2_old (histstresslast_->at(numst*gp + 2));
  Epetra_SerialDenseVector Q_nh_old (artstresslast_->at(numst*gp   + 0));
  Epetra_SerialDenseVector Q_fib1_old (artstresslast_->at(numst*gp + 1));
  Epetra_SerialDenseVector Q_fib2_old (artstresslast_->at(numst*gp + 2));

  // visco parameters
  /*
   * the betas control the ratio serial/parallel elasticity of the Maxwell-body
   * the taus control the relaxation time for each dashpot path with respect to the serial elasticity
   */
  const double beta_nh = matdata_->m.viscoanisotropic->beta[0];
  const double beta_fib = matdata_->m.viscoanisotropic->beta[1];   // assume same beta for both fibers
  const double tau_nh = matdata_->m.viscoanisotropic->relax[0];
  const double tau_fib = matdata_->m.viscoanisotropic->relax[1];   // assume same tau for both fibers

  // get time algorithmic parameters
  double dt = params.get("delta time",-1.0);

  
  /* Time integration according Holzapfel paper
   * with some intrinsic usage of the analytic solution
   * of the underlying ODE (exp-function) */
  /*
  // evaluate exp-factors
  const double expfac_nh = exp(-dt*0.5/tau_nh);
  const double expfac_fib = exp(-dt*0.5/tau_fib);

  // evaluate current Q's
  Epetra_SerialDenseVector Q_nh(SisoEla_nh);
  Q_nh.Scale(beta_nh*expfac_nh);
  Epetra_SerialDenseVector Q_fib1(SisoEla_fib1);
  Q_fib1.Scale(beta_fib*expfac_fib);
  Epetra_SerialDenseVector Q_fib2(SisoEla_fib2);
  Q_fib2.Scale(beta_fib*expfac_fib);

  // evaluate 'H' history summands
  SisoEla_nh_old.Scale(-beta_nh*expfac_nh);
  Q_nh_old.Scale(expfac_nh*expfac_nh);
  SisoEla_fib1_old.Scale(-beta_fib*expfac_fib);
  Q_fib1_old.Scale(expfac_fib*expfac_fib);
  SisoEla_fib2_old.Scale(-beta_fib*expfac_fib);
  Q_fib2_old.Scale(expfac_fib*expfac_fib);
  */
  
  
  /* Time integration according Zien/Taylor and the viscoNeoHooke */
  const double theta = 0.5;
  const double artscalar1_nh=(tau_nh - dt + theta*dt)/tau_nh;
  const double artscalar2_nh=tau_nh/(tau_nh + theta*dt);
  const double artscalar1_fib=(tau_fib - dt + theta*dt)/tau_fib;
  const double artscalar2_fib=tau_fib/(tau_fib + theta*dt);

  
  // evaluate current Q's
  Epetra_SerialDenseVector Q_nh(SisoEla_nh);
  Q_nh.Scale(artscalar2_nh*beta_nh);
  Epetra_SerialDenseVector Q_fib1(SisoEla_fib1);
  Q_fib1.Scale(beta_fib*artscalar2_fib);
  Epetra_SerialDenseVector Q_fib2(SisoEla_fib2);
  Q_fib2.Scale(beta_fib*artscalar2_fib);
  
  // scale history
  SisoEla_nh_old.Scale(-beta_nh*artscalar2_nh);
  Q_nh_old.Scale(artscalar1_nh*artscalar2_nh);
  SisoEla_fib1_old.Scale(-beta_fib*artscalar2_fib);
  Q_fib1_old.Scale(artscalar1_fib*artscalar2_fib);
  SisoEla_fib2_old.Scale(-beta_fib*artscalar2_fib);
  Q_fib2_old.Scale(artscalar1_fib*artscalar2_fib);
  
  
  /* evaluate current stress */

  // elastic part
  (*stress) += SisoEla_nh;    // S_{n+1} = S_vol^ela + S_iso^ela
  (*stress) += SisoEla_fib1;
  (*stress) += SisoEla_fib2;  // S_{n+1} = S_vol^ela + S_iso^ela

  // viscous part
  Q_nh += Q_nh_old;
  Q_nh += SisoEla_nh_old;  // H_nh = Q_nh_old + S_nh_old
  (*stress) += Q_nh;            // S_{n+1} += H_nh + Q_nh_{n+1}

  Q_fib1 += Q_fib1_old;
  Q_fib1 += SisoEla_fib1_old;  // H_fib1 = Q_fib1_old + S_fib1_old
  (*stress) += Q_fib1;            // S_{n+1} += H_fib1 + Q_fib1_{n+1}

  Q_fib2 += Q_fib2_old;
  Q_fib2 += SisoEla_fib2_old;  // H_fib2 = Q_fib2_old + S_fib2_old
  (*stress) += Q_fib2;            // S_{n+1} += H_fib2 + Q_fib2_{n+1}

  /* evaluate current C-mat */
 
  /* Time integration according Holzapfel paper */
  /*
  CisoEla_nh.Scale(1.+beta_nh*expfac_nh);
  CisoEla_fib1.Scale(1.+beta_fib*expfac_fib);
  CisoEla_fib2.Scale(1.+beta_fib*expfac_fib);
  */
  
  
  /* Time integration according Zien/Taylor and the viscoNeoHooke */
  CisoEla_nh.Scale(1+beta_nh*artscalar2_nh);
  CisoEla_fib1.Scale(1+beta_fib*artscalar2_fib);
  CisoEla_fib2.Scale(1+beta_fib*artscalar2_fib);

  
  (*cmat) += CisoEla_nh;
  (*cmat) += CisoEla_fib1;
  (*cmat) += CisoEla_fib2;

  // update history
  histstresscurr_->at(numst*gp + 0) = SisoEla_nh;
  histstresscurr_->at(numst*gp + 1) = SisoEla_fib1;
  histstresscurr_->at(numst*gp + 2) = SisoEla_fib2;
  artstresscurr_->at(numst*gp + 0) = Q_nh;
  artstresscurr_->at(numst*gp + 1) = Q_fib1;
  artstresscurr_->at(numst*gp + 2) = Q_fib2;
  
  return;
}


#endif

