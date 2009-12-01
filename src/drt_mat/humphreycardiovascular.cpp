/*!----------------------------------------------------------------------
\file humphreycardiovascular.cpp
\brief
This file contains routines for an anisotropic material with amorphous elastin,
four fiber families of collagen (axial, circumferential, diagonal) and one
fiber family of smooth muscle (circumferential).
The strain energy function for one fiber family is the same as in holzapfelcardiovascular.

example input line
MAT 1 MAT_HUMPHREYCARDIO KAPPA 0.833 MUE 0.385 DENS 1.0 K1C 1.0 K2C 1.0 K1M 1.0 K2M 1.0 PHIE 0.02 PHIC 0.22 PHIM 0.76

<pre>
Maintainer: Susanna Tinkl
            tinkl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <vector>
#include "humphreycardiovascular.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::HumphreyCardio::HumphreyCardio(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  kappa_(matdata->GetDouble("KAPPA")),
  mue_(matdata->GetDouble("MUE")),
  density_(matdata->GetDouble("DENS")),
  k1c_(matdata->GetDouble("K1C")),
  k2c_(matdata->GetDouble("K2C")),
  k1m_(matdata->GetDouble("K1M")),
  k2m_(matdata->GetDouble("K2M")),
  phie_(matdata->GetDouble("PHIE")),
  phic_(matdata->GetDouble("PHIC")),
  phim_(matdata->GetDouble("PHIM"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         11/09|
 *----------------------------------------------------------------------*/
MAT::HumphreyCardio::HumphreyCardio()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          11/09|
 *----------------------------------------------------------------------*/
MAT::HumphreyCardio::HumphreyCardio(MAT::PAR::HumphreyCardio* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)         11/09|
 *----------------------------------------------------------------------*/
void MAT::HumphreyCardio::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);

  int numgp;
  if (!isinit_)
  {
    numgp = 0; // not initialized -> nothing to pack
  }
  else
  {
    numgp = a1_->size();   // size is number of gausspoints
  }
  AddtoPack(data,numgp);
  // Pack internal variables
  for (int gp = 0; gp < numgp; ++gp)
  {
    AddtoPack(data,a1_->at(gp));
    AddtoPack(data,a2_->at(gp));
    AddtoPack(data,a3_->at(gp));
    AddtoPack(data,a4_->at(gp));
    AddtoPack(data,ca1_->at(gp));
    AddtoPack(data,ca2_->at(gp));
    AddtoPack(data,ca3_->at(gp));
    AddtoPack(data,ca4_->at(gp));
  }
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)         11/09|
 *----------------------------------------------------------------------*/
void MAT::HumphreyCardio::Unpack(const vector<char>& data)
{
  isinit_=true;
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  // in post-process mode we do not have any instance of DRT::Problem
  if (DRT::Problem::NumInstances() > 0)
  {
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
    MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
    if (mat->Type() == MaterialType())
      params_ = static_cast<MAT::PAR::HumphreyCardio*>(mat);
    else
      dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
  }
  else
  {
    params_ = NULL;
  }

  int numgp;
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
  a3_ = rcp(new vector<vector<double> >(numgp));
  a4_ = rcp(new vector<vector<double> >(numgp));
  ca1_ = rcp(new vector<vector<double> >(numgp));
  ca2_ = rcp(new vector<vector<double> >(numgp));
  ca3_ = rcp(new vector<vector<double> >(numgp));
  ca4_ = rcp(new vector<vector<double> >(numgp));
  
  for (int gp = 0; gp < numgp; ++gp) {
    vector<double> a;
    ExtractfromPack(position,data,a);
    a1_->at(gp) = a;
    ExtractfromPack(position,data,a);
    a2_->at(gp) = a;
    ExtractfromPack(position,data,a);
    a3_->at(gp) = a;
    ExtractfromPack(position,data,a);
    a4_->at(gp) = a;
    ExtractfromPack(position,data,a);
    ca1_->at(gp) = a;
    ExtractfromPack(position,data,a);
    ca2_->at(gp) = a;
    ExtractfromPack(position,data,a);
    ca3_->at(gp) = a;
    ExtractfromPack(position,data,a);
    ca4_->at(gp) = a;
  }

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 |  Setup                                         (public)         11/09|
 *----------------------------------------------------------------------*/
void MAT::HumphreyCardio::Setup(const int numgp, DRT::INPUT::LineDefinition* linedef)
{

  /*fiber directions has to be defined in the element line
    Fibers are related to a local element cosy which has to be
    specified in the element line */

  a1_ = rcp(new vector<vector<double> > (numgp));
  a2_ = rcp(new vector<vector<double> > (numgp));
  a3_ = rcp(new vector<vector<double> > (numgp));
  a4_ = rcp(new vector<vector<double> > (numgp));
  ca1_ = rcp(new vector<vector<double> > (numgp));
  ca2_ = rcp(new vector<vector<double> > (numgp));
  ca3_ = rcp(new vector<vector<double> > (numgp));
  ca4_ = rcp(new vector<vector<double> > (numgp));
  
  const double gamma = (45*PI)/180.; //angle for diagonal fibers

  // read local (cylindrical) cosy-directions at current element
  vector<double> rad;
  vector<double> axi;
  vector<double> cir;
  linedef->ExtractDoubleVector("RAD",rad);
  linedef->ExtractDoubleVector("AXI",axi);
  linedef->ExtractDoubleVector("CIR",cir);

  LINALG::Matrix<3,3> locsys;
  // basis is local cosy with third vec e3 = circumferential dir and e2 = axial dir
  double radnorm=0.; double axinorm=0.; double cirnorm=0.;
  for (int i = 0; i < 3; i++) {
    radnorm += rad[i]*rad[i]; axinorm += axi[i]*axi[i]; cirnorm += cir[i]*cir[i];
  }
  radnorm = sqrt(radnorm); axinorm = sqrt(axinorm); cirnorm = sqrt(cirnorm);
  for (int i=0; i<3; i++){
    locsys(i,0) = rad[i]/radnorm;
    locsys(i,1) = axi[i]/axinorm;
    locsys(i,2) = cir[i]/cirnorm;
  }
  for (int gp = 0; gp < numgp; gp++) {
    a1_->at(gp).resize(3);
    a2_->at(gp).resize(3);
    a3_->at(gp).resize(3);
    a4_->at(gp).resize(3);
    ca1_->at(gp).resize(3);
    ca2_->at(gp).resize(3);
    ca3_->at(gp).resize(3);
    ca4_->at(gp).resize(3);
    for (int i = 0; i < 3; i++) {
      // a1 = e3, circumferential direction, used for collagen and smooth muscle
      a1_->at(gp)[i] = locsys(i,2);
      // a2 = e2
      a2_->at(gp)[i] = locsys(i,1);
      // a3 = cos gamma e3 + sin gamma e2
      a3_->at(gp)[i] = cos(gamma)*locsys(i,2) + sin(gamma)*locsys(i,1);
      // a4 = cos gamma e3 - sin gamma e2
      a4_->at(gp)[i] = cos(gamma)*locsys(i,2) - sin(gamma)*locsys(i,1);
      ca1_->at(gp)[i] = a1_->at(gp)[i];
      ca2_->at(gp)[i] = a2_->at(gp)[i];
      ca3_->at(gp)[i] = a3_->at(gp)[i];
      ca4_->at(gp)[i] = a4_->at(gp)[i];
    }
  }

  isinit_ = true;
  return;
}

/*----------------------------------------------------------------------*
 |  UpdateFiberDirs                               (public)         11/09|
 *----------------------------------------------------------------------*/
void MAT::HumphreyCardio::UpdateFiberDirs(const int gp, LINALG::Matrix<3,3>* defgrad)
{
  //Loop over all gp and update fiber directions
  ca1_->at(gp).resize(3);
  ca2_->at(gp).resize(3);
  ca3_->at(gp).resize(3);
  ca4_->at(gp).resize(3);
  LINALG::DENSEFUNCTIONS::multiply<3,3,1>(&((ca1_->at(gp))[0]),defgrad->A(),&((a1_->at(gp))[0]));
  LINALG::DENSEFUNCTIONS::multiply<3,3,1>(&((ca2_->at(gp))[0]),defgrad->A(),&((a2_->at(gp))[0]));
  LINALG::DENSEFUNCTIONS::multiply<3,3,1>(&((ca3_->at(gp))[0]),defgrad->A(),&((a3_->at(gp))[0]));
  LINALG::DENSEFUNCTIONS::multiply<3,3,1>(&((ca4_->at(gp))[0]),defgrad->A(),&((a4_->at(gp))[0]));
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)         11/09|
 *----------------------------------------------------------------------*
 strain energy function

 W    = phie*1/2 mue (I1*J^(-1/3)-3) + phif*(k1/(2.0*k2))*(exp(k2*(I_fib*J^(-2/3) - 1.0)^2-1.0))

 taken from
 L. Cardamone, A. Valentin, J.F. Eberth, J.D. Humphrey: Origin of axial prestretch
 and residual stress in arteries, Biomech Model Mechanobiol 8 (2009) 431-446
 (to be consistent with holzapfelcardiovascular we use here k1 instead of c2 = 2 k1
 which is used in the paper above)

 here

 I1    .. first invariant of right Cauchy-Green tensor C
 I_fib .. invariants accounting for the fiber directions
 J     .. det(F) determinant of the Jacobian matrix
 phie  .. mass fraction of elastin
 phif  .. mass fraction of fiber family

 The volumetric part is done by a volumetric strain energy function

 W_vol = 1/2 kappa (J-1)^2
 
 important: Cardamone et al used a Lagrange multiplier method to enforce
 incompressibility, here we add a volumetric strain energy function
 and use modified invariants
  
 */

void MAT::HumphreyCardio::Evaluate
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  const int gp,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> * cmat,
  LINALG::Matrix<NUM_STRESS_3D,1> * stress
)
{
  const double mue = params_->mue_;
  const double kappa = params_->kappa_;
  const double k1c = params_->k1c_;
  const double k2c = params_->k2c_;
  double k1m = params_->k1m_;
  const double k2m = params_->k2m_;
  const double phie = params_->phie_;
  const double phic = params_->phic_;
  const double phim = params_->phim_;

  //--------------------------------------------------------------------------------------
  // build identity tensor I
  LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
  C.Scale(2.0);
  C += Id;

  // isotropic invariants
  const double I1 = C(0) + C(1) + C(2);  // 1st invariant, trace
  const double I3 = C(0)*C(1)*C(2)
        + 0.25 * C(3)*C(4)*C(5)
        - 0.25 * C(1)*C(5)*C(5)
        - 0.25 * C(2)*C(3)*C(3)
        - 0.25 * C(0)*C(4)*C(4);    // 3rd invariant, determinant
  const double J = sqrt(I3);     // determinant of F
  const double incJ = pow(I3,-1.0/3.0);  // J^{-2/3}

  //--------------------------------------------------------------------------------------
  // invert C
  LINALG::Matrix<NUM_STRESS_3D,1> Cinv(6);

  Cinv(0) = C(1)*C(2) - 0.25*C(4)*C(4);
  Cinv(1) = C(0)*C(2) - 0.25*C(5)*C(5);
  Cinv(2) = C(0)*C(1) - 0.25*C(3)*C(3);
  Cinv(3) = 0.25*C(5)*C(4) - 0.5*C(3)*C(2);
  Cinv(4) = 0.25*C(3)*C(5) - 0.5*C(0)*C(4);
  Cinv(5) = 0.25*C(3)*C(4) - 0.5*C(5)*C(1);
  Cinv.Scale(1.0/I3);
  
  //--------------------------------------------------------------------------------------
  // structural tensors in voigt notation
  // A1 = a1 x a1, A2 = a2 x a2, A3 = a3 x a3, A4 = a4 x a4
  LINALG::Matrix<NUM_STRESS_3D,1>  A1;
  LINALG::Matrix<NUM_STRESS_3D,1>  A2;
  LINALG::Matrix<NUM_STRESS_3D,1>  A3;
  LINALG::Matrix<NUM_STRESS_3D,1>  A4;
  for (int i = 0; i < 3; i++) {
    A1(i) = a1_->at(gp)[i]*a1_->at(gp)[i];
    A2(i) = a2_->at(gp)[i]*a2_->at(gp)[i];
    A3(i) = a3_->at(gp)[i]*a3_->at(gp)[i];
    A4(i) = a4_->at(gp)[i]*a4_->at(gp)[i];
  }
  A1(3) = a1_->at(gp)[0]*a1_->at(gp)[1]; A1(4) = a1_->at(gp)[1]*a1_->at(gp)[2]; A1(5) = a1_->at(gp)[0]*a1_->at(gp)[2];
  A2(3) = a2_->at(gp)[0]*a2_->at(gp)[1]; A2(4) = a2_->at(gp)[1]*a2_->at(gp)[2]; A2(5) = a2_->at(gp)[0]*a2_->at(gp)[2];
  A3(3) = a3_->at(gp)[0]*a3_->at(gp)[1]; A3(4) = a3_->at(gp)[1]*a3_->at(gp)[2]; A3(5) = a3_->at(gp)[0]*a3_->at(gp)[2];
  A4(3) = a4_->at(gp)[0]*a4_->at(gp)[1]; A4(4) = a4_->at(gp)[1]*a4_->at(gp)[2]; A4(5) = a4_->at(gp)[0]*a4_->at(gp)[2];

  // modified (fiber-) invariants J_{4,6,8,10} = J^{-2/3}*I_{4,6,8,10}
  // trace(AB) =  a11 b11 + 2 a12 b12 + 2 a13 b13 + a22 b22 + 2 a23 b23 + a33 b33
  // however factor 2 for shear terms is already in C
  double J4 = incJ * ( A1(0)*C(0) + A1(1)*C(1) + A1(2)*C(2)
                    + 1.*(A1(3)*C(3) + A1(4)*C(4) + A1(5)*C(5))); //J4 = trace(A1:C^dev)
  double J6 = incJ * ( A2(0)*C(0) + A2(1)*C(1) + A2(2)*C(2)
                    + 1.*(A2(3)*C(3) + A2(4)*C(4) + A2(5)*C(5))); //J6 = trace(A2:C^dev)
  double J8 = incJ * ( A3(0)*C(0) + A3(1)*C(1) + A3(2)*C(2)
                    + 1.*(A3(3)*C(3) + A3(4)*C(4) + A3(5)*C(5))); //J8 = trace(A3:C^dev)
  double J10 = incJ * ( A4(0)*C(0) + A4(1)*C(1) + A4(2)*C(2)
                    + 1.*(A4(3)*C(3) + A4(4)*C(4) + A4(5)*C(5))); //J10 = trace(A4:C^dev)
  
  //--------------------------------------------------------------------------------------
  // collagen fibers behave as smooth muscle in compression
  // additional build in of fiber mass fractions
  double k1c_fib1 = k1c*phic;
  double k2c_fib1 = k2c;
  double k1c_fib2 = k1c*phic;
  double k2c_fib2 = k2c;
  double k1c_fib3 = k1c*phic;
  double k2c_fib3 = k2c;
  double k1c_fib4 = k1c*phic;
  double k2c_fib4 = k2c;

  if (J4 < 1)
  {
    k1c_fib1=k1m*phic;
    k2c_fib1=k2m;
  }
  if (J6 < 1)
  {
    k1c_fib2=k1m*phic;
    k2c_fib2=k2m;
  }
  if (J8 < 1)
  {
    k1c_fib3=k1m*phic;
    k2c_fib3=k2m;
  }
  if (J10 < 1)
  {
    k1c_fib4=k1m*phic;
    k2c_fib4=k2m;
  }
  // consider mass fraction of smooth muscle
  k1m=k1m*phim;

  //--- prepare some constants -----------------------------------------------------------
  const double third = 1./3.;
  const double p = kappa*(J-1); // dW_vol/dJ
  
  
  //--- determine 2nd Piola Kirchhoff stresses S -------------------------------------
  // 1st step: isotropic part
  //=========================
  LINALG::Matrix<NUM_STRESS_3D,1> Siso;  // isotropic S
  for (int i = 0; i < 6; i++) {
	// volumetric part J*kappa*(J-1)*Cinv
    Siso(i) = J*p * Cinv(i);
    // isochoric part via projection J^(-2/3)*PP:Sbar, see Holzapfel p. 230
    // with Sbar = mue*Id
    Siso(i) += phie*incJ* (mue*Id(i) - third*mue*I1*Cinv(i));
  }
  
  // 2nd step: anisotropic part
  //===========================
  // PK2 fiber part in splitted formulation, see Holzapfel p. 271
  
  // collagen fiber families
  //------------------------
  LINALG::Matrix<NUM_STRESS_3D,1> Saniso_fib1(A1); // first compute Sfbar1 = 2 dW/dJ4 A1
  LINALG::Matrix<NUM_STRESS_3D,1> Saniso_fib2(A2); // first compute Sfbar2 = 2 dW/dJ6 A2
  LINALG::Matrix<NUM_STRESS_3D,1> Saniso_fib3(A3); // first compute Sfbar3 = 2 dW/dJ8 A3
  LINALG::Matrix<NUM_STRESS_3D,1> Saniso_fib4(A4); // first compute Sfbar4 = 2 dW/dJ10 A4  
  const double exp1 = exp(k2c_fib1*(J4-1.)*(J4-1.));
  const double exp2 = exp(k2c_fib2*(J6-1.)*(J6-1.));
  const double exp3 = exp(k2c_fib3*(J8-1.)*(J8-1.));
  const double exp4 = exp(k2c_fib4*(J10-1.)*(J10-1.));
  const double fib1 = 2.*(k1c_fib1*(J4-1.)*exp1);  // 2 dW/dJ4
  const double fib2 = 2.*(k1c_fib2*(J6-1.)*exp2);  // 2 dW/dJ6
  const double fib3 = 2.*(k1c_fib3*(J8-1.)*exp3);  // 2 dW/dJ8
  const double fib4 = 2.*(k1c_fib4*(J10-1.)*exp4);  // 2 dW/dJ10
  Saniso_fib1.Scale(fib1);  //Sfbar1
  Saniso_fib2.Scale(fib2);  //Sfbar2
  Saniso_fib3.Scale(fib3);  //Sfbar3
  Saniso_fib4.Scale(fib4);  //Sfbar4

  const double traceCSfbar1 =  Saniso_fib1(0)*C(0) + Saniso_fib1(1)*C(1) + Saniso_fib1(2)*C(2)
                 + 1.*(Saniso_fib1(3)*C(3) + Saniso_fib1(4)*C(4) + Saniso_fib1(5)*C(5)); // trace(Sfbar1 C)
  const double traceCSfbar2 =  Saniso_fib2(0)*C(0) + Saniso_fib2(1)*C(1) + Saniso_fib2(2)*C(2)
                 + 1.*(Saniso_fib2(3)*C(3) + Saniso_fib2(4)*C(4) + Saniso_fib2(5)*C(5)); // trace(Sfbar2 C)
  const double traceCSfbar3 =  Saniso_fib3(0)*C(0) + Saniso_fib3(1)*C(1) + Saniso_fib3(2)*C(2)
                 + 1.*(Saniso_fib3(3)*C(3) + Saniso_fib3(4)*C(4) + Saniso_fib3(5)*C(5)); // trace(Sfbar3 C)
  const double traceCSfbar4 =  Saniso_fib4(0)*C(0) + Saniso_fib4(1)*C(1) + Saniso_fib4(2)*C(2)
                 + 1.*(Saniso_fib4(3)*C(3) + Saniso_fib4(4)*C(4) + Saniso_fib4(5)*C(5)); // trace(Sfbar4 C)
  // compute Sfiso_a = J^{-2/3} * (Sfbar_a - 1/3 trace(Sfbar_a C) Cinv)
  for (int i = 0; i < 6; i++) {
    Saniso_fib1(i) = incJ * (Saniso_fib1(i) - third*traceCSfbar1*Cinv(i));
    Saniso_fib2(i) = incJ * (Saniso_fib2(i) - third*traceCSfbar2*Cinv(i));
    Saniso_fib3(i) = incJ * (Saniso_fib3(i) - third*traceCSfbar3*Cinv(i));
    Saniso_fib4(i) = incJ * (Saniso_fib4(i) - third*traceCSfbar4*Cinv(i));
  }
  
  // smooth muscle fiber family
  //---------------------------
  LINALG::Matrix<NUM_STRESS_3D,1> Sm_fib(A1); // first compute Smbar = 2 dWm/dJ4 A1
  const double expm = exp(k2m*(J4-1.)*(J4-1.));
  const double fibm = 2.*(k1m*(J4-1.)*expm);  // 2 dWm/dJ4
  Sm_fib.Scale(fibm); // Smbar
  const double traceCSmbar =  Sm_fib(0)*C(0) + Sm_fib(1)*C(1) + Sm_fib(2)*C(2)
                 + 1.*(Sm_fib(3)*C(3) + Sm_fib(4)*C(4) + Sm_fib(5)*C(5)); // trace(Smbar C)
  for (int i = 0; i < 6; i++) {
	Sm_fib(i) = incJ * (Sm_fib(i) - third*traceCSmbar*Cinv(i));
  }
  
  
  // 3rd step: add everything up
  //============================
  (*stress) = Siso;    
  (*stress) += Saniso_fib1;
  (*stress) += Saniso_fib2;
  (*stress) += Saniso_fib3;
  (*stress) += Saniso_fib4;
  (*stress) += Sm_fib;
  
  
  //--- do elasticity matrix -------------------------------------------------------------
  // ensure that cmat is zero when it enters the computation
  // It is an implicit law that cmat is zero upon input
  //cmat.PutScalar(0.0);

  // 1st step: isotropic part
  //=========================
  // Elasticity =  cmatvol + cmatiso, via projection see Holzapfel p. 255
  // cmatvol = J(p + J dp/dJ) Cinv x Cinv  -  2 J p Cinv o Cinv
  // cmatiso = 0 + 2/3 J^{-2/3} Sbar:C Psl - 2/3 (Cinv x Siso + Siso x Cinv)

  AddtoCmatHolzapfelProduct((*cmat),Cinv,(-2*J*p));
  // the rest of cmatvol is computed during the loop
  
  const double fac = 2*third*incJ*mue*I1;  // 2/3 J^{-2/3} Sbar:C
  
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>  Psl(true);        // Psl = Cinv o Cinv - 1/3 Cinv x Cinv
  AddtoCmatHolzapfelProduct(Psl,Cinv,1.0);  // first part Psl = Cinv o Cinv

  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      // isochoric part of isotropic S (Siso)
      double Siso_i = incJ* (mue*Id(i) - third*mue*I1*Cinv(i));
      double Siso_j = incJ* (mue*Id(j) - third*mue*I1*Cinv(j));
      // missing part of cmatvol J(p + J dp/dJ) Cinv x Cinv
      (*cmat)(i,j) += J*(p+J*kappa) * Cinv(i) * Cinv(j);
      Psl(i,j) += (-third) * Cinv(i) * Cinv(j);    // on the fly complete Psl needed later
      (*cmat)(i,j) += phie*(fac * Psl(i,j)                            // fac Psl
                      - 2*third * Cinv(i) * Siso_j                // -2/3 Cinv x Siso
                      - 2*third * Cinv(j) * Siso_i);               // -2/3 Siso x Cinv
    }
  }

  // 2nd step: anisotropic part
  //===========================
  // Elasticity fiber part in splitted formulation, see Holzapfel p. 255 and 272
  
  // collagen fiber families
  //------------------------  
  const double delta7bar1 = 4.*(k1c_fib1*exp1 + 2.*k1c_fib1*k2c_fib1*(J4-1.)*(J4-1.)*exp1); // 4 d^2Wf/dJ4dJ4
  const double delta7bar2 = 4.*(k1c_fib2*exp2 + 2.*k1c_fib2*k2c_fib2*(J6-1.)*(J6-1.)*exp2); // 4 d^2Wf/dJ6dJ6
  const double delta7bar3 = 4.*(k1c_fib3*exp3 + 2.*k1c_fib3*k2c_fib3*(J8-1.)*(J8-1.)*exp3); // 4 d^2Wf/dJ8dJ8
  const double delta7bar4 = 4.*(k1c_fib4*exp4 + 2.*k1c_fib4*k2c_fib4*(J10-1.)*(J10-1.)*exp4); // 4 d^2Wf/dJ10dJ10

  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> Caniso_fib1; // isochoric elastic C from Fib1
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> Caniso_fib2; // isochoric elastic C from Fib2
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> Caniso_fib3; // isochoric elastic C from Fib3
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> Caniso_fib4; // isochoric elastic C from Fib4

  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      double A1iso_i = incJ*A1(i)-third*J4*Cinv(i);  // A1iso = J^{-2/3} A1 - 1/3 J4 Cinv
      double A1iso_j = incJ*A1(j)-third*J4*Cinv(j);  //       = J^(-2/3)*PP:A1 has no physical meaning
      double A2iso_i = incJ*A2(i)-third*J6*Cinv(i);  // A2iso = J^{-2/3} A2 - 1/3 J6 Cinv
      double A2iso_j = incJ*A2(j)-third*J6*Cinv(j);
      double A3iso_i = incJ*A3(i)-third*J8*Cinv(i);  // A3iso = J^{-2/3} A3 - 1/3 J8 Cinv
      double A3iso_j = incJ*A3(j)-third*J8*Cinv(j);
      double A4iso_i = incJ*A4(i)-third*J10*Cinv(i);  // A4iso = J^{-2/3} A4 - 1/3 J10 Cinv
      double A4iso_j = incJ*A4(j)-third*J10*Cinv(j);
      Caniso_fib1(i,j) = delta7bar1 * A1iso_i * A1iso_j  // delta7bar1 A1iso x A1iso
                        + 2.*third*incJ*traceCSfbar1 * Psl(i,j)  // 2/3 J^{-2/3} trace(Sfbar1 C) Psl
                        - 2.*third* (Cinv(i) * Saniso_fib1(j) + Cinv(j) * Saniso_fib1(i)); // -2/3 (Cinv x Sfiso1 + Sfiso1 x Cinv)
      Caniso_fib2(i,j) = delta7bar2 * A2iso_i * A2iso_j  // delta7bar2 A2iso x A2iso
                        + 2.*third*incJ*traceCSfbar2 * Psl(i,j)  // 2/3 J^{-2/3} trace(Sfbar2 C) Psl
                        - 2.*third* (Cinv(i) * Saniso_fib2(j) + Cinv(j) * Saniso_fib2(i)); // -2/3 (Cinv x Sfiso2 + Sfiso2 x Cinv)
      Caniso_fib3(i,j) = delta7bar3 * A3iso_i * A3iso_j  // delta7bar3 A3iso x A3iso
                        + 2.*third*incJ*traceCSfbar3 * Psl(i,j)  // 2/3 J^{-2/3} trace(Sfbar3 C) Psl
                        - 2.*third* (Cinv(i) * Saniso_fib3(j) + Cinv(j) * Saniso_fib3(i)); // -2/3 (Cinv x Sfiso3 + Sfiso3 x Cinv)
      Caniso_fib4(i,j) = delta7bar4 * A4iso_i * A4iso_j  // delta7bar4 A4iso x A4iso
                        + 2.*third*incJ*traceCSfbar4 * Psl(i,j)  // 2/3 J^{-2/3} trace(Sfbar4 C) Psl
                        - 2.*third* (Cinv(i) * Saniso_fib4(j) + Cinv(j) * Saniso_fib4(i)); // -2/3 (Cinv x Sfiso4 + Sfiso4 x Cinv)
    }
  }
  
  (*cmat) += Caniso_fib1;
  (*cmat) += Caniso_fib2;
  (*cmat) += Caniso_fib3;
  (*cmat) += Caniso_fib4;
  
  // smooth muscle fiber family
  //---------------------------
  const double delta7barm = 4.*(k1m*exp1 + 2.*k1m*k2m*(J4-1.)*(J4-1.)*expm); // 4 d^2Wm/dJ4dJ4
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> Cm_fib; // isochoric elastic C from Fib1
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      double A1isom_i = incJ*A1(i)-third*J4*Cinv(i);  // A1isom = J^{-2/3} A1 - 1/3 J4 Cinv
      double A1isom_j = incJ*A1(j)-third*J4*Cinv(j);  //       = J^(-2/3)*PP:A1 has no physical meaning
      Cm_fib(i,j) = delta7barm * A1isom_i * A1isom_j  // delta7barm A1isom x A1isom
                   + 2.*third*incJ*traceCSmbar * Psl(i,j)  // 2/3 J^{-2/3} trace(Smbar C) Psl
                   - 2.*third* (Cinv(i) * Sm_fib(j) + Cinv(j) * Sm_fib(i)); // -2/3 (Cinv x Smiso + Smiso x Cinv)
    }
  }
  
  (*cmat) += Cm_fib;
  
  return;
}


#endif
