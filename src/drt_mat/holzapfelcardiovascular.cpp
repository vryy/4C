/*!----------------------------------------------------------------------
\file holzapfelcardiovascular.cpp
\brief
This file contains routines for an anisotropic material with two fiber families.
example input line
MAT 1 MAT_HOLZAPFELCARDIO KAPPA 0.833 MUE 0.385 DENS 1.0 K1 1.0 K2 1.0 GAMMA 45.0 MINSTRETCH 1.0 INIT 1

<pre>
Maintainer: Susanna Tinkl
            tinkl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/

#include <vector>
#include "holzapfelcardiovascular.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::HolzapfelCardio::HolzapfelCardio(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  kappa_(matdata->GetDouble("KAPPA")),
  mue_(matdata->GetDouble("MUE")),
  density_(matdata->GetDouble("DENS")),
  k1_(matdata->GetDouble("K1")),
  k2_(matdata->GetDouble("K2")),
  gamma_(matdata->GetDouble("GAMMA")),
  minstretch_(matdata->GetDouble("MINSTRETCH")),
  init_(matdata->GetInt("INIT"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::HolzapfelCardio::CreateMaterial()
{
  return Teuchos::rcp(new MAT::HolzapfelCardio(this));
}



MAT::HolzapfelCardioType MAT::HolzapfelCardioType::instance_;


DRT::ParObject* MAT::HolzapfelCardioType::Create( const std::vector<char> & data )
{
  MAT::HolzapfelCardio* holzapfelcard = new MAT::HolzapfelCardio();
  holzapfelcard->Unpack(data);
  return holzapfelcard;
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         11/09|
 *----------------------------------------------------------------------*/
MAT::HolzapfelCardio::HolzapfelCardio()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          11/09|
 *----------------------------------------------------------------------*/
MAT::HolzapfelCardio::HolzapfelCardio(MAT::PAR::HolzapfelCardio* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)         11/09|
 *----------------------------------------------------------------------*/
void MAT::HolzapfelCardio::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

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
    AddtoPack(data,ca1_->at(gp));
    AddtoPack(data,ca2_->at(gp));
  }
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)         11/09|
 *----------------------------------------------------------------------*/
void MAT::HolzapfelCardio::Unpack(const std::vector<char>& data)
{
  isinit_=true;
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::HolzapfelCardio*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  int numgp;
  ExtractfromPack(position,data,numgp);
  if (numgp == 0){ // no history data to unpack
    isinit_=false;
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d",data.size(),position);
    return;
  }
  // unpack fiber internal variables
  a1_ = Teuchos::rcp(new std::vector<std::vector<double> >(numgp));
  a2_ = Teuchos::rcp(new std::vector<std::vector<double> >(numgp));
  ca1_ = Teuchos::rcp(new std::vector<std::vector<double> >(numgp));
  ca2_ = Teuchos::rcp(new std::vector<std::vector<double> >(numgp));

  for (int gp = 0; gp < numgp; ++gp) {
    std::vector<double> a;
    ExtractfromPack(position,data,a);
    a1_->at(gp) = a;
    ExtractfromPack(position,data,a);
    a2_->at(gp) = a;
    ExtractfromPack(position,data,a);
    ca1_->at(gp) = a;
    ExtractfromPack(position,data,a);
    ca2_->at(gp) = a;
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 |  Setup                                         (public)         11/09|
 *----------------------------------------------------------------------*/
void MAT::HolzapfelCardio::Setup(const int numgp, DRT::INPUT::LineDefinition* linedef)
{

  /*fiber directions has to be defined in the element line
    Fibers are related to a local element cosy which has to be
    specified in the element line */

  a1_ = Teuchos::rcp(new std::vector<std::vector<double> > (numgp));
  a2_ = Teuchos::rcp(new std::vector<std::vector<double> > (numgp));
  ca1_ = Teuchos::rcp(new std::vector<std::vector<double> > (numgp));
  ca2_ = Teuchos::rcp(new std::vector<std::vector<double> > (numgp));

  for (int gp = 0; gp < numgp; gp++) {
    a1_->at(gp).resize(3);
    a2_->at(gp).resize(3);
    ca1_->at(gp).resize(3);
    ca2_->at(gp).resize(3);
  }

  int initflag = params_->init_;

  if ((params_->gamma_<0) || (params_->gamma_ >90)) dserror("Fiber angle not in [0,90]");
  const double gamma = (params_->gamma_*PI)/180.; //convert

  if (initflag==0){
    // fibers aligned in YZ-plane with gamma around Z in global cartesian cosy
    LINALG::Matrix<3,3> id(true);
    // basis is identity
    for (int i=0; i<3; ++i) id(i,i) = 1.0;

    for (int gp = 0; gp < numgp; ++gp) EvaluateFiberVecs(gp,gamma,id,id);

  } else if (initflag==1){
	// read local (cylindrical) cosy-directions at current element
    std::vector<double> rad;
    std::vector<double> axi;
    std::vector<double> cir;
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
    LINALG::Matrix<3,3> Id(true);
    for (int i = 0; i < 3; i++) Id(i,i) = 1.0;

    for (int gp = 0; gp < numgp; gp++) EvaluateFiberVecs(gp,gamma,locsys,Id);

  } else if (initflag==3){
    // start with isotropic computation, thus fiber directions are set to zero
    // nothing has to be done
  } else dserror("INIT type not implemented");

  isinit_ = true;
  return;
}

/*----------------------------------------------------------------------*
 |  UpdateFiberDirs                               (public)         11/09|
 *----------------------------------------------------------------------*/
void MAT::HolzapfelCardio::UpdateFiberDirs(const int gp, LINALG::Matrix<3,3>* defgrad)
{
  //Loop over all gp and update fiber directions
  ca1_->at(gp).resize(3);
  ca2_->at(gp).resize(3);
  LINALG::DENSEFUNCTIONS::multiply<double,3,3,1>(&((ca1_->at(gp))[0]),defgrad->A(),&((a1_->at(gp))[0]));
  LINALG::DENSEFUNCTIONS::multiply<double,3,3,1>(&((ca2_->at(gp))[0]),defgrad->A(),&((a2_->at(gp))[0]));
  //cout << (ca1_->at(gp))[0] << ",  " << (ca1_->at(gp))[1] << ",  " << (ca1_->at(gp))[2] << endl;
  //cout <<  (a1_->at(gp))[0] << ",  " <<  (a1_->at(gp))[1] << ",  " <<  (a1_->at(gp))[2] << endl;
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)         11/09|
 *----------------------------------------------------------------------*
 strain energy function

 W    = 1/2 mue (I1*J^(-2/3)-3) + (k1/(2.0*k2))*(exp(k2*(I_{4,6}*J^(-2/3) - 1.0)^2-1.0))

 taken from
 G.A. Holzapfel, T.C. Gasser, R.W. Ogden: A new constitutive framework for arterial wall mechanics
 and a comparative study of material models, J. of Elasticity 61 (2000) 1-48.


 here

 I1      .. first invariant of right Cauchy-Green tensor C
 I_{4,6} .. invariants accounting for the fiber directions
 J       .. det(F) determinant of the Jacobian matrix

 The volumetric part is done by a volumetric strain energy function

 W_vol = 1/2 kappa (J-1)^2

 */

void MAT::HolzapfelCardio::Evaluate
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  const int gp,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> * cmat,
  LINALG::Matrix<NUM_STRESS_3D,1> * stress
)
{
  const double mue = params_->mue_;
  const double kappa = params_->kappa_;
  const double k1 = params_->k1_;
  const double k2 = params_->k2_;
  const double minstretch = params_->minstretch_;

  // here one can make a difference between output and normal call of material
  LINALG::Matrix<3,1> a1(true);
  LINALG::Matrix<3,1> a2(true);
  for (int i = 0; i < 3; i++) {
    a1(i) = a1_->at(gp)[i];
    a2(i) = a2_->at(gp)[i];
  }

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
  const double incJ = std::pow(I3,-1.0/3.0);  // J^{-2/3}

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
  // A1 = a1 x a1, A2 = a2 x a2
  LINALG::Matrix<NUM_STRESS_3D,1>  A1;
  LINALG::Matrix<NUM_STRESS_3D,1>  A2;
  for (int i = 0; i < 3; i++) {
    A1(i) = a1(i)*a1(i);
    A2(i) = a2(i)*a2(i);
  }
  A1(3) = a1(0)*a1(1); A1(4) = a1(1)*a1(2); A1(5) = a1(0)*a1(2);
  A2(3) = a2(0)*a2(1); A2(4) = a2(1)*a2(2); A2(5) = a2(0)*a2(2);

  // modified (fiber-) invariants J_{4,6} = J^{-2/3}*I_{4,6}
  // trace(AB) =  a11 b11 + 2 a12 b12 + 2 a13 b13 + a22 b22 + 2 a23 b23 + a33 b33
  // however factor 2 for shear terms is already in C
  double J4 = incJ * ( A1(0)*C(0) + A1(1)*C(1) + A1(2)*C(2)
                    + 1.*(A1(3)*C(3) + A1(4)*C(4) + A1(5)*C(5))); //J4 = trace(A1:C^dev)
  double J6 = incJ * ( A2(0)*C(0) + A2(1)*C(1) + A2(2)*C(2)
                    + 1.*(A2(3)*C(3) + A2(4)*C(4) + A2(5)*C(5))); //J6 = trace(A2:C^dev)

  //--------------------------------------------------------------------------------------
  // fibers can only stretch/compress down to a minimal value
  // for no compression minstretch has to be 1
  // also minstretch = 0 makes sense, if the fibers can be compressed
  double fib1_tension = 1.;
  double fib2_tension = 1.;

  if (J4 < minstretch)
  {
    J4 = minstretch;
    fib1_tension = 0.;
  }
  if (J6 < minstretch)
  {
    J6 = minstretch;
    fib2_tension = 0.;
  }

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
    Siso(i) += incJ* (mue*Id(i) - third*mue*I1*Cinv(i));
  }

  // 2nd step: anisotropic part
  //===========================
  // PK2 fiber part in splitted formulation, see Holzapfel p. 271
  LINALG::Matrix<NUM_STRESS_3D,1> Saniso_fib1(A1); // first compute Sfbar1 = 2 dW/dJ4 A1
  LINALG::Matrix<NUM_STRESS_3D,1> Saniso_fib2(A2); // first compute Sfbar2 = 2 dW/dJ6 A2
  const double exp1 = exp(k2*(J4-1.)*(J4-1.));
  const double exp2 = exp(k2*(J6-1.)*(J6-1.));
  const double fib1 = 2.*(k1*(J4-1.)*exp1);  // 2 dW/dJ4
  const double fib2 = 2.*(k1*(J6-1.)*exp2);  // 2 dW/dJ6
  Saniso_fib1.Scale(fib1);  //Sfbar1
  Saniso_fib2.Scale(fib2);  //Sfbar2

  const double traceCSfbar1 =  Saniso_fib1(0)*C(0) + Saniso_fib1(1)*C(1) + Saniso_fib1(2)*C(2)
                 + 1.*(Saniso_fib1(3)*C(3) + Saniso_fib1(4)*C(4) + Saniso_fib1(5)*C(5)); // trace(Sfbar1 C)
  const double traceCSfbar2 =  Saniso_fib2(0)*C(0) + Saniso_fib2(1)*C(1) + Saniso_fib2(2)*C(2)
                 + 1.*(Saniso_fib2(3)*C(3) + Saniso_fib2(4)*C(4) + Saniso_fib2(5)*C(5)); // trace(Sfbar2 C)
  // compute Sfiso_a = J^{-2/3} * (Sfbar_a - 1/3 trace(Sfbar_a C) Cinv)
  for (int i = 0; i < 6; i++) {
    Saniso_fib1(i) = incJ * (Saniso_fib1(i) - third*traceCSfbar1*Cinv(i));
    Saniso_fib2(i) = incJ * (Saniso_fib2(i) - third*traceCSfbar2*Cinv(i));
  }

  if (a1(0)==0 && a1(1)==0 && a1(2)==0){
    // isotropic fiber part for initial iteration step
    // W=(k1/(2.0*k2))*(exp(k2*pow((Ibar_1 - 3.0),2)-1.0));
    // the stress which is computed until now in the anisotropic part is zero
    const double expiso = exp(k2*(I1*incJ-3.)*(I1*incJ-3.));
    const double faciso = 2.*k1*(I1*incJ-3.)*expiso;
    for (int i = 0; i < 6; i++)
      Saniso_fib1(i) += incJ* faciso* (Id(i) - third*I1*Cinv(i));
  }

  // 3rd step: add everything up
  //============================
  (*stress) = Siso;
  (*stress) += Saniso_fib1;
  (*stress) += Saniso_fib2;


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
      (*cmat)(i,j) += fac * Psl(i,j)                            // fac Psl
                     - 2*third * Cinv(i) * Siso_j                // -2/3 Cinv x Siso
                     - 2*third * Cinv(j) * Siso_i;               // -2/3 Siso x Cinv
    }
  }

  // 2nd step: anisotropic part
  //===========================
  // check wether an initial isotropic step is needed
  if (a1(0)==0 && a1(1)==0 && a1(2)==0){
    // isotropic fiber part for initial iteration step
    // W=(k1/(2.0*k2))*(exp(k2*pow((Ibar_1 - 3.0),2)-1.0));
    const double expiso = exp(k2*(I1*incJ-3.)*(I1*incJ-3.));
    const double faciso = 2.*k1*(I1*incJ-3.)*expiso;
    const double delta7iso = incJ*incJ* 4.*(k1 + 2.*k1*k2*(I1*incJ-3.)*(I1*incJ-3.))*expiso;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        double Siso_i = incJ* faciso* (Id(i) - third*I1*Cinv(i));
        double Siso_j = incJ* faciso* (Id(j) - third*I1*Cinv(j));
        double Aiso_i = Id(i) - third* I1* Cinv(i);
        double Aiso_j = Id(j) - third* I1* Cinv(j);
        (*cmat)(i,j) += 2*third*incJ*faciso* I1 * Psl(i,j)
             - 2*third * Cinv(i) * Siso_j         // -2/3 Cinv x Siso
             - 2*third * Cinv(j) * Siso_i         // -2/3 Siso x Cinv
             + delta7iso * Aiso_i * Aiso_j;       // part with 4 d^2W/dC^2
      }
    }
  } else {
    // Elasticity fiber part in splitted formulation, see Holzapfel p. 255 and 272
    const double delta7bar1 = fib1_tension* 4.*(k1*exp1 + 2.*k1*k2*(J4-1.)*(J4-1.)*exp1); // 4 d^2Wf/dJ4dJ4
    const double delta7bar2 = fib2_tension* 4.*(k1*exp2 + 2.*k1*k2*(J6-1.)*(J6-1.)*exp2); // 4 d^2Wf/dJ6dJ6

    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> Caniso_fib1; // isochoric elastic C from Fib1
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> Caniso_fib2; // isochoric elastic C from Fib2

    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        double A1iso_i = incJ*A1(i)-third*J4*Cinv(i);  // A1iso = J^{-2/3} A1 - 1/3 J4 Cinv
        double A1iso_j = incJ*A1(j)-third*J4*Cinv(j);  //       = J^(-2/3)*PP:A1 has no physical meaning
        double A2iso_i = incJ*A2(i)-third*J6*Cinv(i);  // A2iso = J^{-2/3} A2 - 1/3 J6 Cinv
        double A2iso_j = incJ*A2(j)-third*J6*Cinv(j);
        Caniso_fib1(i,j) = delta7bar1 * A1iso_i * A1iso_j  // delta7bar1 A1iso x A1iso
                          + 2.*third*incJ*traceCSfbar1 * Psl(i,j)  // 2/3 J^{-2/3} trace(Sfbar C) Psl
                          - 2.*third* (Cinv(i) * Saniso_fib1(j) + Cinv(j) * Saniso_fib1(i)); // -2/3 (Cinv x Sfiso1 + Sfiso1 x Cinv)
        Caniso_fib2(i,j) = delta7bar2 * A2iso_i * A2iso_j  // delta7bar2 A2iso x A2iso
                          + 2.*third*incJ*traceCSfbar2 * Psl(i,j)  // 2/3 J^{-2/3} trace(Sfbar2 C) Psl
                          - 2.*third* (Cinv(i) * Saniso_fib2(j) + Cinv(j) * Saniso_fib2(i)); // -2/3 (Cinv x Sfiso2 + Sfiso2 x Cinv)
      }
    }

    (*cmat) += Caniso_fib1;
    (*cmat) += Caniso_fib2;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  EvaluateFiberVecs                             (public)         01/10|
 *----------------------------------------------------------------------*/
void MAT::HolzapfelCardio::EvaluateFiberVecs
(const int gp, const double gamma, const LINALG::Matrix<3,3>& locsys, const LINALG::Matrix<3,3>& defgrd)
{
  // gamma is the angle between the fibers and locsys(,2), locsys holds the principal directions
  // The deformation gradient (defgrd) is needed in remodeling as then locsys is given in the
  // spatial configuration and thus the fiber vectors have to be pulled back in the reference
  // configuration as the material is evaluated there.
  // If this function is called during Setup defgrd should be replaced by the Identity.

  for (int i = 0; i < 3; i++) {
    // a1 = cos gamma e1 + sin gamma e2 with e1 related to maximal princ stress, e2 2nd largest
    ca1_->at(gp)[i] = cos(gamma)*locsys(i,2) + sin(gamma)*locsys(i,1);
    // a2 = cos gamma e1 - sin gamma e2 with e1 related to maximal princ stress, e2 2nd largest
    ca2_->at(gp)[i] = cos(gamma)*locsys(i,2) - sin(gamma)*locsys(i,1);
  }

  // pull back in reference configuration
  std::vector<double> a1_0(3);
  std::vector<double> a2_0(3);
  LINALG::Matrix<3,3> idefgrd(false);
  idefgrd.Invert(defgrd);
  for (int i = 0; i < 3; i++) {
    a1_0[i] = idefgrd(i,0)*ca1_->at(gp)[0] + idefgrd(i,1)*ca1_->at(gp)[1] + idefgrd(i,2)*ca1_->at(gp)[2];
    a2_0[i] = idefgrd(i,0)*ca2_->at(gp)[0] + idefgrd(i,1)*ca2_->at(gp)[1] + idefgrd(i,2)*ca2_->at(gp)[2];
  }
  double a1_0norm = sqrt(a1_0[0]*a1_0[0] + a1_0[1]*a1_0[1] + a1_0[2]*a1_0[2]);
  double a2_0norm = sqrt(a2_0[0]*a2_0[0] + a2_0[1]*a2_0[1] + a2_0[2]*a2_0[2]);
  for (int i = 0; i < 3; i++) {
    a1_->at(gp)[i] = a1_0[i]/a1_0norm;
    a2_->at(gp)[i] = a2_0[i]/a2_0norm;
  }

  return;
}

