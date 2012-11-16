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


#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "Epetra_SerialDenseSolver.h"
#include "artwallremod.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io_gmsh.H"
#include "contchainnetw.H" // for debug plotting
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material_service.H"

//#define PI        (asin(1.0)*2.0)


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::ArtWallRemod::ArtWallRemod(Teuchos::RCP<MAT::PAR::Material> matdata)
  : Parameter(matdata),
    kappa_(matdata->GetDouble("KAPPA")),
    mue_(matdata->GetDouble("MUE")),
    density_(matdata->GetDouble("DENS")),
    k1_(matdata->GetDouble("K1")),
    k2_(matdata->GetDouble("K2")),
    gamma_(matdata->GetDouble("GAMMA")),
    init_(matdata->GetInt("INIT")),
    rembegt_(matdata->GetDouble("REMBEGT")),
    tensonly_(matdata->GetInt("TENSION_ONLY"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::ArtWallRemod::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ArtWallRemod(this));
}

MAT::ArtWallRemodType MAT::ArtWallRemodType::instance_;


DRT::ParObject* MAT::ArtWallRemodType::Create( const std::vector<char> & data )
{
  MAT::ArtWallRemod* remod = new MAT::ArtWallRemod();
  remod->Unpack(data);
  return remod;
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         06/08|
 *----------------------------------------------------------------------*/
MAT::ArtWallRemod::ArtWallRemod()
  : params_(NULL)
{
  dserror("This material law - ARTWALLREMOD - is not maintained anymore.");
  isinit_=false;
  gamma_ = Teuchos::rcp(new std::vector<double>);
  lambda_ = Teuchos::rcp(new std::vector<std::vector<double> >);
  a1_ = Teuchos::rcp(new std::vector<std::vector<double> >);
  a2_ = Teuchos::rcp(new std::vector<std::vector<double> >);
  phi_ = Teuchos::rcp(new std::vector<Epetra_SerialDenseMatrix>);
  stresses_ = Teuchos::rcp(new std::vector<Epetra_SerialDenseMatrix>);
  remtime_ = Teuchos::rcp(new std::vector<double>);
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          06/08|
 *----------------------------------------------------------------------*/
MAT::ArtWallRemod::ArtWallRemod(MAT::PAR::ArtWallRemod* params)
  : params_(params)
{
  dserror("This material law - ARTWALLREMOD - is not maintained anymore.");
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)         06/08|
 *----------------------------------------------------------------------*/
void MAT::ArtWallRemod::Pack(DRT::PackBuffer& data) const
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
  if (!Initialized())
  {
    numgp=0;
  }
  else
  {
    numgp = remtime_->size();   // size is number of gausspoints
  }
  AddtoPack(data,numgp);  // Length of history vector(s)
  // Pack internal variables independent of remodeling
  for (int gp = 0; gp < numgp; ++gp)
  {
    AddtoPack(data,remtime_->at(gp));
    AddtoPack(data,a1_->at(gp));
    AddtoPack(data,a2_->at(gp));
  }
  // Pack internal variables only for remodeling
  if (params_->rembegt_ != -1.){
    for (int gp = 0; gp < numgp; ++gp)
    {
      AddtoPack(data,gamma_->at(gp));
      AddtoPack(data,phi_->at(gp));
      AddtoPack(data,stresses_->at(gp));
      AddtoPack(data,lambda_->at(gp));
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)         06/08|
 *----------------------------------------------------------------------*/
void MAT::ArtWallRemod::Unpack(const std::vector<char>& data)
{
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
        params_ = static_cast<MAT::PAR::ArtWallRemod*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // history data
  isinit_ = true;
  int numgp;
  ExtractfromPack(position,data,numgp);

  if (numgp == 0) isinit_=false;

  // unpack internal variables independent of remodeling
  remtime_ = Teuchos::rcp(new std::vector<double>(numgp));
  a1_ = Teuchos::rcp(new std::vector<std::vector<double> >(numgp));
  a2_ = Teuchos::rcp(new std::vector<std::vector<double> >(numgp));
  bool haveremodeldata = false;
  for (int gp = 0; gp < numgp; ++gp) {
    double mytime;
    ExtractfromPack(position,data,mytime);
    remtime_->at(gp) = mytime;
    if (mytime != -1.) haveremodeldata = true;
    std::vector<double> a;
    ExtractfromPack(position,data,a);
    a1_->at(gp) = a;
    ExtractfromPack(position,data,a);
    a2_->at(gp) = a;
  }

  // post_drt wants to unpack but has no input variables!
  if (params_ == NULL){  // we are in post-process mode
    if (haveremodeldata){
      // read data into nowhere
      for (int gp = 0; gp < numgp; ++gp) {
        double gamma;
        Epetra_SerialDenseMatrix tmp(3,3);
        ExtractfromPack(position,data,gamma);
        ExtractfromPack(position,data,tmp);
        ExtractfromPack(position,data,tmp);
        std::vector<double> a;
        ExtractfromPack(position,data,a);
      }
    }
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d",data.size(),position);
  }
  else  // we are in analysis mode
  {
    // check whether we currently want remodeling
    // because remodeling might be switched after restart
    if (params_->rembegt_ != -1.){
      // initialize internal variables of remodeling
      gamma_ = Teuchos::rcp(new std::vector<double>(numgp));
      phi_ = Teuchos::rcp(new std::vector<Epetra_SerialDenseMatrix>(numgp));
      stresses_ = Teuchos::rcp(new std::vector<Epetra_SerialDenseMatrix>(numgp));
      lambda_ = Teuchos::rcp(new std::vector<std::vector<double> >(numgp));
      if (haveremodeldata){ // unpack remodel data
        for (int gp = 0; gp < numgp; ++gp) {
          double gamma;
          Epetra_SerialDenseMatrix tmp(3,3);
          ExtractfromPack(position,data,gamma);
          ExtractfromPack(position,data,tmp);
          gamma_->at(gp) = gamma;
          phi_->at(gp) = tmp;
          ExtractfromPack(position,data,tmp);
          stresses_->at(gp) = tmp;
          std::vector<double> a;
          ExtractfromPack(position,data,a);
          lambda_->at(gp) = a;
        }
      } else { // assign input variables or zeros, respectively
        for (int gp = 0; gp < numgp; ++gp) {
          gamma_->at(gp) = (params_->gamma_ * M_PI)/180.0;  // convert to radians
          lambda_->at(gp).resize(3);
          for (int i = 0; i < 3; ++i) lambda_->at(gp)[i] = 0.0;
          Epetra_SerialDenseMatrix emptymat(3,3);
          phi_->at(gp) = emptymat;
          stresses_->at(gp) = emptymat;
          remtime_->at(gp) = params_->rembegt_; // overwrite restart data with input when switching
        }
      }
    }

    // correct position in case remodeling is switched off
    if ((params_->rembegt_ == -1.) && (haveremodeldata)){
      // read data into nowhere
      for (int gp = 0; gp < numgp; ++gp) {
        double gamma;
        Epetra_SerialDenseMatrix tmp(3,3);
        ExtractfromPack(position,data,gamma);
        ExtractfromPack(position,data,tmp);
        ExtractfromPack(position,data,tmp);
        std::vector<double> a;
        ExtractfromPack(position,data,a);
        remtime_->at(gp) = params_->rembegt_; // overwrite restart data with input when switching
      }
    }

    // check if everything was savely unpacked
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d",data.size(),position);
  }

  // get away from this
  return;
}

void MAT::ArtWallRemod::Setup(const int numgp, const int eleid, DRT::INPUT::LineDefinition* linedef)
{
  a1_ = Teuchos::rcp(new std::vector<std::vector<double> > (numgp));
  a2_ = Teuchos::rcp(new std::vector<std::vector<double> > (numgp));
  int initflag = params_->init_;
  double gamma = params_->gamma_;
  gamma = (gamma * M_PI)/180.0;  // convert to radians
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
    std::vector<double> rad;
    std::vector<double> axi;
    std::vector<double> cir;
    // read local (cylindrical) cosy-directions at current element
    linedef->ExtractDoubleVector("RAD",rad);
    linedef->ExtractDoubleVector("AXI",axi);
    linedef->ExtractDoubleVector("CIR",cir);
    Epetra_SerialDenseMatrix locsys(3,3);
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
      EvaluateFiberVecs(gp,gamma,locsys);
    }
  } else if (initflag==2){
    dserror("Random init not yet implemented for ARTWALLREMOD");
  } else if (initflag==3){
  // start with isotropic computation, thus fiber directions are set to zero
    for (int gp = 0; gp < numgp; ++gp) {
      Epetra_SerialDenseMatrix id(3,3);
      a1_->at(gp).resize(3);
      a2_->at(gp).resize(3);
      EvaluateFiberVecs(gp,gamma,id);
    }
  } else dserror("Unknown init for ARTWALLREMOD");

  // check for remodelling option and initialize
  if ((params_->rembegt_ > 0.)){
    // history
    gamma_ = Teuchos::rcp(new std::vector<double> (numgp));  // of alignment angles
    lambda_ = Teuchos::rcp(new std::vector<std::vector<double> > (numgp)); // of eigenvalues
    phi_ = Teuchos::rcp(new std::vector<Epetra_SerialDenseMatrix> (numgp)); // of eigenvectors
    stresses_ = Teuchos::rcp(new std::vector<Epetra_SerialDenseMatrix> (numgp)); // of stresses
    for (int gp = 0; gp < numgp; ++gp) {
      gamma_->at(gp) = gamma;
      lambda_->at(gp).resize(3);
      for (int i = 0; i < 3; ++i) lambda_->at(gp)[i] = 0.0;
      Epetra_SerialDenseMatrix emptymat(3,3);
      phi_->at(gp) = emptymat;
      stresses_->at(gp) = emptymat;
    }
  } else {
    // no remodeling
    params_->rembegt_ = -1.0;
  }
  remtime_ = Teuchos::rcp(new std::vector<double> (numgp)); // of remodelling time
  for (int gp = 0; gp < numgp; ++gp) remtime_->at(gp) = params_->rembegt_;

  isinit_ = true;
  return;
}



/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)         06/08|
 *----------------------------------------------------------------------*

*/

void MAT::ArtWallRemod::Evaluate(
        const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
        const int gp,
        Teuchos::ParameterList& params,
        LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> * cmat,
        LINALG::Matrix<NUM_STRESS_3D,1> * stress,
        const LINALG::Matrix<3,3>& defgrd)

{
  const double mue = params_->mue_;
  const double kappa = params_->kappa_;
  const double k1 = params_->k1_;
  const double k2 = params_->k2_;
  const int tensonly = params_->tensonly_;

  // right Cauchy-Green Tensor  C = 2 * E + I
  // build identity tensor I
  LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
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
  const double incJ = std::pow(I3,-1.0/3.0);  // J^{-2/3}

  // invert C
  LINALG::Matrix<NUM_STRESS_3D,1> Cinv(false);

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

  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> Psl(true);
  // Psl = Cinv o Cinv - 1/3 Cinv x Cinv
  AddtoCmatHolzapfelProduct(Psl,Cinv,1.0);  // first part Psl = Cinv o Cinv

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

  // decide whether its time to remodel
  const double time = params.get("total time",-1.0);
//  const double dt = params.get("delta time",-1.0);
//  if ((remtime_->at(gp) != -1.) && (time > dt+remtime_->at(gp))){
  if ((remtime_->at(gp) != -1.) && (time > remtime_->at(gp))){
    Remodel(gp,time,defgrd);
  }

  // switch between isotropic and anisotropic fiber material for initial iteration step
  if (a1_->at(gp)[0]==0 && a1_->at(gp)[1]==0 && a1_->at(gp)[2]==0) {
    // isotropic fiber part: ***********************************************
    // W=(k1/(2.0*k2))*(exp(k2*pow((Ibar_1 - 3.0),2)-1.0)); fiber SEF
    const double expiso = exp(k2*(I1*incJ-3.)*(I1*incJ-3.));
    const double faciso = 2.*k1*(I1*incJ-3.)*expiso;
    const double delta7iso = incJ*incJ* 4.*(k1 + 2.*k1*k2*(I1*incJ-3.)*(I1*incJ-3.))*expiso;
    for (int i = 0; i < 6; ++i) {
      (*stress)(i) += incJ* faciso* (Id(i) - third*I1*Cinv(i));
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
    // anisotropic part: ***************************************************
    // W_aniso=(k1/(2.0*k2))*(exp(k2*pow((Ibar_{4,6} - 1.0),2)-1.0)); fiber SEF
    // structural tensors in voigt notation
    LINALG::Matrix<NUM_STRESS_3D,1> A1;
    LINALG::Matrix<NUM_STRESS_3D,1> A2;
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

    // 'tensonly' determines if fibers can only take tension or not
    double fib1_tension = 1.;
    double fib2_tension = 1.;
    if (tensonly==1)
    {
      if (J4 < 1.0) fib1_tension = 0.;
      if (J6 < 1.0) fib2_tension = 0.;
    }

    // PK2 fiber part in splitted formulation, see Holzapfel p. 271
    LINALG::Matrix<NUM_STRESS_3D,1> Sfiso(A1); // first compute Sfbar = dWf/dJ4 A1 + dWf/dJ6 A2
    const double fib1 = fib1_tension* 2.*(k1*(J4-1.)*exp1);  // 2 dWf/dJ4
    const double fib2 = fib2_tension* 2.*(k1*(J6-1.)*exp2);  // 2 dWf/dJ6
    Sfiso.Scale(fib1);
    Sfiso.Update(fib2,A2,1.0);

    const double traceCSfbar =  Sfiso(0)*C(0) + Sfiso(1)*C(1) + Sfiso(2)*C(2)
                   + 1.*(Sfiso(3)*C(3) + Sfiso(4)*C(4) + Sfiso(5)*C(5)); // trace(Sfbar C)
    // compute Sfiso = J^{-2/3} * (Sfbar - 1/3 trace(Sfbar C) Cinv
    for (int i = 0; i < 6; ++i) {
      Sfiso(i) = incJ * (Sfiso(i) - third*traceCSfbar*Cinv(i));
    }
    (*stress) += Sfiso;

    // Elasticity fiber part in splitted formulation, see Holzapfel p. 255 and 272
    const double delta7bar1 = fib1_tension* 4.*(k1*exp1 + 2.*k1*k2*(J4-1.)*(J4-1.)*exp1); // 4 d^2Wf/dJ4dJ4
    const double delta7bar2 = fib2_tension* 4.*(k1*exp2 + 2.*k1*k2*(J6-1.)*(J6-1.)*exp2); // 4 d^2Wf/dJ6dJ6

    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
        double A1iso_i = incJ*A1(i)-third*J4*Cinv(i);  // A1iso = J^{-2/3} A1 - 1/3 J4 Cinv
        double A1iso_j = incJ*A1(j)-third*J4*Cinv(j);
        double A2iso_i = incJ*A2(i)-third*J6*Cinv(i);  // A2iso = J^{-2/3} A2 - 1/3 J6 Cinv
        double A2iso_j = incJ*A2(j)-third*J6*Cinv(j);
        (*cmat)(i,j)    += delta7bar1 * A1iso_i * A1iso_j  // delta7bar1 A1iso x A1iso
                      + delta7bar2 * A2iso_i * A2iso_j  // delta7bar2 A2iso x A2iso
                      + 2.*third*incJ*traceCSfbar * Psl(i,j)  // 2/3 J^{-2/3} trace(Sfbar C) Psl
                      - 2.*third* (Cinv(i) * Sfiso(j) + Cinv(j) * Sfiso(i)); // -2/3 (Cinv x Sfiso + Sfiso x Cinv)
      }
    }
  } // end of switch between isotropic and anisotropic fiber strain energy function

  // store current stress in case of remodeling
  if (remtime_->at(gp) != -1.){
    for (int i = 0; i < 3; ++i) stresses_->at(gp)(i,i) = (*stress)(i);
    stresses_->at(gp)(0,1) = (*stress)(3); stresses_->at(gp)(1,0) = (*stress)(3);
    stresses_->at(gp)(1,2) = (*stress)(4); stresses_->at(gp)(2,1) = (*stress)(4);
    stresses_->at(gp)(0,2) = (*stress)(5); stresses_->at(gp)(2,0) = (*stress)(5);

    // store Cauchy stresses and use those for remodeling driver
    double detF = defgrd(0,0)*defgrd(1,1)*defgrd(2,2) +
                  defgrd(0,1)*defgrd(1,2)*defgrd(2,0) +
                  defgrd(0,2)*defgrd(1,0)*defgrd(2,1) -
                  defgrd(0,2)*defgrd(1,1)*defgrd(2,0) -
                  defgrd(0,0)*defgrd(1,2)*defgrd(2,1) -
                  defgrd(0,1)*defgrd(1,0)*defgrd(2,2);
    LINALG::Matrix<3,3> pk1(true);
    //pk1.Multiply('N','N',1.0,defgrd,stresses_->at(gp),0.);
    LINALG::DENSEFUNCTIONS::multiply<double,3,3,3>(pk1.A(),defgrd.A(),stresses_->at(gp).A());
    LINALG::DENSEFUNCTIONS::multiplyNT<double,3,3,3>(stresses_->at(gp).A(),1./detF,pk1.A(),defgrd.A());
  }

  return;
}

void MAT::ArtWallRemod::Remodel(const int gp, const double time, const LINALG::Matrix<3,3>& defgrd)
{
  // evaluate eigenproblem based on stress of previous step
  Epetra_SerialDenseVector lambda(3);
  // watch out! stress matrix will temporarily hold eigenvectors!
  LINALG::SymmetricEigenProblem(stresses_->at(gp),lambda);

#if DEBUG1
  cout << "eigenvectors: " << stresses_->at(gp);
  cout << "eigenvalues: " << lambda << endl;
#endif
  // modulation function acc. Hariton: tan g = 2nd max lambda / max lambda
  double newgamma = atan(lambda(1)/lambda(2));
  //compression in 2nd max direction, thus fibers are alligned to max principal direction
  if (lambda(1) < 0) newgamma = 0.0;

//  // check whether delta gamma is larger than tolerance
//  const double gammatol = 0.01;
//  if (abs( (newgamma - gamma_->at(gp)) / newgamma) < gammatol){
//    //remtime_->at(gp) = -1.;  // switch off future remodeling for this gp
//    remtime_->at(gp) = time; // no remodelling, but update time
//    return; // get out here
//  }

  EvaluateFiberVecs(gp,newgamma,stresses_->at(gp)); // remember! stresses holds eigenvectors

  // pull-back of new fiber vecs
  std::vector<double> a1_0(3);
  std::vector<double> a2_0(3);
  LINALG::Matrix<3,3> idefgrd(false);
  idefgrd.Invert(defgrd);
  for (int i = 0; i < 3; ++i) {
    a1_0[i] = idefgrd(i,0)*a1_->at(gp)[0] + idefgrd(i,1)*a1_->at(gp)[1] + idefgrd(i,2)*a1_->at(gp)[2];
    a2_0[i] = idefgrd(i,0)*a2_->at(gp)[0] + idefgrd(i,1)*a2_->at(gp)[1] + idefgrd(i,2)*a2_->at(gp)[2];
  }
  double a1_0norm = sqrt(a1_0[0]*a1_0[0] + a1_0[1]*a1_0[1] + a1_0[2]*a1_0[2]);
  double a2_0norm = sqrt(a2_0[0]*a2_0[0] + a2_0[1]*a2_0[1] + a2_0[2]*a2_0[2]);
  for (int i = 0; i < 3; ++i) {
    a1_->at(gp)[i] = a1_0[i]/a1_0norm;
    a2_->at(gp)[i] = a2_0[i]/a2_0norm;
  }

  // update
  gamma_->at(gp) = newgamma;
  remtime_->at(gp) = time; // we remodel only once per timestep, not during iteration

  // debug/plotting storage
  phi_->at(gp) = stresses_->at(gp);
  for (int i = 0; i < 3; ++i) lambda_->at(gp)[i] = lambda(i);

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

/// Debug Output to txt-file

/* this needs to be copied to STR::TimInt::OutputStep() to enable debug output
{
  discret_->SetState("displacement",Dis());
  MAT::ArtWallRemodOutputToGmsh(discret_, GetStep(), 1);
  MAT::ArtWallRemodOutputToTxt(discret_, GetStep(), 1);
}
don't forget to include artwallremod.H */

void MAT::ArtWallRemodOutputToTxt(const Teuchos::RCP<DRT::Discretization> dis,
    const double time,
    const int iter)
{
    std::stringstream filename;
    const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
    filename << filebase << "_rem" << ".txt";
    std::ofstream outfile;
    outfile.open(filename.str().c_str(),ios_base::app);
    //int nele = dis->NumMyColElements();
    int endele = 200; //nele;
    for (int iele=0; iele<endele; iele+=10) //++iele) iele+=10)
    {
      const DRT::Element* actele = dis->lColElement(iele);
      RCP<MAT::Material> mat = actele->Material();
      if (mat->MaterialType() != INPAR::MAT::m_artwallremod) return;
      MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(mat.get());
      //int ngp = remo->Geta1()->size();  // unused
      int endgp = 1; //ngp;
      for (int gp = 0; gp < endgp; ++gp){ //gp+=4
        double gamma = remo->Getgammas()->at(gp);
        double remtime = remo->Getremtimes()->at(gp);
        std::vector<double> lamb = remo->Getlambdas()->at(gp);
        Epetra_SerialDenseMatrix phi = remo->Getphis()->at(gp);
        std::vector<double> a1s = remo->Geta1()->at(gp);
        std::vector<double> a2s = remo->Geta2()->at(gp);

        // time
        outfile << time << ",";
        // iter
        outfile << iter << ",";
        // eleId
        outfile << iele << ",";
        // gp
        outfile << gp << ",";

        outfile << gamma << ",";
        outfile << gamma*180./PI << ",";
        outfile << remtime << ",";
        // eigenvalues
        for (int i=0;i<3;++i) outfile << lamb[i] << ",";
        // eigenvectors
        for (int i=0;i<3;++i)
          for (int j=0;j<3;++j)
            outfile << phi(j,i) << ",";
        // end
        // fiber directions
        for (int i=0;i<3;++i) outfile << a1s[i] << ",";
        for (int i=0;i<3;++i) outfile << a2s[i] << ",";
        outfile << endl;
      }
    }
    outfile.close();
    return;
}

/// Debug output to gmsh-file
void MAT::ArtWallRemodOutputToGmsh(const Teuchos::RCP<DRT::Discretization> dis,
                                      const double time,
                                      const int iter)
{
  std::stringstream filename;
  const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
  filename << filebase << "_rem" << std::setw(3) << std::setfill('0') << time << std::setw(2) << std::setfill('0') << iter << ".pos";
  std::ofstream f_system(filename.str().c_str());

  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" Time: " << time << " Iter: " << iter << " \" {" << endl;
  for (int iele=0; iele<dis->NumMyColElements(); ++iele)
  {
    const DRT::Element* actele = dis->lColElement(iele);

    // build current configuration
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    actele->LocationVector(*dis,lm,lmowner,lmstride);
    RCP<const Epetra_Vector> disp = dis->GetState("displacement");
    std::vector<double> mydisp(lm.size(),0);
    //DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
    const int numnode = actele->NumNode();
    const int numdof = 3;
    LINALG::SerialDenseMatrix xyze(3, numnode);
    for (int inode = 0; inode < numnode; ++inode)
    {
      xyze(0, inode) = actele->Nodes()[inode]->X()[0]+ mydisp[inode*numdof+0];
      xyze(1, inode) = actele->Nodes()[inode]->X()[1]+ mydisp[inode*numdof+1];
      xyze(2, inode) = actele->Nodes()[inode]->X()[2]+ mydisp[inode*numdof+2];
    }
    gmshfilecontent << IO::GMSH::cellWithScalarToString(actele->Shape(),
        1.0, xyze) << endl;

    std::vector<double> elecenter = MAT::MatPointCoords(actele,mydisp);
    RCP<MAT::Material> mat = actele->Material();
    MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(mat.get());
    RCP<vector<std::vector<double> > > a1s = remo->Geta1();
    RCP<vector<std::vector<double> > > a2s = remo->Geta2();

    // material plot at gauss points
    int ngp = remo->Geta1()->size();
    for (int gp = 0; gp < ngp; ++gp){
      std::vector<double> point = MAT::MatPointCoords(actele,mydisp,gp);

      vector<std::vector<double> > fibgp(2);
      fibgp.at(0) = a1s->at(gp);
      fibgp.at(1) = a2s->at(gp);

      for (int k=0; k<2; ++k){

        gmshfilecontent << "VP(" << std::scientific << point[0] << ",";
        gmshfilecontent << std::scientific << point[1] << ",";
        gmshfilecontent << std::scientific << point[2] << ")";
        gmshfilecontent << "{" << std::scientific
        <<        fibgp.at(k)[0]
        << "," << fibgp.at(k)[1]
        << "," << fibgp.at(k)[2]
        << "};" << endl;

        // draw also negative direction to avoid "jumping"
        gmshfilecontent << "VP(" << std::scientific << point[0] << ",";
        gmshfilecontent << std::scientific << point[1] << ",";
        gmshfilecontent << std::scientific << point[2] << ")";
        gmshfilecontent << "{" << std::scientific
        <<        -fibgp.at(k)[0]
        << "," << -fibgp.at(k)[1]
        << "," << -fibgp.at(k)[2]
        << "};" << endl;

      }

//      Epetra_SerialDenseMatrix Phi= remo->Getphis()->at(gp);
//      std::vector<double> lambda = remo->Getlambdas()->at(gp);
//      const double scale = 100.0;
//      for (int k=0; k<3; ++k){
//        gmshfilecontent << "VP(" << std::scientific << point[0] << ",";
//        gmshfilecontent << std::scientific << point[1] << ",";
//        gmshfilecontent << std::scientific << point[2] << ")";
//        gmshfilecontent << "{" << std::scientific
//        <<        scale*lambda[k]*Phi(0,k)
//        << "," << scale*lambda[k]*Phi(1,k)
//        << "," << scale*lambda[k]*Phi(2,k)
//        << "};" << endl;
//      }

    }
  }
  gmshfilecontent << "};" << endl;

  f_system << gmshfilecontent.str();
  f_system.close();

  return;
}



