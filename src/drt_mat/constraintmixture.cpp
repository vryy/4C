/*!----------------------------------------------------------------------
\file constraintmixture.cpp
\brief
This file contains routines for constraint mixture growth and remodeling.
example input line
MAT 1 MAT_ConstraintMixture DENS 0.001 MUE 1.0E3 PHIE 0.08 PREELA 1.0
K1 1.0 K2 1.0 PRECOLL 1.0 KAPPA 1.0E4 LIFETIME 5.0 HOMSTR 0.0
GROWTHFAC 0.0 STARTTIME 100.0 INTEGRATION Explicit TOL 1.0E-8
GROWTHFORCE Single

Here an approach for growth and remodeling of an artery is modeled.
For a detailed description see:
- Humphrey, J. & Rajagopal, K.: A constrained mixture model for arterial
  adaptations to a sustained step change in blood flow,
  Biomechanics and Modeling in Mechanobiology, 2003, 2, 109-126

<pre>
Maintainer: Susanna Tinkl
            tinkl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <vector>
#include "constraintmixture.H"
#include "../drt_lib/drt_globalproblem.H"
#include "matpar_bundle.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_io/io_gmsh.H" // for debug plotting with gmsh
#include "../drt_io/io_control.H" // for debug plotting with gmsh
#include "contchainnetw.H" // for debug plotting with gmsh

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::ConstraintMixture::ConstraintMixture(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  density_(matdata->GetDouble("DENS")),
  mue_(matdata->GetDouble("MUE")),
  phielastin_(matdata->GetDouble("PHIE")),
  prestretchelastin_(matdata->GetDouble("PREELA")),
  k1_(matdata->GetDouble("K1")),
  k2_(matdata->GetDouble("K2")),
  prestretchcollagen_(matdata->GetDouble("PRECOLL")),
  kappa_(matdata->GetDouble("KAPPA")),
  lifetime_(matdata->GetDouble("LIFETIME")),
  homstress_(matdata->GetDouble("HOMSTR")),
  growthfactor_(matdata->GetDouble("GROWTHFAC")),
  starttime_(matdata->GetDouble("STARTTIME")),
  integration_(matdata->Get<string>("INTEGRATION")),
  abstol_(matdata->GetDouble("TOL")),
  growthforce_(matdata->Get<string>("GROWTHFORCE"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::ConstraintMixture::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ConstraintMixture(this));
}


MAT::ConstraintMixtureType MAT::ConstraintMixtureType::instance_;

DRT::ParObject* MAT::ConstraintMixtureType::Create( const std::vector<char> & data )
{
  MAT::ConstraintMixture* comix = new MAT::ConstraintMixture();
  comix->Unpack(data);
  return comix;
}

/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         12/10|
 *----------------------------------------------------------------------*/
MAT::ConstraintMixture::ConstraintMixture()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          12/10|
 *----------------------------------------------------------------------*/
MAT::ConstraintMixture::ConstraintMixture(MAT::PAR::ConstraintMixture* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*
 |  Pack                                          (public)         12/10|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::Pack(DRT::PackBuffer& data) const
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
  } else {
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
    AddtoPack(data,vismassstress_->at(gp));
    AddtoPack(data,refmassdens_->at(gp));
  }
  if (numgp > 0)
  {
    AddtoPack(data,massprodbasal_);

    // Pack history
    int sizehistory = history_->size();
    AddtoPack(data,sizehistory);
    for (int idpast = 0; idpast < sizehistory; idpast++)
      history_->at(idpast).Pack(data);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)         12/10|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::Unpack(const vector<char>& data)
{
  isinit_=true;
  vector<char>::size_type position = 0;
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
        params_ = static_cast<MAT::PAR::ConstraintMixture*>(mat);
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
  a1_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> >(numgp));
  a2_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> >(numgp));
  a3_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> >(numgp));
  a4_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> >(numgp));
  vismassstress_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> > (numgp));
  refmassdens_ = Teuchos::rcp(new vector<double> (numgp));

  for (int gp = 0; gp < numgp; ++gp) {
    LINALG::Matrix<3,1> alin;
    ExtractfromPack(position,data,alin);
    a1_->at(gp) = alin;
    ExtractfromPack(position,data,alin);
    a2_->at(gp) = alin;
    ExtractfromPack(position,data,alin);
    a3_->at(gp) = alin;
    ExtractfromPack(position,data,alin);
    a4_->at(gp) = alin;
    ExtractfromPack(position,data,alin);
    vismassstress_->at(gp) = alin;
    double a;
    ExtractfromPack(position,data,a);
    refmassdens_->at(gp) = a;
  }
  double basal;
  ExtractfromPack(position,data,basal);
  massprodbasal_=basal;

  // unpack history
  int sizehistory;
  ExtractfromPack(position,data,sizehistory);
  history_ = Teuchos::rcp(new vector<ConstraintMixtureHistory> (sizehistory));
  vector<char> datahistory;
  for (int idpast = 0; idpast < sizehistory; idpast++)
  {
    ExtractfromPack(position,data,datahistory);
    history_->at(idpast).Unpack(datahistory);
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 |  Setup                                         (public)         12/10|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::Setup (const int numgp, DRT::INPUT::LineDefinition* linedef)
{
  if (*params_->integration_ != "Implicit" && *params_->integration_ != "Explicit")
    dserror("unknown option for integration");
  if (*params_->growthforce_ != "Single" && *params_->growthforce_ != "All")
    dserror("unknown driving force for growth");

  // visualization
  vismassstress_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> > (numgp));
  refmassdens_ = Teuchos::rcp(new vector<double> (numgp));
  for (int gp = 0; gp < numgp; gp++)
  {
    vismassstress_->at(gp)(0) = 0.0;
    vismassstress_->at(gp)(1) = 0.0;
    vismassstress_->at(gp)(2) = 0.0;
    refmassdens_->at(gp) = params_->density_;
  }

  // basal mass production rate determined by DENS, PHIE and LIFETIME
  massprodbasal_ = (1 - params_->phielastin_) * params_->density_ / 4.0 / params_->lifetime_;

  // history
  SetupHistory(numgp);

  // fiber vectors
  a1_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> > (numgp));
  a2_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> > (numgp));
  a3_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> > (numgp));
  a4_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> > (numgp));

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

  const double gamma = (45*PI)/180.; //angle for diagonal fibers

  for (int gp = 0; gp < numgp; gp++)
  {
    for (int i = 0; i < 3; i++) {
      // a1 = e3, circumferential direction
      a1_->at(gp)(i) = locsys(i,2);
      // a2 = e2, axial direction
      a2_->at(gp)(i) = locsys(i,1);
      // a3 = cos gamma e3 + sin gamma e2
      a3_->at(gp)(i) = cos(gamma)*locsys(i,2) + sin(gamma)*locsys(i,1);
      // a4 = cos gamma e3 - sin gamma e2
      a4_->at(gp)(i) = cos(gamma)*locsys(i,2) - sin(gamma)*locsys(i,1);
    }
  }

  isinit_ = true;
}

/*----------------------------------------------------------------------*
 |  SetupHistory                                  (public)         03/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::SetupHistory (const int numgp)
{
  const Teuchos::ParameterList& timeintegr = DRT::Problem::Instance()->StructuralDynamicParams();
  double dt = timeintegr.get<double>("TIMESTEP");
  int numint = 0;
  if (*params_->integration_ == "Explicit")
    numint = 1;
  int numpast = 0;
  if (abs(round(params_->lifetime_ / dt) - params_->lifetime_ / dt) < 1.0e-8){
    numpast = static_cast<int>(round(params_->lifetime_ / dt)) + numint;
  } else {
    numpast = static_cast<int>(ceil(params_->lifetime_ / dt)) + numint;
  }
  // history
  history_ = Teuchos::rcp(new vector<ConstraintMixtureHistory> (numpast));
  for (int idpast = 0; idpast < numpast; idpast++)
  {
    history_->at(idpast).Setup(numgp,massprodbasal_);
    history_->at(idpast).SetTime(dt-(numpast-1-idpast)*dt,dt);
  }
}

/*----------------------------------------------------------------------*
 |  Update internal variables                     (public)         12/10|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::Update()
{
  // erase oldest collagen
  int numgp = history_->at(0).NumGP();
  int sizehistory = history_->size();
  double deptime = 0.0;
  double depdt = 0.0;
  double eps = 1.0e-12;
  // delete just the steps that surely won't be needed, especially with a smaller timestep later
  // thus reference time is deposition time of last collagen fibers
  history_->back().GetTime(&deptime,&depdt);
  double erasetime = deptime - params_->lifetime_;
  int eraseiter = 0;
  history_->at(eraseiter).GetTime(&deptime,&depdt);
  while (deptime < erasetime-eps && eraseiter < sizehistory)
  {
    eraseiter +=1;
    history_->at(eraseiter).GetTime(&deptime,&depdt);
  }
  if (eraseiter > 0)
  {
    history_->erase(history_->begin(),history_->begin()+eraseiter);
    //cout << "erased " << eraseiter << " history variables" << endl;
  }

  // append new collagen
  ConstraintMixtureHistory newhis;
  newhis.Setup(numgp,massprodbasal_);
  // it is very important to set time and dt to 0.0
  // this makes it clear that this step was created in update and has no reliable content
  // they are not known here either
  // in EvaluateStress this is important, if called after Update
  newhis.SetTime(0.0,0.0);
  history_->push_back(newhis);
}

/*----------------------------------------------------------------------*
 |  Evaluate                                      (public)         12/10|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::Evaluate
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  const int gp,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress,
  double dt,
  double time,
  bool output
)
{
  // differentiate between explicit and implicit description
  int firstiter = 0;
  if (*params_->integration_ == "Explicit")
    firstiter = 1;

  if (!output) {
    // set actual time as it might have changed after an restart etc. but just once
    // in remodeling time might be wrong depending on the time integration used
    // correct this for the computation but do not store it
    double temptime = 0.0;
    double tempdt = 0.0;
    history_->back().GetTime(&temptime,&tempdt);
    if (temptime == 0.0 && tempdt == 0.0)
      history_->back().SetTime(time,dt);
    else if (time > temptime)
      time = temptime;

    //--------------------------------------------------------------------------------------
    // build identity tensor I
    LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
    for (int i = 0; i < 3; i++) Id(i) = 1.0;
    // right Cauchy-Green Tensor  C = 2 * E + I
    LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
    C.Scale(2.0);
    C += Id;

    //--------------------------------------------------------------------------------------
    // compute actual collagen stretches
    LINALG::Matrix<4,1> actstretch(true);
    double actcollagenstretch =  a1_->at(gp)(0)*a1_->at(gp)(0)*C(0) + a1_->at(gp)(1)*a1_->at(gp)(1)*C(1)
          + a1_->at(gp)(2)*a1_->at(gp)(2)*C(2) + a1_->at(gp)(0)*a1_->at(gp)(1)*C(3)
          + a1_->at(gp)(1)*a1_->at(gp)(2)*C(4) + a1_->at(gp)(0)*a1_->at(gp)(2)*C(5); // = trace(A1:C)
    actcollagenstretch = sqrt(actcollagenstretch);
    actstretch(0) = actcollagenstretch;
    actcollagenstretch =  a2_->at(gp)(0)*a2_->at(gp)(0)*C(0) + a2_->at(gp)(1)*a2_->at(gp)(1)*C(1)
          + a2_->at(gp)(2)*a2_->at(gp)(2)*C(2) + a2_->at(gp)(0)*a2_->at(gp)(1)*C(3)
          + a2_->at(gp)(1)*a2_->at(gp)(2)*C(4) + a2_->at(gp)(0)*a2_->at(gp)(2)*C(5); // = trace(A2:C)
    actcollagenstretch = sqrt(actcollagenstretch);
    actstretch(1) = actcollagenstretch;
    actcollagenstretch =  a3_->at(gp)(0)*a3_->at(gp)(0)*C(0) + a3_->at(gp)(1)*a3_->at(gp)(1)*C(1)
          + a3_->at(gp)(2)*a3_->at(gp)(2)*C(2) + a3_->at(gp)(0)*a3_->at(gp)(1)*C(3)
          + a3_->at(gp)(1)*a3_->at(gp)(2)*C(4) + a3_->at(gp)(0)*a3_->at(gp)(2)*C(5); // = trace(A3:C)
    actcollagenstretch = sqrt(actcollagenstretch);
    actstretch(2) = actcollagenstretch;
    actcollagenstretch =  a4_->at(gp)(0)*a4_->at(gp)(0)*C(0) + a4_->at(gp)(1)*a4_->at(gp)(1)*C(1)
          + a4_->at(gp)(2)*a4_->at(gp)(2)*C(2) + a4_->at(gp)(0)*a4_->at(gp)(1)*C(3)
          + a4_->at(gp)(1)*a4_->at(gp)(2)*C(4) + a4_->at(gp)(0)*a4_->at(gp)(2)*C(5); // = trace(A4:C)
    actcollagenstretch = sqrt(actcollagenstretch);
    actstretch(3) = actcollagenstretch;

    // store them
    // time <= starttime needs further considerations
    double eps = 1.0e-12;
    if (time > params_->starttime_ + eps)
    {
      history_->back().SetStretches(gp,actstretch);
    }
//    } else {  // this is not working for all material parameters
//      int numsteps = history_->size();
//      for (int i = 0; i < numsteps; i++)
//        history_->at(i).SetStretches(gp,actstretch);
//    }

    EvaluateStress(glstrain, gp, cmat, stress, firstiter, dt, time);

    //--------------------------------------------------------------------------------------
    // compute new deposition rates
    // either for future use or just for visualization
    LINALG::Matrix<4,1> massstress(true);
    LINALG::Matrix<4,1> massprodcomp(true);
    if (*params_->growthforce_ == "All")
      MassProduction(gp, C, *stress, &massstress, &massprodcomp);
    else {
      double massstresstemp = 0.0;
      double massprodtemp = 0.0;
      LINALG::Matrix<NUM_STRESS_3D,1> stresstemp(true);
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmattemp(true);
      double masstemp = 0.0;
      EvaluateFiberFamily(glstrain, gp, &cmattemp, &stresstemp, a1_->at(gp), &masstemp, firstiter, dt, time, 0);
      MassProductionSingleFiber(gp, C, stresstemp, &massstresstemp, &massprodtemp, a1_->at(gp), 0);
      massstress(0) = massstresstemp;
      massprodcomp(0) = massprodtemp;
      stresstemp.Scale(0.0);
      EvaluateFiberFamily(glstrain,gp,&cmattemp,&stresstemp,a2_->at(gp),&masstemp,firstiter,dt,time,1);
      MassProductionSingleFiber(gp, C, stresstemp, &massstresstemp, &massprodtemp, a2_->at(gp), 1);
      massstress(1) = massstresstemp;
      massprodcomp(1) = massprodtemp;
      stresstemp.Scale(0.0);
      EvaluateFiberFamily(glstrain,gp,&cmattemp,&stresstemp,a3_->at(gp),&masstemp,firstiter,dt,time,2);
      MassProductionSingleFiber(gp, C, stresstemp, &massstresstemp, &massprodtemp, a3_->at(gp), 2);
      massstress(2) = massstresstemp;
      massprodcomp(2) = massprodtemp;
      stresstemp.Scale(0.0);
      EvaluateFiberFamily(glstrain,gp,&cmattemp,&stresstemp,a4_->at(gp),&masstemp,firstiter,dt,time,3);
      MassProductionSingleFiber(gp, C, stresstemp, &massstresstemp, &massprodtemp, a4_->at(gp), 3);
      massstress(3) = massstresstemp;
      massprodcomp(3) = massprodtemp;
    }

    if (time > params_->starttime_ + eps && params_->growthfactor_ != 0.0)
    {
      // start values for local Newton iteration are computed
      // distinguish between explicit and implicit integration
      if (*params_->integration_ == "Explicit")
      {
        history_->back().SetMass(gp,massprodcomp);
        vismassstress_->at(gp)(0) = massstress(0);
        vismassstress_->at(gp)(1) = massstress(1);
        vismassstress_->at(gp)(2) = massstress(2);
      } else {
        if (*params_->growthforce_ == "All")
          EvaluateImplicitAll(glstrain, gp, cmat, stress, dt, time, massprodcomp, massstress);
        else
          EvaluateImplicitSingle(glstrain, gp, cmat, stress, dt, time, massprodcomp, massstress);
      }

    } else {
      // Visualization of massstresss also for time <= starttime
      if (time > params_->starttime_ + eps)
      {
        vismassstress_->at(gp)(0) = massstress(0);
        vismassstress_->at(gp)(1) = massstress(1);
        vismassstress_->at(gp)(2) = massstress(2);
      } else {
        vismassstress_->at(gp)(0) = massstress(0)/actstretch(0);
        vismassstress_->at(gp)(1) = massstress(1)/actstretch(1);
        vismassstress_->at(gp)(2) = massstress(2)/actstretch(2);//actstretch(0);
      // current stretch in Massproduction is 1.0, thus we have to consider it here
      // this scaling is not needed if the actual stretch is stored for t <= starttime
      }
    }
  } else {
    // in case of output everything is fully converged, we just have to evaluate stress etc.
    // should be independent of order of update and output, as new steps are set with dt = 0.0
    // and oldest fibers are carefully erased
    EvaluateStress(glstrain, gp, cmat, stress, firstiter, dt, time);
  }
} // Evaluate

/*----------------------------------------------------------------------*
 |  EvaluateStress                                (private)        01/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::EvaluateStress
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  const int gp,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress,
  const int firstiter,
  double dt,
  double time
)
{
  //--------------------------------------------------------------------------------------
  // some variables
  double prestretchelastin = params_->prestretchelastin_;
  double density = params_->density_;
  double currmassdens = 0.0;

  //--------------------------------------------------------------------------------------
  // calculate stress and elasticity matrix

  // 1st step: elastin
  //==========================
  double refmassdenselastin = params_->phielastin_ * density;
  if (refmassdenselastin < 0.0 || refmassdenselastin > density) dserror("mass fraction of elastin not in [0;1]");
  // account for isotropic prestretch of elastin
  LINALG::Matrix<NUM_STRESS_3D,1> glstrainiso(*glstrain);
  glstrainiso.Scale(pow(prestretchelastin,2));
  LINALG::Matrix<NUM_STRESS_3D,1> Siso(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatiso(true);
  EvaluateElastin(& glstrainiso, & cmatiso, & Siso);
  Siso.Scale(refmassdenselastin/density*pow(prestretchelastin,2));
  (*stress) = Siso;
  cmatiso.Scale(refmassdenselastin/density*pow(prestretchelastin,4));
  (*cmat) = cmatiso;
  currmassdens += refmassdenselastin;

  // 2nd step: collagen
  //==========================
  EvaluateFiberFamily(glstrain,gp,cmat,stress,a1_->at(gp),&currmassdens,firstiter,dt,time,0);

  EvaluateFiberFamily(glstrain,gp,cmat,stress,a2_->at(gp),&currmassdens,firstiter,dt,time,1);

  EvaluateFiberFamily(glstrain,gp,cmat,stress,a3_->at(gp),&currmassdens,firstiter,dt,time,2);

  EvaluateFiberFamily(glstrain,gp,cmat,stress,a4_->at(gp),&currmassdens,firstiter,dt,time,3);

  // 3rd step: volumetric part
  //==========================
  LINALG::Matrix<NUM_STRESS_3D,1> Svol(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatvol(true);
  EvaluateVolumetric(glstrain, & cmatvol, & Svol, currmassdens, density);
  (*stress) += Svol;
  (*cmat) += cmatvol;

  // set actual mass density
  refmassdens_->at(gp) = currmassdens;
}

/*----------------------------------------------------------------------*
 |  EvaluateFiberFamily                           (private)        05/11|
 *----------------------------------------------------------------------*
 strain energy function

 W^k = M^k(0)/rho*Q^k(0)*Psi^k + \int_0^t m^k(tau)/rho*q^k(t-tau)*Psi^k dtau

 */
void MAT::ConstraintMixture::EvaluateFiberFamily
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  const int gp,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress,
  LINALG::Matrix<3,1> a,
  double* currmassdens,
  const int firstiter,
  double dt,
  double time,
  const int idfiber
 )
{
  //--------------------------------------------------------------------------------------
  // some variables
  double prestretchcollagen = params_->prestretchcollagen_;
//  if (time < params_->starttime_/5) //could be problematic in output!!
//    prestretchcollagen = prestretchcollagen / 2; //* time / (params_->starttime_ / 5);
  double density = params_->density_;
  int sizehistory = history_->size();
  double eps = 1.0e-12;

  //--------------------------------------------------------------------------------------
  // prestress time
//  const ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
//  INPAR::STR::PreStress pstype = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(pslist,"PRESTRESS");
//  double pstime = -1.0 * params_->lifetime_ - dt;
//  if (pstype == INPAR::STR::prestress_mulf)
//    pstime = pslist.get<double>("PRESTRESSTIME");

  //--------------------------------------------------------------------------------------
  // calculate stress and elasticity matrix
  for (int idpast = 0; idpast < sizehistory - firstiter; idpast++)
  {
    double deptime = 0.0;
    double depdt = 0.0;
    history_->at(idpast).GetTime(&deptime,&depdt);
    // right dt for explicit integration & adaptation for special implicit step
    if (firstiter == 1)
    {
      double timeloc = 0.0;
      double dtloc = 0.0;
      history_->at(idpast+1).GetTime(&timeloc,&dtloc);
      depdt = dtloc;
    } else if (deptime >= (time - params_->lifetime_ - eps) &&
               (deptime - depdt) < (time - params_->lifetime_ - eps))
    {
      depdt = deptime - time + params_->lifetime_;
    }
    LINALG::Matrix<4,1> collstretch(true);
    history_->at(idpast).GetStretches(gp,&collstretch);
    double stretch = prestretchcollagen / collstretch(idfiber);
    // we already have prestretch in prestress time, thus the prestretch from the inputfile is not applied
//    if (deptime <= pstime + 1.0e-12)
//      stretch = 1.0 / collstretch(idfiber);
    LINALG::Matrix<NUM_STRESS_3D,1> Saniso(true);
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmataniso(true);
    EvaluateSingleFiber(glstrain, & cmataniso, & Saniso, a, stretch);
    double qdegrad = 0.0;
    Degradation(time - deptime, &qdegrad);
    LINALG::Matrix<4,1> collmass(true);
    history_->at(idpast).GetMass(gp,&collmass);
    double facS = pow(stretch,2) * qdegrad * collmass(idfiber) / density * depdt;
    Saniso.Scale(facS);
    (*stress) += Saniso;
    double faccmat = pow(stretch,4) * qdegrad * collmass(idfiber) / density * depdt;
    cmataniso.Scale(faccmat);
    // these lines are needed if the actual stretch is stored for time <= starttime
//    if (time <= params_->starttime_ + eps){
//      LINALG::Matrix<NUM_STRESS_3D,1>  A;
//      for (int i = 0; i < 3; i++) {
//        A(i) = a(i)*a(i);
//      }
//      A(3) = a(0)*a(1); A(4) = a(1)*a(2); A(5) = a(0)*a(2);
//      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatanisoadd(true);
//      cmatanisoadd.MultiplyNT(Saniso,A);
//      cmatanisoadd.Scale(-2.0/pow(collstretch(idfiber),2));
//      cmataniso = cmatanisoadd;
//    }
    (*cmat) += cmataniso;
    (*currmassdens) += qdegrad * collmass(idfiber) * depdt;
  }
}

/*----------------------------------------------------------------------*
 |  EvaluateSingleFiber                           (private)        12/10|
 *----------------------------------------------------------------------*
 strain energy function

 Psi    = k1/(2.0*k2)*(exp(k2*(I_4 - 1.0)^2)-1.0)

 I_4 .. invariant accounting for the fiber direction
 */
void MAT::ConstraintMixture::EvaluateSingleFiber
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress,
  LINALG::Matrix<3,1> a,
  double stretch
)
{
  const double k1 = params_->k1_;
  const double k2 = params_->k2_;

  //--------------------------------------------------------------------------------------
  // build identity tensor I
  LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
  C.Scale(2.0);
  C += Id;

  //--------------------------------------------------------------------------------------
  // structural tensors in voigt notation
  // A = a x a
  LINALG::Matrix<NUM_STRESS_3D,1>  A;
  for (int i = 0; i < 3; i++) {
    A(i) = a(i)*a(i);
  }
  A(3) = a(0)*a(1); A(4) = a(1)*a(2); A(5) = a(0)*a(2);

  double I4 =  A(0)*C(0) + A(1)*C(1) + A(2)*C(2)
             + 1.*(A(3)*C(3) + A(4)*C(4) + A(5)*C(5)); // I4 = trace(A C)
  I4 = I4 * pow(stretch,2);  // account for prestretch and stretch at deposition time

  //--------------------------------------------------------------------------------------
  // fibers can only stretch/compress down to a minimal value
  // look this up, I am not sure. There is nothing stated explicitly.
  // take material parameters of muscle cells
  double fib1_tension = 1.0;
  if (I4 < 1.0)
  {
    I4 = 1.0;
    fib1_tension = 0.0;
  }

  //--- determine 2nd Piola Kirchhoff stresses S -----------------------------------------
  LINALG::Matrix<NUM_STRESS_3D,1> Saniso(A); // first compute S = 2 dW/dI4 A
  const double exp1 = exp(k2*(I4-1.)*(I4-1.));
  if (isinf(exp1))
    dserror("stretch in fiber direction is too high");
  const double fib1 = 2.*(k1*(I4-1.)*exp1);  // 2 dW/dI4
  Saniso.Scale(fib1);  //S

  *stress = Saniso;

  //--- do elasticity matrix -------------------------------------------------------------
  if (cmat != NULL)
  {
    const double delta7 = fib1_tension* 4.0*(k1*exp1 + 2.0*k1*k2*(I4-1.0)*(I4-1.0)*exp1); // 4 d^2Wf/dI4dI4
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmataniso; // isochoric elastic C from Fib
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
          cmataniso(i,j) = delta7 * A(i) * A(j);  // delta7 A x A
      }
    }
    *cmat = cmataniso;
  }
}

/*----------------------------------------------------------------------*
 |  EvaluateElastin                               (private)         12/10|
 *----------------------------------------------------------------------*
 strain energy function

 W    = 1/2 mue (I1-3)

*/
void MAT::ConstraintMixture::EvaluateElastin
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress
)
{
  double mue = params_->mue_;
  LINALG::Matrix<NUM_STRESS_3D,1> Siso(true);
  for (int i = 0; i < 3; i++) Siso(i) = mue;

  *stress = Siso;
}

/*----------------------------------------------------------------------*
 |  EvaluateVolumetric                            (private)        12/10|
 *----------------------------------------------------------------------*
 strain energy function

 W    = 1/2 kappa (J-M(t)/M(0))^2

*/
void MAT::ConstraintMixture::EvaluateVolumetric
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress,
  double currMassDens,
  double refMassDens
)
{
  double kappa = params_->kappa_;

  //--------------------------------------------------------------------------------------
  // build identity tensor I
  LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
  C.Scale(2.0);
  C += Id;

  // isotropic invariant
  const double I3 = C(0)*C(1)*C(2)
        + 0.25 * C(3)*C(4)*C(5)
        - 0.25 * C(1)*C(5)*C(5)
        - 0.25 * C(2)*C(3)*C(3)
        - 0.25 * C(0)*C(4)*C(4);    // 3rd invariant, determinant
  if (I3 < 0.0)
    dserror("fatal failure in constraint mixture artery wall material");
  const double J = sqrt(I3);     // determinant of F
  const double p = kappa*(J- currMassDens/refMassDens); // dW_vol/dJ

  //--------------------------------------------------------------------------------------
  // invert C
  LINALG::Matrix<NUM_STRESS_3D,1> Cinv(true);
  Cinv(0) = C(1)*C(2) - 0.25*C(4)*C(4);
  Cinv(1) = C(0)*C(2) - 0.25*C(5)*C(5);
  Cinv(2) = C(0)*C(1) - 0.25*C(3)*C(3);
  Cinv(3) = 0.25*C(5)*C(4) - 0.5*C(3)*C(2);
  Cinv(4) = 0.25*C(3)*C(5) - 0.5*C(0)*C(4);
  Cinv(5) = 0.25*C(3)*C(4) - 0.5*C(5)*C(1);
  Cinv.Scale(1.0/I3);

  //--- determine 2nd Piola Kirchhoff stresses S -----------------------------------------
  LINALG::Matrix<NUM_STRESS_3D,1> Svol(true);
  // volumetric part J*kappa*(J-M(t)/M(0))*Cinv
  for (int i = 0; i < 6; i++) Svol(i) = J*p * Cinv(i);

  *stress = Svol;

  //--- do elasticity matrix -------------------------------------------------------------
  // cmatvol = J(p + J dp/dJ) Cinv x Cinv  -  2 J p Cinv o Cinv
  AddtoCmatHolzapfelProduct((*cmat),Cinv,(-2*J*p));
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 6; j++) {
      (*cmat)(i,j) += J*(p+J*kappa) * Cinv(i) * Cinv(j);
    }
  }
}

/*----------------------------------------------------------------------*
 |  MassProduction                                (private)        01/11|
 *----------------------------------------------------------------------*
compute new deposition rates for all fiber families
driving force depends on input, may be S or S^k
I derived a strange form for the driving factor just with S and C:

m^k = mbasal * (1 + K * ( sqrt(a_i^T*C*S*C*S*C*a_i/detC)/lambda_k(t) /homeo - 1))

therefore we need C and S as matrices
*/
void MAT::ConstraintMixture::MassProduction
(
  const int gp,
  LINALG::Matrix<NUM_STRESS_3D,1> C,
  LINALG::Matrix<NUM_STRESS_3D,1> S,
  LINALG::Matrix<4,1>* massstress,
  LINALG::Matrix<4,1>* massprodcomp
)
{
  LINALG::Matrix<3,3> Smatrix(true);
  Smatrix(0,0) = S(0);
  Smatrix(0,1) = S(3);
  Smatrix(0,2) = S(5);
  Smatrix(1,0) = Smatrix(0,1);
  Smatrix(1,1) = S(1);
  Smatrix(1,2) = S(4);
  Smatrix(2,0) = Smatrix(0,2);
  Smatrix(2,1) = Smatrix(1,2);
  Smatrix(2,2) = S(2);
  LINALG::Matrix<3,3> Cmatrix(true);
  Cmatrix(0,0) = C(0);
  Cmatrix(0,1) = 0.5 * C(3);
  Cmatrix(0,2) = 0.5 * C(5);
  Cmatrix(1,0) = Cmatrix(0,1);
  Cmatrix(1,1) = C(1);
  Cmatrix(1,2) = 0.5 * C(4);
  Cmatrix(2,0) = Cmatrix(0,2);
  Cmatrix(2,1) = Cmatrix(1,2);
  Cmatrix(2,2) = C(2);
  double detC = C(0)*C(1)*C(2) + 0.25 * C(3)*C(4)*C(5) - 0.25 * C(1)*C(5)*C(5)
        - 0.25 * C(2)*C(3)*C(3) - 0.25 * C(0)*C(4)*C(4);
  LINALG::Matrix<3,3> temp1(true);
  LINALG::Matrix<3,3> temp2(true);
  temp1.Multiply(Smatrix,Cmatrix);
  temp2.Multiply(Cmatrix,temp1);
  temp1.Scale(0.0);
  temp1.Multiply(Smatrix,temp2);
  temp2.Scale(0.0);
  temp2.Multiply(1.0/detC,Cmatrix,temp1);

  LINALG::Matrix<4,1> currentstretch(true);
  history_->back().GetStretches(gp,&currentstretch);

  double massstress1 = 0.0;
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      massstress1 += a1_->at(gp)(i) * temp2(i,j) * a1_->at(gp)(j);
    }
  }
  massstress1 = sqrt(massstress1) / currentstretch(0);
  if (params_->homstress_ != 0.0)
    (*massprodcomp)(0) = massprodbasal_ * (1.0 + params_->growthfactor_ * (massstress1 / params_->homstress_ - 1.0));
  else
    (*massprodcomp)(0) = massprodbasal_ * (1.0 + params_->growthfactor_ * massstress1);

  double massstress2 = 0.0;
  for (int i=0; i < 3; i++){
    for (int j=0; j < 3; j++){
      massstress2 += a2_->at(gp)(i) * temp2(i,j) * a2_->at(gp)(j);
    }
  }
  massstress2 = sqrt(massstress2) / currentstretch(1);
  if (params_->homstress_ != 0.0)
    (*massprodcomp)(1) = massprodbasal_ * (1.0 + params_->growthfactor_ * (massstress2 / params_->homstress_ - 1.0));
  else
    (*massprodcomp)(1) = massprodbasal_ * (1.0 + params_->growthfactor_ * massstress2);

  double massstress3 = 0.0;
  for (int i=0; i < 3; i++){
    for (int j=0; j < 3; j++){
      massstress3 += a3_->at(gp)(i) * temp2(i,j) * a3_->at(gp)(j);
    }
  }
  massstress3 = sqrt(massstress3) / currentstretch(2);
  if (params_->homstress_ != 0.0)
    (*massprodcomp)(2) = massprodbasal_ * (1.0 + params_->growthfactor_ * (massstress3 / params_->homstress_ - 1.0));
  else
    (*massprodcomp)(2) = massprodbasal_ * (1.0 + params_->growthfactor_ * massstress3);

  double massstress4 = 0.0;
  for (int i=0; i < 3; i++){
    for (int j=0; j < 3; j++){
      massstress4 += a4_->at(gp)(i) * temp2(i,j) * a4_->at(gp)(j);
    }
  }
  massstress4 = sqrt(massstress4) / currentstretch(3);
  if (params_->homstress_ != 0.0)
    (*massprodcomp)(3) = massprodbasal_ * (1.0 + params_->growthfactor_ * (massstress4 / params_->homstress_ - 1.0));
  else
    (*massprodcomp)(3) = massprodbasal_ * (1.0 + params_->growthfactor_ * massstress4);

  (*massstress)(0) = massstress1;
  (*massstress)(1) = massstress2;
  (*massstress)(2) = massstress3;
  (*massstress)(3) = massstress4;

//  if (((*massprodcomp)(0) < 0) || ((*massprodcomp)(1) < 0) ||
//      ((*massprodcomp)(2) < 0) || ((*massprodcomp)(3) < 0))
//  {
//    cout << "1: " << (*massprodcomp)(0) << endl;
//    cout << "2: " << (*massprodcomp)(1) << endl;
//    cout << "3: " << (*massprodcomp)(2) << endl;
//    cout << "4: " << (*massprodcomp)(3) << endl;
//    dserror("negative mass production computed for at least one collagen fiber family!");
//  }
}

/*----------------------------------------------------------------------*
 |  MassProductionSingleFiber                     (private)        05/11|
 *----------------------------------------------------------------------*
compute new deposition rate for one fiber family
driving force depends on input, may be S or S^k
I derived a strange form for the driving factor just with S and C:

m^k = mbasal * (1 + K * ( sqrt(a_i^T*C*S*C*S*C*a_i/detC)/lambda_k(t) /homeo - 1))

therefore we need C and S as matrices
*/
void MAT::ConstraintMixture::MassProductionSingleFiber
(
  const int gp,
  LINALG::Matrix<NUM_STRESS_3D,1> C,
  LINALG::Matrix<NUM_STRESS_3D,1> S,
  double* massstress,
  double* massprodcomp,
  LINALG::Matrix<3,1> a,
  const int idfiber
)
{
  LINALG::Matrix<3,3> Smatrix(true);
  Smatrix(0,0) = S(0);
  Smatrix(0,1) = S(3);
  Smatrix(0,2) = S(5);
  Smatrix(1,0) = Smatrix(0,1);
  Smatrix(1,1) = S(1);
  Smatrix(1,2) = S(4);
  Smatrix(2,0) = Smatrix(0,2);
  Smatrix(2,1) = Smatrix(1,2);
  Smatrix(2,2) = S(2);
  LINALG::Matrix<3,3> Cmatrix(true);
  Cmatrix(0,0) = C(0);
  Cmatrix(0,1) = 0.5 * C(3);
  Cmatrix(0,2) = 0.5 * C(5);
  Cmatrix(1,0) = Cmatrix(0,1);
  Cmatrix(1,1) = C(1);
  Cmatrix(1,2) = 0.5 * C(4);
  Cmatrix(2,0) = Cmatrix(0,2);
  Cmatrix(2,1) = Cmatrix(1,2);
  Cmatrix(2,2) = C(2);
  double detC = C(0)*C(1)*C(2) + 0.25 * C(3)*C(4)*C(5) - 0.25 * C(1)*C(5)*C(5)
        - 0.25 * C(2)*C(3)*C(3) - 0.25 * C(0)*C(4)*C(4);

  LINALG::Matrix<3,1> CSCa(true);
  LINALG::Matrix<3,1> SCa(true);
  LINALG::Matrix<1,1> temp(true);
  CSCa.Multiply(Cmatrix,a);
  SCa.Multiply(Smatrix,CSCa);
  CSCa.Multiply(Cmatrix,SCa);
  temp.MultiplyTN(1.0/detC,SCa,CSCa);

  LINALG::Matrix<4,1> currentstretch(true);
  history_->back().GetStretches(gp,&currentstretch);

  (*massstress) = temp(0);
  (*massstress) = sqrt(*massstress) / currentstretch(idfiber);
  if (params_->homstress_ != 0.0)
    (*massprodcomp) = massprodbasal_ * (1.0 + params_->growthfactor_ * ((*massstress) / params_->homstress_ - 1.0));
  else
    (*massprodcomp) = massprodbasal_ * (1.0 + params_->growthfactor_ * (*massstress));
}

/*----------------------------------------------------------------------*
 |  EvaluateImplicitAll                           (private)        05/11|
 *----------------------------------------------------------------------*
 evaluate stress and cmat for implicit integration
 driving force of massproduction is the total stress S
 */
void MAT::ConstraintMixture::EvaluateImplicitAll
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  const int gp,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress,
  double dt,
  double time,
  LINALG::Matrix<4,1> massprodcomp,
  LINALG::Matrix<4,1> massstress
)
{
  //--------------------------------------------------------------------------------------
  // build identity tensor I
  LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
  C.Scale(2.0);
  C += Id;

  //--------------------------------------------------------------------------------------
  // store stress and stretch in a matrix to do matrix-matrix multiplications
  LINALG::Matrix<3,3> Smatrix(true);
  Smatrix(0,0) = (*stress)(0);
  Smatrix(0,1) = (*stress)(3);
  Smatrix(0,2) = (*stress)(5);
  Smatrix(1,0) = Smatrix(0,1);
  Smatrix(1,1) = (*stress)(1);
  Smatrix(1,2) = (*stress)(4);
  Smatrix(2,0) = Smatrix(0,2);
  Smatrix(2,1) = Smatrix(1,2);
  Smatrix(2,2) = (*stress)(2);
  LINALG::Matrix<3,3> Cmatrix(true);
  Cmatrix(0,0) = C(0);
  Cmatrix(0,1) = 0.5 * C(3);
  Cmatrix(0,2) = 0.5 * C(5);
  Cmatrix(1,0) = Cmatrix(0,1);
  Cmatrix(1,1) = C(1);
  Cmatrix(1,2) = 0.5 * C(4);
  Cmatrix(2,0) = Cmatrix(0,2);
  Cmatrix(2,1) = Cmatrix(1,2);
  Cmatrix(2,2) = C(2);

  //--------------------------------------------------------------------------------------
  // some variables
  const int firstiter = 0;
  double prestretchcollagen = params_->prestretchcollagen_;
  // store actual collagen stretches, do not change anymore
  LINALG::Matrix<4,1> actcollstretch(true);
  history_->back().GetStretches(gp,&actcollstretch);

  //--------------------------------------------------------------------------------------
  // prestress time
//        const ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
//        INPAR::STR::PreStress pstype = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(pslist,"PRESTRESS");
//        double pstime = -1.0 * params_->lifetime_ - dt;
//        if (pstype == INPAR::STR::prestress_mulf)
//          pstime = pslist.get<double>("PRESTRESSTIME");
  // we already have prestretch in prestress time, thus the prestretch from the inputfile is not applied
//        if (time <= pstime + 1.0e-12)
//          prestretchcollagen = 1.0;

  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatelastic(*cmat);
  // isotropic invariant
  const double I3 = C(0)*C(1)*C(2)
        + 0.25 * C(3)*C(4)*C(5)
        - 0.25 * C(1)*C(5)*C(5)
        - 0.25 * C(2)*C(3)*C(3)
        - 0.25 * C(0)*C(4)*C(4);    // 3rd invariant, determinant
  const double J = sqrt(I3);     // determinant of F

  //-------------------------------------
  // invert C
  LINALG::Matrix<NUM_STRESS_3D,1> Cinv(true);
  Cinv(0) = C(1)*C(2) - 0.25*C(4)*C(4);
  Cinv(1) = C(0)*C(2) - 0.25*C(5)*C(5);
  Cinv(2) = C(0)*C(1) - 0.25*C(3)*C(3);
  Cinv(3) = 0.25*C(5)*C(4) - 0.5*C(3)*C(2);
  Cinv(4) = 0.25*C(3)*C(5) - 0.5*C(0)*C(4);
  Cinv(5) = 0.25*C(3)*C(4) - 0.5*C(5)*C(1);
  Cinv.Scale(1.0/I3);

  // determine residual
  LINALG::Matrix<10,1> Residual(true);
  // first six components are zero as the start value for stress is the computed value
  LINALG::Matrix<4,1> massprod(true);
  history_->back().GetMass(gp,&massprod);
  Residual(6) = massprod(0) - massprodcomp(0);
  Residual(7) = massprod(1) - massprodcomp(1);
  Residual(8) = massprod(2) - massprodcomp(2);
  Residual(9) = massprod(3) - massprodcomp(3);

  //--------------------------------------------------------------------------------------
  // local Newton iteration
  int localistep = 0;
  int maxstep = 50;
  while (Residual.Norm2() > params_->abstol_ && localistep < maxstep)
  {
    localistep += 1;
    //--------------------------------------------------------------------------------------
    // derivative of residual
    LINALG::Matrix<10,10> DResidual(true);
    for (int id = 0; id < 10; id++)
      DResidual(id,id) = 1.0;

    // derivative of stress part with respect to mass production
    // do it for all 4 fiber families
    double stretch = prestretchcollagen/actcollstretch(0);
    LINALG::Matrix<NUM_STRESS_3D,1> temp(true);
    GradStressDMass(glstrain, &temp, Cinv, a1_->at(gp), stretch, J, dt, true);
    for (int id = 0; id < 6; id++)
      DResidual(id,6) = - temp(id);

    stretch = prestretchcollagen/actcollstretch(1);
    GradStressDMass(glstrain, &temp, Cinv, a2_->at(gp), stretch, J, dt, true);
    for (int id = 0; id < 6; id++)
      DResidual(id,7) = - temp(id);

    stretch = prestretchcollagen/actcollstretch(2);
    GradStressDMass(glstrain, &temp, Cinv, a3_->at(gp), stretch, J, dt, true);
    for (int id = 0; id < 6; id++)
      DResidual(id,8) = - temp(id);

    stretch = prestretchcollagen/actcollstretch(3);
    GradStressDMass(glstrain, &temp, Cinv, a4_->at(gp), stretch, J, dt, true);
    for (int id = 0; id < 6; id++)
      DResidual(id,9) = - temp(id);

    // derivative of mass production with respect to stress
    GradMassDStress(&temp, Cmatrix, Smatrix, a1_->at(gp), J, massstress(0), actcollstretch(0));
    for (int id = 0; id < 6; id++)
      DResidual(6,id) = - temp(id);

    GradMassDStress(&temp, Cmatrix, Smatrix, a2_->at(gp), J, massstress(1), actcollstretch(1));
    for (int id = 0; id < 6; id++)
      DResidual(7,id) = - temp(id);

    GradMassDStress(&temp, Cmatrix, Smatrix, a3_->at(gp), J, massstress(2), actcollstretch(2));
    for (int id = 0; id < 6; id++)
      DResidual(8,id) = - temp(id);

    GradMassDStress(&temp, Cmatrix, Smatrix, a4_->at(gp), J, massstress(3), actcollstretch(3));
    for (int id = 0; id < 6; id++)
      DResidual(9,id) = - temp(id);

    //----------------------------------------------------
    // solve linear system of equations: gradF * incr = -F
    //----------------------------------------------------
    // F = F*-1.0
    Residual.Scale(-1.0);
    LINALG::Matrix<10,1> increment(true);
    // solve A.X=B
    LINALG::FixedSizeSerialDenseSolver<10,10,1> solver;
    solver.SetMatrix(DResidual);              // set A=DResidual
    solver.SetVectors(increment, Residual);           // set X=increment, B=Residual
    solver.FactorWithEquilibration(true); // "some easy type of preconditioning" (Michael)
    int err2 = solver.Factor();           // ?
    int err = solver.Solve();             // X = A^-1 B
    if ((err!=0) || (err2!=0))
      dserror("solving linear system in Newton-Raphson method for implicit integration failed");

    // damping strategy
    double omega = 2.0;
    LINALG::Matrix<6,1> stepstress(true);
    LINALG::Matrix<4,1> oldmass(massprod);
    LINALG::Matrix<10,1> Residualtemp(Residual);
    double omegamin = 1.0/64.0;
    while (Residualtemp.Norm2() > (1.0-0.5*omega)*Residual.Norm2() && omega > omegamin)
    {
      // update of theta
      omega = omega/2.0;
      for (int id = 0; id < 6; id++) stepstress(id) = (*stress)(id) + omega*increment(id);
      LINALG::Matrix<4,1> newmass(true);
      newmass(0) = increment(6); newmass(1) = increment(7); newmass(2) = increment(8); newmass(3) = increment(9);
      newmass.Update(1.0,oldmass,omega);
      history_->back().SetMass(gp,newmass);
      massprod = newmass;

      LINALG::Matrix<NUM_STRESS_3D,1> Scomp(true);
      cmatelastic.Scale(0.0);
      EvaluateStress(glstrain, gp, &cmatelastic, &Scomp, firstiter, dt, time);
      MassProduction(gp,C,stepstress,&massstress,&massprodcomp);

      for (int id = 0; id < 6; id++)
        Residualtemp(id) = stepstress(id) - Scomp(id);
      Residualtemp(6) = massprod(0) - massprodcomp(0);
      Residualtemp(7) = massprod(1) - massprodcomp(1);
      Residualtemp(8) = massprod(2) - massprodcomp(2);
      Residualtemp(9) = massprod(3) - massprodcomp(3);
    }
    if (omega <= omegamin && Residualtemp.Norm2() > (1.0-0.5*omega)*Residual.Norm2())
      dserror("no damping coefficient found");

    for (int id = 0; id < 6; id++)
      (*stress)(id) = stepstress(id);
    Residual = Residualtemp;

    Smatrix(0,0) = (*stress)(0);
    Smatrix(0,1) = (*stress)(3);
    Smatrix(0,2) = (*stress)(5);
    Smatrix(1,0) = Smatrix(0,1);
    Smatrix(1,1) = (*stress)(1);
    Smatrix(1,2) = (*stress)(4);
    Smatrix(2,0) = Smatrix(0,2);
    Smatrix(2,1) = Smatrix(1,2);
    Smatrix(2,2) = (*stress)(2);

    if ((massprod(0) < 0.0) || (massprod(1) < 0.0) ||
        (massprod(2) < 0.0) || (massprod(3) < 0.0))
    {
      cout << "1: " << massprod(0) << endl;
      cout << "2: " << massprod(1) << endl;
      cout << "3: " << massprod(2) << endl;
      cout << "4: " << massprod(3) << endl;
      dserror("negative mass production computed for at least one collagen fiber family!");
    }

  } //while loop
  if (localistep == maxstep && Residual.Norm2() > params_->abstol_)
    dserror("local Newton iteration did not converge %e", Residual.Norm2());

  //--------------------------------------------------------------------------------------
  // compute cmat
  // right handside of the linear equations
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> RHS(cmatelastic);
  // left matrix of the linear equations
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> LM(true);
  for (int id = 0; id < 6; id++) LM(id,id) = 1.0;

  // Fiber1
  double stretch = prestretchcollagen/actcollstretch(0);
  LINALG::Matrix<6,1> dmassdstress(true);
  LINALG::Matrix<6,1> dmassdstretch(true);
  LINALG::Matrix<6,1> dstressdmass(true);
  GradMassDStretch(&dmassdstretch, Cmatrix, Smatrix, Cinv, a1_->at(gp), J, massstress(0), actcollstretch(0), dt);
  GradMassDStress(&dmassdstress, Cmatrix, Smatrix, a1_->at(gp), J, massstress(0), actcollstretch(0));
  GradStressDMass(glstrain, &dstressdmass, Cinv, a1_->at(gp), stretch, J, dt, true);
  RHS.MultiplyNT(2.0,dstressdmass,dmassdstretch,1.0);
  LM.MultiplyNT(-1.0,dstressdmass,dmassdstress,1.0);

  // Fiber2
  stretch = prestretchcollagen/actcollstretch(1);
  GradMassDStretch(&dmassdstretch, Cmatrix, Smatrix, Cinv, a2_->at(gp), J, massstress(1), actcollstretch(1), dt);
  GradMassDStress(&dmassdstress, Cmatrix, Smatrix, a2_->at(gp), J, massstress(1), actcollstretch(1));
  GradStressDMass(glstrain, &dstressdmass, Cinv, a2_->at(gp), stretch, J, dt, true);
  RHS.MultiplyNT(2.0,dstressdmass,dmassdstretch,1.0);
  LM.MultiplyNT(-1.0,dstressdmass,dmassdstress,1.0);

  // Fiber3
  stretch = prestretchcollagen/actcollstretch(2);
  GradMassDStretch(&dmassdstretch, Cmatrix, Smatrix, Cinv, a3_->at(gp), J, massstress(2), actcollstretch(2), dt);
  GradMassDStress(&dmassdstress, Cmatrix, Smatrix, a3_->at(gp), J, massstress(2), actcollstretch(2));
  GradStressDMass(glstrain, &dstressdmass, Cinv, a3_->at(gp), stretch, J, dt, true);
  RHS.MultiplyNT(2.0,dstressdmass,dmassdstretch,1.0);
  LM.MultiplyNT(-1.0,dstressdmass,dmassdstress,1.0);

  // Fiber4
  stretch = prestretchcollagen/actcollstretch(3);
  GradMassDStretch(&dmassdstretch, Cmatrix, Smatrix, Cinv, a4_->at(gp), J, massstress(3), actcollstretch(3), dt);
  GradMassDStress(&dmassdstress, Cmatrix, Smatrix, a4_->at(gp), J, massstress(3), actcollstretch(3));
  GradStressDMass(glstrain, &dstressdmass, Cinv, a4_->at(gp), stretch, J, dt, true);
  RHS.MultiplyNT(2.0,dstressdmass,dmassdstretch,1.0);
  LM.MultiplyNT(-1.0,dstressdmass,dmassdstress,1.0);

  (*cmat).Scale(0.0);
  //----------------------------------------------------
  // solve linear system of equations: A.X=B
  //----------------------------------------------------
  LINALG::FixedSizeSerialDenseSolver<6,6,6> solver;
  solver.SetMatrix(LM);              // set A=LM
  solver.SetVectors(*cmat, RHS);           // set X=increment, B=RHS
  solver.FactorWithEquilibration(true); // "some easy type of preconditioning" (Michael)
  int err2 = solver.Factor();           // ?
  int err = solver.Solve();             // X = A^-1 B
  if ((err!=0) || (err2!=0))
    dserror("solving linear system for cmat failed");

  vismassstress_->at(gp)(0) = massstress(0);
  vismassstress_->at(gp)(1) = massstress(1);
  vismassstress_->at(gp)(2) = massstress(2);
}

/*----------------------------------------------------------------------*
 |  EvaluateImplicitSingle                        (private)        05/11|
 *----------------------------------------------------------------------*
evaluate stress and cmat for implicit integration
driving force of massproduction is the fiber stress S^k
*/
void MAT::ConstraintMixture::EvaluateImplicitSingle
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  const int gp,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress,
  double dt,
  double time,
  LINALG::Matrix<4,1> massprodcomp,
  LINALG::Matrix<4,1> massstress
)
{
  //--------------------------------------------------------------------------------------
  // build identity tensor I
  LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
  C.Scale(2.0);
  C += Id;

  //--------------------------------------------------------------------------------------
  // store stretch in a matrix to do matrix-matrix multiplications
  LINALG::Matrix<3,3> Cmatrix(true);
  Cmatrix(0,0) = C(0);
  Cmatrix(0,1) = 0.5 * C(3);
  Cmatrix(0,2) = 0.5 * C(5);
  Cmatrix(1,0) = Cmatrix(0,1);
  Cmatrix(1,1) = C(1);
  Cmatrix(1,2) = 0.5 * C(4);
  Cmatrix(2,0) = Cmatrix(0,2);
  Cmatrix(2,1) = Cmatrix(1,2);
  Cmatrix(2,2) = C(2);

  // isotropic invariant
  const double I3 = C(0)*C(1)*C(2)
        + 0.25 * C(3)*C(4)*C(5)
        - 0.25 * C(1)*C(5)*C(5)
        - 0.25 * C(2)*C(3)*C(3)
        - 0.25 * C(0)*C(4)*C(4);    // 3rd invariant, determinant
  const double J = sqrt(I3);     // determinant of F

  //-------------------------------------
  // invert C
  LINALG::Matrix<NUM_STRESS_3D,1> Cinv(true);
  Cinv(0) = C(1)*C(2) - 0.25*C(4)*C(4);
  Cinv(1) = C(0)*C(2) - 0.25*C(5)*C(5);
  Cinv(2) = C(0)*C(1) - 0.25*C(3)*C(3);
  Cinv(3) = 0.25*C(5)*C(4) - 0.5*C(3)*C(2);
  Cinv(4) = 0.25*C(3)*C(5) - 0.5*C(0)*C(4);
  Cinv(5) = 0.25*C(3)*C(4) - 0.5*C(5)*C(1);
  Cinv.Scale(1.0/I3);

  //--------------------------------------------------------------------------------------
  // some variables
  const int firstiter = 0;
  double prestretchcollagen = params_->prestretchcollagen_;
  double currmassdens = 0.0;
  double qdegrad = 0.0;
  Degradation(0.0,&qdegrad);
  double density = params_->density_;
  // store actual collagen stretches, do not change anymore
  LINALG::Matrix<4,1> actcollstretch(true);
  history_->back().GetStretches(gp,&actcollstretch);
  LINALG::Matrix<4,1> massprod(true);
  history_->back().GetMass(gp,&massprod);

  //--------------------------------------------------------------------------------------
  // prestress time
//        const ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
//        INPAR::STR::PreStress pstype = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(pslist,"PRESTRESS");
//        double pstime = -1.0 * params_->lifetime_ - dt;
//        if (pstype == INPAR::STR::prestress_mulf)
//          pstime = pslist.get<double>("PRESTRESSTIME");
  // we already have prestretch in prestress time, thus the prestretch from the inputfile is not applied
//        if (time <= pstime + 1.0e-12)
//          prestretchcollagen = 1.0;

  // set stress and cmat to zero, as they are not zero here
  (*stress).Scale(0.0);
  (*cmat).Scale(0.0);

  // everything related to the fiber families
  for (int idfiber = 0; idfiber < 4; idfiber++)
  {
    LINALG::Matrix<3,1> a(true);
    if (idfiber == 0)
      a = a1_->at(gp);
    else if (idfiber == 1)
      a = a2_->at(gp);
    else if (idfiber == 2)
      a = a3_->at(gp);
    else
      a = a4_->at(gp);

    LINALG::Matrix<NUM_STRESS_3D,1> stressfiber(true);
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatfiber(true);
    double currmassdensfiber = 0.0;
    EvaluateFiberFamily(glstrain,gp,&cmatfiber,&stressfiber,a,&currmassdensfiber,firstiter,dt,time,idfiber);
    LINALG::Matrix<3,3> Smatrix(true);
    Smatrix(0,0) = stressfiber(0);
    Smatrix(0,1) = stressfiber(3);
    Smatrix(0,2) = stressfiber(5);
    Smatrix(1,0) = Smatrix(0,1);
    Smatrix(1,1) = stressfiber(1);
    Smatrix(1,2) = stressfiber(4);
    Smatrix(2,0) = Smatrix(0,2);
    Smatrix(2,1) = Smatrix(1,2);
    Smatrix(2,2) = stressfiber(2);

    // determine residual
    LINALG::Matrix<7,1> Residual(true);
    // first six components are zero as the start value for stress is the computed value
    Residual(6) = massprod(idfiber) - massprodcomp(idfiber);

    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatelastic(cmatfiber);

    //--------------------------------------------------------------------------------------
    // local Newton iteration
    int localistep = 0;
    int maxstep = 50;
    while (Residual.Norm2() > params_->abstol_ && localistep < maxstep)
    {
      localistep += 1;
      //--------------------------------------------------------------------------------------
      // derivative of residual
      LINALG::Matrix<7,7> DResidual(true);
      for (int id = 0; id < 7; id++)
        DResidual(id,id) = 1.0;

      // derivative of stress part with respect to mass production
      double stretch = prestretchcollagen/actcollstretch(idfiber);
      LINALG::Matrix<NUM_STRESS_3D,1> temp(true);
      GradStressDMass(glstrain, &temp, Cinv, a, stretch, J, dt,false);
      for (int id = 0; id < 6; id++)
        DResidual(id,6) = - temp(id);

      // derivative of mass production with respect to stress
      GradMassDStress(&temp, Cmatrix, Smatrix, a, J, massstress(idfiber), actcollstretch(idfiber));
      for (int id = 0; id < 6; id++)
        DResidual(6,id) = - temp(id);

      //----------------------------------------------------
      // solve linear system of equations: gradF * incr = -F
      //----------------------------------------------------
      // F = F*-1.0
      Residual.Scale(-1.0);
      LINALG::Matrix<7,1> increment(true);
      // solve A.X=B
      LINALG::FixedSizeSerialDenseSolver<7,7,1> solver;
      solver.SetMatrix(DResidual);              // set A=DResidual
      solver.SetVectors(increment, Residual);           // set X=increment, B=Residual
      solver.FactorWithEquilibration(true); // "some easy type of preconditioning" (Michael)
      int err2 = solver.Factor();           // ?
      int err = solver.Solve();             // X = A^-1 B
      if ((err!=0) || (err2!=0))
        dserror("solving linear system in Newton-Raphson method for implicit integration failed");

      // damping strategy
      double omega = 2.0;
      LINALG::Matrix<6,1> stepstress(true);
      double oldmass = massprod(idfiber);
      LINALG::Matrix<7,1> Residualtemp(Residual);
      double omegamin = 1.0/64.0;
      while (Residualtemp.Norm2() > (1.0-0.5*omega)*Residual.Norm2() && omega > omegamin)
      {
        // update of theta
        omega = omega/2.0;
        for (int id = 0; id < 6; id++) stepstress(id) = stressfiber(id) + omega*increment(id);
        double newmass = omega * increment(6) + oldmass;
        massprod(idfiber) = newmass;
        history_->back().SetMass(gp,massprod);

        LINALG::Matrix<NUM_STRESS_3D,1> Scomp(true);
        cmatelastic.Scale(0.0);
        currmassdensfiber = 0.0;
        EvaluateFiberFamily(glstrain, gp, &cmatelastic, &Scomp, a, &currmassdensfiber, firstiter, dt, time, idfiber);
        double massprodcompfiber;
        double massstressfiber;
        MassProductionSingleFiber(gp, C, stepstress, &massstressfiber, &massprodcompfiber, a, idfiber);
        massprodcomp(idfiber) = massprodcompfiber;
        massstress(idfiber) = massstressfiber;

        for (int id = 0; id < 6; id++)
          Residualtemp(id) = stepstress(id) - Scomp(id);
        Residualtemp(6) = massprod(idfiber) - massprodcomp(idfiber);
      }
      if (omega <= omegamin && Residualtemp.Norm2() > (1.0-0.5*omega)*Residual.Norm2())
        dserror("no damping coefficient found");

      stressfiber = stepstress;
      Residual = Residualtemp;

      Smatrix(0,0) = stressfiber(0);
      Smatrix(0,1) = stressfiber(3);
      Smatrix(0,2) = stressfiber(5);
      Smatrix(1,0) = Smatrix(0,1);
      Smatrix(1,1) = stressfiber(1);
      Smatrix(1,2) = stressfiber(4);
      Smatrix(2,0) = Smatrix(0,2);
      Smatrix(2,1) = Smatrix(1,2);
      Smatrix(2,2) = stressfiber(2);

      if ((massprod(idfiber) < 0.0))
      {
        cout << idfiber+1 << ": " << massprod(idfiber) << endl;
        dserror("negative mass production computed for one collagen fiber family!");
      }

    } //while loop
    if (localistep == maxstep && Residual.Norm2() > params_->abstol_)
      dserror("local Newton iteration did not converge %e", Residual.Norm2());

    //--------------------------------------------------------------------------------------
    // compute cmat
    // right handside of the linear equations
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> RHS(cmatelastic);
    // left matrix of the linear equations
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> LM(true);
    for (int id = 0; id < 6; id++) LM(id,id) = 1.0;

    double stretch = prestretchcollagen/actcollstretch(idfiber);
    LINALG::Matrix<6,1> dmassdstress(true);
    LINALG::Matrix<6,1> dmassdstretch(true);
    LINALG::Matrix<6,1> dstressdmass(true);
    GradMassDStretch(&dmassdstretch, Cmatrix, Smatrix, Cinv, a, J, massstress(idfiber), actcollstretch(idfiber), dt);
    GradMassDStress(&dmassdstress, Cmatrix, Smatrix, a, J, massstress(idfiber), actcollstretch(idfiber));
    GradStressDMass(glstrain, &dstressdmass, Cinv, a, stretch, J, dt, false);
    RHS.MultiplyNT(2.0,dstressdmass,dmassdstretch,1.0);
    LM.MultiplyNT(-1.0,dstressdmass,dmassdstress,1.0);

    cmatfiber.Scale(0.0);
    //----------------------------------------------------
    // solve linear system of equations: A.X=B
    //----------------------------------------------------
    LINALG::FixedSizeSerialDenseSolver<6,6,6> solver;
    solver.SetMatrix(LM);              // set A=LM
    solver.SetVectors(cmatfiber, RHS);           // set X=increment, B=RHS
    solver.FactorWithEquilibration(true); // "some easy type of preconditioning" (Michael)
    int err2 = solver.Factor();           // ?
    int err = solver.Solve();             // X = A^-1 B
    if ((err!=0) || (err2!=0))
      dserror("solving linear system for cmat failed");

    (*stress) += stressfiber;
    (*cmat) += cmatfiber;
    // volumetric part, that is related to this fiber family
    (*cmat).MultiplyNT(-2.0*dt/density*qdegrad*params_->kappa_*J,Cinv,dmassdstretch,1.0);
    currmassdens += currmassdensfiber;
  }

  // elastin
  double refmassdenselastin = params_->phielastin_ * density;
  double prestretchelastin = params_->prestretchelastin_;
  // account for isotropic prestretch of elastin
  LINALG::Matrix<NUM_STRESS_3D,1> glstrainiso(*glstrain);
  glstrainiso.Scale(pow(prestretchelastin,2));
  LINALG::Matrix<NUM_STRESS_3D,1> Siso(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatiso(true);
  EvaluateElastin(& glstrainiso, & cmatiso, & Siso);
  Siso.Scale(refmassdenselastin/density*pow(prestretchelastin,2));
  (*stress) += Siso;
  cmatiso.Scale(refmassdenselastin/density*pow(prestretchelastin,4));
  (*cmat) += cmatiso;
  currmassdens += refmassdenselastin;

  // volumetric part
  LINALG::Matrix<NUM_STRESS_3D,1> Svol(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatvol(true);
  EvaluateVolumetric(glstrain, &cmatvol, &Svol, currmassdens, density);
  (*stress) += Svol;
  (*cmat) += cmatvol;

  vismassstress_->at(gp)(0) = massstress(0);
  vismassstress_->at(gp)(1) = massstress(1);
  vismassstress_->at(gp)(2) = massstress(2);
}

/*----------------------------------------------------------------------*
 |  GradStressDMass                               (private)        05/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::GradStressDMass
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  LINALG::Matrix<NUM_STRESS_3D,1>* derivative,
  LINALG::Matrix<NUM_STRESS_3D,1> Cinv,
  LINALG::Matrix<3,1> a,
  double stretch,
  double J,
  double dt,
  bool option
)
{
  double density = params_->density_;
  double qdegrad;
  Degradation(0.0,&qdegrad);
  LINALG::Matrix<NUM_STRESS_3D,1> Saniso(true);
  EvaluateSingleFiber(glstrain, NULL, & Saniso, a, stretch);
  double facS = pow(stretch,2) * qdegrad / density * dt;
  Saniso.Scale(facS);
  if (option)
    (*derivative).Update(-dt/density*qdegrad*params_->kappa_*J,Cinv,1.0,Saniso);
  else
    *derivative = Saniso;
}

/*----------------------------------------------------------------------*
 |  GradMassDStress                               (private)        05/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::GradMassDStress
(
  LINALG::Matrix<NUM_STRESS_3D,1>* derivative,
  LINALG::Matrix<3,3> Cmatrix,
  LINALG::Matrix<3,3> Smatrix,
  LINALG::Matrix<3,1> a,
  double J,
  double massstress,
  double actcollstretch
)
{
  // include the case homstress = 0.0!!
  double homstressfixed = 1.0;
  if (params_->homstress_ != 0.0) homstressfixed = params_->homstress_;

  LINALG::Matrix<3,1> Ca(true);
  LINALG::Matrix<3,1> CSCa(true);
  Ca.Multiply(Cmatrix,a);
  LINALG::Matrix<3,1> temp1(true);
  temp1.Multiply(Smatrix,Ca);
  CSCa.Multiply(Cmatrix,temp1);
  double fac = massprodbasal_ * params_->growthfactor_ / homstressfixed /pow(J,2)
             / massstress / pow(actcollstretch,2);
  (*derivative)(0) = fac * Ca(0) * CSCa(0);
  (*derivative)(1) = fac * Ca(1) * CSCa(1);
  (*derivative)(2) = fac * Ca(2) * CSCa(2);
  (*derivative)(3) = fac * 0.5 * (Ca(0) * CSCa(1) + Ca(1) * CSCa(0));
  (*derivative)(4) = fac * 0.5 * (Ca(1) * CSCa(2) + Ca(2) * CSCa(1));
  (*derivative)(5) = fac * 0.5 * (Ca(0) * CSCa(2) + Ca(2) * CSCa(0));
}

/*----------------------------------------------------------------------*
 |  GradMassDStretch                              (private)        05/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::GradMassDStretch
(
  LINALG::Matrix<NUM_STRESS_3D,1>* derivative,
  LINALG::Matrix<3,3> Cmatrix,
  LINALG::Matrix<3,3> Smatrix,
  LINALG::Matrix<NUM_STRESS_3D,1> Cinv,
  LINALG::Matrix<3,1> a,
  double J,
  double massstress,
  double actcollstretch,
  double dt
)
{
  // include the case homstress = 0.0!!
  double homstressfixed = 1.0;
  if (params_->homstress_ != 0.0) homstressfixed = params_->homstress_;

  LINALG::Matrix<3,1> SCa(true);
  LINALG::Matrix<3,1> SCSCa(true);
  LINALG::Matrix<3,1> temp(true);
  temp.Multiply(Cmatrix,a);
  SCa.Multiply(Smatrix,temp);
  temp.Multiply(Cmatrix,SCa);
  SCSCa.Multiply(Smatrix,temp);

  (*derivative)(0) = a(0)*a(0); (*derivative)(1) = a(1)*a(1);
  (*derivative)(2) = a(2)*a(2); (*derivative)(3) = a(0)*a(1);
  (*derivative)(4) = a(1)*a(2); (*derivative)(5) = a(0)*a(2);
  (*derivative).Update(1.0,Cinv,1.0/pow(actcollstretch,2));
  (*derivative).Scale(- 1.0 * massstress);
  double fac = 1.0 / massstress / pow(actcollstretch,2) / pow(J,2);
  (*derivative)(0) += fac * (2.0 * a(0) * SCSCa(0) + SCa(0) * SCa(0));
  (*derivative)(1) += fac * (2.0 * a(1) * SCSCa(1) + SCa(1) * SCa(1));
  (*derivative)(2) += fac * (2.0 * a(2) * SCSCa(2) + SCa(2) * SCa(2));
  (*derivative)(3) += fac * (a(0) * SCSCa(1) + a(1) * SCSCa(0) + SCa(0) * SCa(1));
  (*derivative)(4) += fac * (a(1) * SCSCa(2) + a(2) * SCSCa(1) + SCa(1) * SCa(2));
  (*derivative)(5) += fac * (a(0) * SCSCa(2) + a(2) * SCSCa(0) + SCa(0) * SCa(2));
  (*derivative).Scale(0.5*massprodbasal_*params_->growthfactor_/homstressfixed);
}

/*----------------------------------------------------------------------*
 |  EvaluateFiberVecs                             (public)         02/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::EvaluateFiberVecs
(const int gp, const LINALG::Matrix<3,3>& locsys, const LINALG::Matrix<3,3>& defgrd)
{
  // locsys holds the principal directions
  // The deformation gradient (defgrd) is needed in remodeling as then locsys is given in the
  // spatial configuration and thus the fiber vectors have to be pulled back in the reference
  // configuration as the material is evaluated there.
  // If this function is called during Setup defgrd should be replaced by the Identity.

  const double gamma = (45*PI)/180.; //angle for diagonal fibers
  LINALG::Matrix<3,1> ca1(true);
  LINALG::Matrix<3,1> ca2(true);
  LINALG::Matrix<3,1> ca3(true);
  LINALG::Matrix<3,1> ca4(true);

  for (int i = 0; i < 3; i++) {
    // store old fiber directions
//    olda1_->at(gp)[i] = a1_->at(gp)[i];
//    olda2_->at(gp)[i] = a2_->at(gp)[i];
//    olda3_->at(gp)[i] = a3_->at(gp)[i];
//    olda4_->at(gp)[i] = a4_->at(gp)[i];
    // a1 = e3, circumferential direction, used for collagen and smooth muscle
    ca1(i) = locsys(i,2);
    // a2 = e2
    ca2(i) = locsys(i,1);
    // a3 = cos gamma e3 + sin gamma e2
    ca3(i) = cos(gamma)*locsys(i,2) + sin(gamma)*locsys(i,1);
    // a4 = cos gamma e3 - sin gamma e2
    ca4(i) = cos(gamma)*locsys(i,2) - sin(gamma)*locsys(i,1);
  }

  // pull back in reference configuration
  LINALG::Matrix<3,1> a1_0(true);
  LINALG::Matrix<3,1> a2_0(true);
  LINALG::Matrix<3,1> a3_0(true);
  LINALG::Matrix<3,1> a4_0(true);
  LINALG::Matrix<3,3> idefgrd(false);
  idefgrd.Invert(defgrd);
  a1_0.Multiply(idefgrd,ca1);
  a2_0.Multiply(idefgrd,ca2);
  a3_0.Multiply(idefgrd,ca3);
  a4_0.Multiply(idefgrd,ca4);

  // normalize vectors
  double a1_0norm = a1_0.Norm2();
  a1_->at(gp).Update(1.0/a1_0norm,a1_0);
  double a2_0norm = a2_0.Norm2();
  a2_->at(gp).Update(1.0/a2_0norm,a2_0);
  double a3_0norm = a3_0.Norm2();
  a3_->at(gp).Update(1.0/a3_0norm,a3_0);
  double a4_0norm = a4_0.Norm2();
  a4_->at(gp).Update(1.0/a4_0norm,a4_0);

  return;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
//History
//----------------------------------------------------------------------
//----------------------------------------------------------------------

MAT::ConstraintMixtureHistoryType MAT::ConstraintMixtureHistoryType::instance_;

DRT::ParObject* MAT::ConstraintMixtureHistoryType::Create( const std::vector<char> & data )
{
  MAT::ConstraintMixtureHistory* cmhis = new MAT::ConstraintMixtureHistory();
  cmhis->Unpack(data);
  return cmhis;
}

/*----------------------------------------------------------------------*
 |  History: Pack                                 (public)         03/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixtureHistory::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // Pack internal variables
  AddtoPack(data,depositiontime_);
  AddtoPack(data,dt_);
  AddtoPack(data,numgp_);
  for (int gp = 0; gp < numgp_; ++gp)
  {
    AddtoPack(data,collagenstretch1_->at(gp));
    AddtoPack(data,collagenstretch2_->at(gp));
    AddtoPack(data,collagenstretch3_->at(gp));
    AddtoPack(data,collagenstretch4_->at(gp));
    AddtoPack(data,massprod1_->at(gp));
    AddtoPack(data,massprod2_->at(gp));
    AddtoPack(data,massprod3_->at(gp));
    AddtoPack(data,massprod4_->at(gp));
  }

  return;
}

/*----------------------------------------------------------------------*
 |  History: Unpack                               (public)         03/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixtureHistory::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // unpack internal variables
  double a;
  ExtractfromPack(position,data,a);
  depositiontime_ = a;
  ExtractfromPack(position,data,a);
  dt_ = a;
  int b;
  ExtractfromPack(position,data,b);
  numgp_ = b;

  collagenstretch1_ = Teuchos::rcp(new vector<double> (numgp_));
  collagenstretch2_ = Teuchos::rcp(new vector<double> (numgp_));
  collagenstretch3_ = Teuchos::rcp(new vector<double> (numgp_));
  collagenstretch4_ = Teuchos::rcp(new vector<double> (numgp_));
  massprod1_ = Teuchos::rcp(new vector<double> (numgp_));
  massprod2_ = Teuchos::rcp(new vector<double> (numgp_));
  massprod3_ = Teuchos::rcp(new vector<double> (numgp_));
  massprod4_ = Teuchos::rcp(new vector<double> (numgp_));
  for (int gp = 0; gp < numgp_; ++gp) {
    ExtractfromPack(position,data,a);
    collagenstretch1_->at(gp) = a;
    ExtractfromPack(position,data,a);
    collagenstretch2_->at(gp) = a;
    ExtractfromPack(position,data,a);
    collagenstretch3_->at(gp) = a;
    ExtractfromPack(position,data,a);
    collagenstretch4_->at(gp) = a;
    ExtractfromPack(position,data,a);
    massprod1_->at(gp) = a;
    ExtractfromPack(position,data,a);
    massprod2_->at(gp) = a;
    ExtractfromPack(position,data,a);
    massprod3_->at(gp) = a;
    ExtractfromPack(position,data,a);
    massprod4_->at(gp) = a;
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 |  History: Setup                                (public)         03/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixtureHistory::Setup(const int ngp,const double massprodbasal)
{
  dt_ = 0.0;
  depositiontime_ = 0.0;

  numgp_=ngp;
  // history variables
  collagenstretch1_ = Teuchos::rcp(new vector<double> (numgp_));
  collagenstretch2_ = Teuchos::rcp(new vector<double> (numgp_));
  collagenstretch3_ = Teuchos::rcp(new vector<double> (numgp_));
  collagenstretch4_ = Teuchos::rcp(new vector<double> (numgp_));
  massprod1_ = Teuchos::rcp(new vector<double> (numgp_));
  massprod2_ = Teuchos::rcp(new vector<double> (numgp_));
  massprod3_ = Teuchos::rcp(new vector<double> (numgp_));
  massprod4_ = Teuchos::rcp(new vector<double> (numgp_));
  for (int gp = 0; gp < numgp_; gp++)
  {
    collagenstretch1_->at(gp) = 1.0;
    collagenstretch2_->at(gp) = 1.0;
    collagenstretch3_->at(gp) = 1.0;
    collagenstretch4_->at(gp) = 1.0;
    massprod1_->at(gp) = massprodbasal;
    massprod2_->at(gp) = massprodbasal;
    massprod3_->at(gp) = massprodbasal;
    massprod4_->at(gp) = massprodbasal;
  }
}

/*----------------------------------------------------------------------*
 |  History: SetStretches                         (private)        03/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixtureHistory::SetStretches(int gp, LINALG::Matrix<4,1> stretches)
{
  if (gp < numgp_) {
    collagenstretch1_->at(gp) = stretches(0);
    collagenstretch2_->at(gp) = stretches(1);
    collagenstretch3_->at(gp) = stretches(2);
    collagenstretch4_->at(gp) = stretches(3);
  } else dserror("gp out of range in SetStretches");
}
/*----------------------------------------------------------------------*
 |  History: GetStretches                         (private)        03/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixtureHistory::GetStretches(int gp, LINALG::Matrix<4,1>* stretches)
{
  if (gp < numgp_) {
    (*stretches)(0) = collagenstretch1_->at(gp);
    (*stretches)(1) = collagenstretch2_->at(gp);
    (*stretches)(2) = collagenstretch3_->at(gp);
    (*stretches)(3) = collagenstretch4_->at(gp);
  } else dserror("gp out of range in GetStretches");
}

/*----------------------------------------------------------------------*
 |  History: SetMass                              (private)        03/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixtureHistory::SetMass(int gp, LINALG::Matrix<4,1> massprod)
{
  if (gp < numgp_) {
    massprod1_->at(gp) = massprod(0);
    massprod2_->at(gp) = massprod(1);
    massprod3_->at(gp) = massprod(2);
    massprod4_->at(gp) = massprod(3);
  } else dserror("gp out of range in SetMass");
}
/*----------------------------------------------------------------------*
 |  History: GetMass                              (private)        03/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixtureHistory::GetMass(int gp, LINALG::Matrix<4,1>* massprod)
{
  if (gp < numgp_) {
    (*massprod)(0) = massprod1_->at(gp);
    (*massprod)(1) = massprod2_->at(gp);
    (*massprod)(2) = massprod3_->at(gp);
    (*massprod)(3) = massprod4_->at(gp);
  } else dserror("gp out of range in GetMass");
}

/*----------------------------------------------------------------------*
 |  Debug output to gmsh-file                                      01/11|
 *----------------------------------------------------------------------*
 this needs to be copied to STR::TimInt::OutputStep() to enable debug output
 {
   discret_->SetState("displacement",Dis());
   MAT::ConstraintMixtureOutputToGmsh(discret_, GetStep(), 1);
 }
 just works with strtimint!
 don't forget to include constraintmixture.H */
void MAT::ConstraintMixtureOutputToGmsh
(
  const Teuchos::RCP<DRT::Discretization> dis,
  const int timestep,
  const int iter
)
{
  const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
  // file for stress
  std::stringstream filename;
  filename << filebase << "_visual" << std::setw(3) << setfill('0') << timestep << std::setw(2) << setfill('0') << iter << ".pos";
  std::ofstream f_system(filename.str().c_str());
  stringstream gmshfilecontent;
  gmshfilecontent << "View \" Time: " << timestep << " Iter: " << iter << " \" {" << endl;


  for (int iele=0; iele<dis->NumMyColElements(); ++iele)
  {
    const DRT::Element* actele = dis->lColElement(iele);

    // build current configuration
    vector<int> lm;
    vector<int> lmowner;
    vector<int> lmstride;
    actele->LocationVector(*dis,lm,lmowner,lmstride);
    Teuchos::RCP<const Epetra_Vector> disp = dis->GetState("displacement");
    vector<double> mydisp(lm.size(),0);

    Teuchos::RefCountPtr<MAT::Material> mat = actele->Material();
    MAT::ConstraintMixture* grow = static_cast <MAT::ConstraintMixture*>(mat.get());
    // Teuchos::RCP<vector<double> > mandel = grow->GetVis();

    // material plot at gauss points
    int ngp = grow->Geta1()->size();
    for (int gp = 0; gp < ngp; ++gp){
      vector<double> point = MAT::MatPointCoords(actele,mydisp,gp); //defined in contchainnetw

      // write mandel stress
      LINALG::Matrix<3,1> mandelgp = grow->GetVis(gp);
      gmshfilecontent << "SP(" << scientific << point[0] << ",";
      gmshfilecontent << scientific << point[1] << ",";
      gmshfilecontent << scientific << point[2] << ")";
      gmshfilecontent << "{" << scientific
      << mandelgp(0)
      << "};" << endl;

    }
  }

  gmshfilecontent << "};" << endl;
  f_system << gmshfilecontent.str();
  f_system.close();

  return;
}

#endif // CCADISCRET
