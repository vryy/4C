/*!----------------------------------------------------------------------
\file constraintmixture.cpp
\brief
This file contains routines for constraint mixture growth and remodeling.
example input line
MAT 1 MAT_ConstraintMixture DENS 0.001 MUE 1.0 PHIE 0.08 PREELA 1.0
K1 1.0 K2 1.0 PRECOLL 1.06 K1M 1.0 K2M 1.0 PHIM 1.0 PREMUS 1.0
SMAX 0.0 KAPPA 1.0E6 LIFETIME 5.0 HOMSTR 6.75E4 GROWTHFAC 0.5
STARTTIME 5.0 INTEGRATION Explicit TOL 1.0E-4 GROWTHFORCE Single
INITSTRETCH None DEGOPTION Cos

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


#include <vector>
#include "constraintmixture.H"
#include "../drt_lib/drt_globalproblem.H"
#include "matpar_bundle.H"
#include "../drt_mat/material_service.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_io/io_gmsh.H" // for debug plotting with gmsh
#include "../drt_io/io_control.H" // for debug plotting with gmsh
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H" // for debug plotting with gmsh
#include "../drt_fem_general/drt_utils_integration.H" // for debug plotting with gmsh
#include "../drt_lib/drt_utils.H" // for debug plotting with gmsh
#include "../drt_inpar/inpar_structure.H"  // for pstime

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
  k1muscle_(matdata->GetDouble("K1M")),
  k2muscle_(matdata->GetDouble("K2M")),
  phimuscle_(matdata->GetDouble("PHIM")),
  prestretchmuscle_(matdata->GetDouble("PREMUS")),
  Smax_(matdata->GetDouble("SMAX")),
  kappa_(matdata->GetDouble("KAPPA")),
  lifetime_(matdata->GetDouble("LIFETIME")),
  homstress_(matdata->GetDouble("HOMSTR")),
  growthfactor_(matdata->GetDouble("GROWTHFAC")),
  starttime_(matdata->GetDouble("STARTTIME")),
  integration_(matdata->Get<string>("INTEGRATION")),
  abstol_(matdata->GetDouble("TOL")),
  growthforce_(matdata->Get<string>("GROWTHFORCE")),
  initstretch_(matdata->Get<string>("INITSTRETCH")),
  degoption_(*(matdata->Get<string>("DEGOPTION"))),
  degtol_(1.0e-4)
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
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  int numgp;
  if (!isinit_)
  {
    numgp = 0; // not initialized -> nothing to pack
  } else {
    numgp = a1_->size();   // size is number of gausspoints
  }
  AddtoPack(data, numgp);
  // Pack internal variables
  for (int gp = 0; gp < numgp; ++gp)
  {
    AddtoPack(data, a1_->at(gp));
    AddtoPack(data, a2_->at(gp));
    AddtoPack(data, a3_->at(gp));
    AddtoPack(data, a4_->at(gp));
    AddtoPack(data, vismassstress_->at(gp));
    AddtoPack(data, refmassdens_->at(gp));
    AddtoPack(data, visrefmassdens_->at(gp));
  }
  if (numgp > 0)
  {
    AddtoPack(data,massprodbasal_);

    // Pack history
    int sizehistory = history_->size();
    AddtoPack(data, sizehistory);
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
  ExtractfromPack(position, data, numgp);
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
  visrefmassdens_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> > (numgp));

  for (int gp = 0; gp < numgp; ++gp) {
    LINALG::Matrix<3,1> alin;
    ExtractfromPack(position, data, alin);
    a1_->at(gp) = alin;
    ExtractfromPack(position, data, alin);
    a2_->at(gp) = alin;
    ExtractfromPack(position, data, alin);
    a3_->at(gp) = alin;
    ExtractfromPack(position, data, alin);
    a4_->at(gp) = alin;
    ExtractfromPack(position, data, alin);
    vismassstress_->at(gp) = alin;
    double a;
    ExtractfromPack(position,data,a);
    refmassdens_->at(gp) = a;
    ExtractfromPack(position, data, alin);
    visrefmassdens_->at(gp) = alin;
  }
  double basal;
  ExtractfromPack(position,data,basal);
  massprodbasal_ = basal;

  // unpack history
  int sizehistory;
  ExtractfromPack(position, data, sizehistory);
  history_ = Teuchos::rcp(new vector<ConstraintMixtureHistory> (sizehistory));
  vector<char> datahistory;
  for (int idpast = 0; idpast < sizehistory; idpast++)
  {
    ExtractfromPack(position, data, datahistory);
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
  visrefmassdens_ = Teuchos::rcp(new vector<LINALG::Matrix<3,1> > (numgp));
  for (int gp = 0; gp < numgp; gp++)
  {
    vismassstress_->at(gp)(0) = 0.0;
    vismassstress_->at(gp)(1) = 0.0;
    vismassstress_->at(gp)(2) = 0.0;
    refmassdens_->at(gp) = params_->density_;
    visrefmassdens_->at(gp)(0) = params_->density_*(1.0 - params_->phielastin_ - params_->phimuscle_)/4.0;
    visrefmassdens_->at(gp)(1) = params_->density_*(1.0 - params_->phielastin_ - params_->phimuscle_)/4.0;
    visrefmassdens_->at(gp)(2) = params_->density_*(1.0 - params_->phielastin_ - params_->phimuscle_)/4.0;
//    visrefmassdens_->at(gp)(0) = params_->density_*(1.0 - params_->phielastin_ - params_->phimuscle_)/10.0;
//    visrefmassdens_->at(gp)(1) = params_->density_*(1.0 - params_->phielastin_ - params_->phimuscle_)/10.0;
//    visrefmassdens_->at(gp)(2) = params_->density_*(1.0 - params_->phielastin_ - params_->phimuscle_)/5.0*2.0;
  }

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
  int firstiter = 0;
  if (*params_->integration_ == "Explicit")
    firstiter = 1;
  int numpast = 0;
  if (abs(round(params_->lifetime_ / dt) - params_->lifetime_ / dt) < 1.0e-8){
    numpast = static_cast<int>(round(params_->lifetime_ / dt)) + firstiter;
  } else {
    numpast = static_cast<int>(ceil(params_->lifetime_ / dt)) + firstiter;
  }
  if (params_->degoption_ == "Exp")
  {
    double taumax = - log(params_->degtol_) * params_->lifetime_;
    numpast = static_cast<int>(round(taumax / dt)) + firstiter;
  }

  // prestress time
  INPAR::STR::PreStress pstype = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(timeintegr,"PRESTRESS");
  if (pstype == INPAR::STR::prestress_mulf)
  {
    double pstime = timeintegr.get<double>("PRESTRESSTIME");
    if (pstime > params_->starttime_ + dt)
      dserror("MULF is only working for PRESTRESSTIME smaller than STARTTIME!");
  }

  // basal mass production rate determined by DENS, PHIE and degradation function
  double intdegr = 0.0;
  for (int idpast = 0; idpast < numpast - firstiter; idpast++)
  {
    double degr = 0.0;
    Degradation((numpast-1-idpast)*dt, &degr);
    intdegr += degr * dt;
  }
  massprodbasal_ = (1.0 - params_->phimuscle_ - params_->phielastin_) * params_->density_ / 4.0 / intdegr;
//  massprodbasal_ = (1.0 - params_->phimuscle_ - params_->phielastin_) * params_->density_ / 10.0 / intdegr;

  // history
  history_ = Teuchos::rcp(new vector<ConstraintMixtureHistory> (numpast));
  for (int idpast = 0; idpast < numpast; idpast++)
  {
    history_->at(idpast).Setup(numgp, massprodbasal_);
    history_->at(idpast).SetTime(dt-(numpast-1-idpast)*dt, dt);
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
  history_->back().GetTime(&deptime, &depdt);

  // just do update in case of growth
  if (deptime > params_->starttime_ + 1.0e-12)
  {
    // delete just the steps that surely won't be needed, especially with a smaller timestep later
    // thus reference time is deposition time of last collagen fibers
    double acttime = deptime;
    int eraseiter = 0;
    history_->at(eraseiter).GetTime(&deptime, &depdt);
    double degrad = 0.0;
    Degradation(acttime-deptime, &degrad);
    while (degrad < params_->degtol_ && eraseiter < sizehistory)
    {
      eraseiter +=1;
      history_->at(eraseiter).GetTime(&deptime, &depdt);
      Degradation(acttime-deptime, &degrad);
    }
    if (eraseiter > 0)
    {
      history_->erase(history_->begin(),history_->begin()+eraseiter);
      //cout << "erased " << eraseiter << " history variables" << endl;
    }

    // append new collagen
    ConstraintMixtureHistory newhis;
    newhis.Setup(numgp, massprodbasal_);
    // it is very important to set time and dt to 0.0
    // this makes it clear that this step was created in update and has no reliable content
    // they are not known here either
    // in EvaluateStress this is important, if called after Update
    newhis.SetTime(0.0, 0.0);
    history_->push_back(newhis);

  } else {
    // just adopt deposition time, the rest stays the same
    double newtime = 0.0;
    double newdt = 0.0;
    history_->at(0).GetTime(&newtime, &newdt);
    double degrad = 0.0;
    Degradation(deptime-newtime, &degrad);
    if (degrad < params_->degtol_)
    {
      for (int iter = 0; iter < sizehistory-1; iter++)
      {
        history_->at(iter+1).GetTime(&newtime, &newdt);
        history_->at(iter).SetTime(newtime, newdt);
      }
      history_->back().SetTime(0.0, 0.0);
    } else if (deptime == depdt || (deptime == 2*depdt && *params_->integration_ == "Implicit"))
    {
      // special case of first time step
      ConstraintMixtureHistory newhis;
      newhis.Setup(numgp, massprodbasal_);
      newhis.SetTime(0.0, 0.0);
      history_->push_back(newhis);
    } else {
      dserror("You should not change your timestep size in the case time < starttime! %f", deptime);
    }
  }
}

/*----------------------------------------------------------------------*
 |  Reset internal variables                      (public)         01/12|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::Reset()
{
  history_->back().SetTime(0.0, 0.0);
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
    double temptime = 0.0;
    double tempdt = 0.0;
    history_->back().GetTime(&temptime, &tempdt);
    if (temptime == 0.0 && tempdt == 0.0)
    {
      int numpast = history_->size();
      history_->at(numpast-2).GetTime(&temptime, &tempdt);
      // for restart the function ApplyForceInternal calls the material with the old time
      // (i.e. time = temptime) thus make sure not to store it
      if (time > temptime + 1.0e-12)
      {
        history_->back().SetTime(time, dt);
        // if you change your time step size the basal mass production rate changes
        // basal mass production rate determined by DENS, PHIE and degradation function
        double intdegr = 0.0;
        double degrtime = 0.0;
        double degrdt = 0.0;
        for (int idpast = 0; idpast < numpast - firstiter; idpast++)
        {
          double degr = 0.0;
          history_->at(idpast).GetTime(&degrtime, &degrdt);
          if (firstiter == 1)
          {
            double timeloc = 0.0;
            double dtloc = 0.0;
            history_->at(idpast+1).GetTime(&timeloc, &dtloc);
            degrdt = dtloc;
          }
          Degradation(time-degrtime, &degr);
          intdegr += degr * degrdt;
        }
        massprodbasal_ = (1.0 - params_->phimuscle_ - params_->phielastin_) * params_->density_ / 4.0 / intdegr;
//        massprodbasal_ = (1.0 - params_->phimuscle_ - params_->phielastin_) * params_->density_ / 10.0 / intdegr;
      }
    }
    else if (time > temptime)
    {
      // might be the case in prestressing as Update is not called during prestress
      // thus time has to be adapted and nothing else, as there is no growth & remodeling during prestress
      const ParameterList& pslist = DRT::Problem::Instance()->StructuralDynamicParams();
      INPAR::STR::PreStress pstype = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(pslist,"PRESTRESS");
      if (pstype == INPAR::STR::prestress_mulf)
      {
        // overwrite time
        int sizehistory = history_->size();
        for (int idpast = 0; idpast < sizehistory; idpast++)
        {
          history_->at(idpast).SetTime(time - (sizehistory-1-idpast)*dt, dt);
        }
      }
      else
      {
        // in remodeling time might be wrong depending on the time integration used
        // correct this for the computation but do not store it
        time = temptime;
      }
    }

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
      history_->back().SetStretches(gp, actstretch);
    }
    else if (*params_->initstretch_ == "Homeo")
    {  // this is not working for all material parameters
      int numsteps = history_->size();
      for (int i = 0; i < numsteps; i++)
        history_->at(i).SetStretches(gp, actstretch);
    }

    // start in every iteration from the original value, this is important for implicit only
    LINALG::Matrix<4,1> massprodstart(true);
    for (int id = 0; id < 4; id++)
      massprodstart(id) = massprodbasal_;
//    massprodstart(2) = massprodstart(2)*4;
//    massprodstart(3) = massprodstart(3)*4;
    history_->back().SetMass(gp, massprodstart);

    EvaluateStress(glstrain, gp, cmat, stress, firstiter, time);

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
      EvaluateFiberFamily(glstrain, gp, &cmattemp, &stresstemp, a1_->at(gp), &masstemp, firstiter, time, 0);
      MassProductionSingleFiber(gp, C, stresstemp, &massstresstemp, &massprodtemp, a1_->at(gp), 0);
      massstress(0) = massstresstemp;
      massprodcomp(0) = massprodtemp;
      stresstemp.Scale(0.0);
      EvaluateFiberFamily(glstrain, gp, &cmattemp, &stresstemp, a2_->at(gp), &masstemp, firstiter, time, 1);
      MassProductionSingleFiber(gp, C, stresstemp, &massstresstemp, &massprodtemp, a2_->at(gp), 1);
      massstress(1) = massstresstemp;
      massprodcomp(1) = massprodtemp;
      stresstemp.Scale(0.0);
      EvaluateFiberFamily(glstrain, gp, &cmattemp, &stresstemp, a3_->at(gp), &masstemp, firstiter, time, 2);
      MassProductionSingleFiber(gp, C, stresstemp, &massstresstemp, &massprodtemp, a3_->at(gp), 2);
      massstress(2) = massstresstemp;
      massprodcomp(2) = massprodtemp;
      stresstemp.Scale(0.0);
      EvaluateFiberFamily(glstrain, gp, &cmattemp, &stresstemp, a4_->at(gp), &masstemp, firstiter, time, 3);
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
      }
      else if (*params_->growthforce_ == "All")
      {
        EvaluateImplicitAll(glstrain, gp, cmat, stress, dt, time, massprodcomp, massstress);
      }
      else
      {
        EvaluateImplicitSingle(glstrain, gp, cmat, stress, dt, time);
      }

    } else {
      // visualization of massstresss for the other cases
      if (time > params_->starttime_ + eps || *params_->initstretch_ == "Homeo")
      {
        vismassstress_->at(gp)(0) = massstress(0);
        vismassstress_->at(gp)(1) = massstress(1);
        vismassstress_->at(gp)(2) = massstress(2);
      } else {
        // current stretch in Massproduction is 1.0, thus we have to consider it here
        vismassstress_->at(gp)(0) = massstress(0)/actstretch(0);
        vismassstress_->at(gp)(1) = massstress(1)/actstretch(1);
        vismassstress_->at(gp)(2) = massstress(2)/actstretch(2);
      }
    }
  } else {
    // in case of output everything is fully converged, we just have to evaluate stress etc.
    // should be independent of order of update and output, as new steps are set with dt = 0.0
    // and oldest fibers are carefully erased
    double temptime = 0.0;
    double tempdt = 0.0;
    history_->back().GetTime(&temptime,&tempdt);
    if (time != temptime)
    {
      if (temptime == 0.0)
      {
        int size = history_->size();
        history_->at(size-2).GetTime(&temptime,&tempdt);
        EvaluateStress(glstrain, gp, cmat, stress, firstiter, temptime);
        dserror("has to be checked, update called before output");
      } else
        dserror("times do not match: %f actual time, %f deposition time of last fiber",time,temptime);
    }
    else
      EvaluateStress(glstrain, gp, cmat, stress, firstiter, time);
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
  glstrainiso.Scale(prestretchelastin*prestretchelastin);
  LINALG::Matrix<NUM_STRESS_3D,1> Siso(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatiso(true);
  EvaluateElastin(&glstrainiso, &cmatiso, &Siso);
  Siso.Scale(refmassdenselastin/density*prestretchelastin*prestretchelastin);
  (*stress) = Siso;
  cmatiso.Scale(refmassdenselastin/density*prestretchelastin*prestretchelastin*prestretchelastin*prestretchelastin);
  (*cmat) = cmatiso;
  currmassdens += refmassdenselastin;

  // 2nd step: collagen
  //==========================
  EvaluateFiberFamily(glstrain, gp, cmat, stress, a1_->at(gp), &currmassdens, firstiter, time, 0);

  EvaluateFiberFamily(glstrain, gp, cmat, stress, a2_->at(gp), &currmassdens, firstiter, time, 1);

  EvaluateFiberFamily(glstrain, gp, cmat, stress, a3_->at(gp), &currmassdens, firstiter, time, 2);

  EvaluateFiberFamily(glstrain, gp, cmat, stress, a4_->at(gp), &currmassdens, firstiter, time, 3);

  // 3rd step: smooth muscle
  //==========================
  LINALG::Matrix<NUM_STRESS_3D,1> Smus(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatmus(true);
  EvaluateMuscle(glstrain, &cmatmus, &Smus, gp, &currmassdens);
  (*stress) += Smus;
  (*cmat) += cmatmus;

  // 4th step: volumetric part
  //==========================
  LINALG::Matrix<NUM_STRESS_3D,1> Svol(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatvol(true);
  EvaluateVolumetric(glstrain, &cmatvol, &Svol, currmassdens, density);
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
  double time,
  const int idfiber
 )
{
  //--------------------------------------------------------------------------------------
  // some variables
  double prestretchcollagen = params_->prestretchcollagen_;
  double density = params_->density_;
  int sizehistory = history_->size();
  double eps = 1.0e-12;
  if (idfiber != 3)
    visrefmassdens_->at(gp)(idfiber) = 0.0;

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
    history_->at(idpast).GetTime(&deptime, &depdt);
    // right dt for explicit integration & adaptation for special implicit step
    if (firstiter == 1)
    {
      double timeloc = 0.0;
      double dtloc = 0.0;
      history_->at(idpast+1).GetTime(&timeloc, &dtloc);
      depdt = dtloc;
    } else if (deptime >= (time - params_->lifetime_ - eps) &&
               (deptime - depdt) < (time - params_->lifetime_ - eps) &&
               (params_->degoption_ != "Exp"))
    {
      depdt = deptime - time + params_->lifetime_;
    }
    LINALG::Matrix<4,1> collstretch(true);
    history_->at(idpast).GetStretches(gp, &collstretch);
    double stretch = prestretchcollagen / collstretch(idfiber);

    // prestretch of collagen fibers is not apllied, might be reasonable combined with prestress
    if (*params_->initstretch_ == "experimental" && deptime <= params_->starttime_ + eps)
      stretch = 1.0 / collstretch(idfiber);

    LINALG::Matrix<NUM_STRESS_3D,1> Saniso(true);
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmataniso(true);
    EvaluateSingleFiber(glstrain, & cmataniso, & Saniso, a, stretch);
    double qdegrad = 0.0;
    Degradation(time - deptime, &qdegrad);
    LINALG::Matrix<4,1> collmass(true);
    history_->at(idpast).GetMass(gp, &collmass);
    double facS = stretch*stretch * qdegrad * collmass(idfiber) / density * depdt;
    Saniso.Scale(facS);
    (*stress) += Saniso;

    if (*params_->initstretch_ == "Homeo" || (idpast == sizehistory - 1 && time > params_->starttime_ + 1.0e-12))
    {
      LINALG::Matrix<NUM_STRESS_3D,1>  A;
      for (int i = 0; i < 3; i++)
        A(i) = a(i)*a(i);

      A(3) = a(0)*a(1); A(4) = a(1)*a(2); A(5) = a(0)*a(2);
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatanisoadd(true);
      cmatanisoadd.MultiplyNT(A, Saniso);
      cmatanisoadd.Scale(-2.0/(collstretch(idfiber)*collstretch(idfiber)));
      (*cmat) += cmatanisoadd;
    } else {
      double faccmat = stretch*stretch*stretch*stretch * qdegrad * collmass(idfiber) / density * depdt;
      cmataniso.Scale(faccmat);
      (*cmat) += cmataniso;
    }

    (*currmassdens) += qdegrad * collmass(idfiber) * depdt;
    if (idfiber != 3)
      visrefmassdens_->at(gp)(idfiber) += qdegrad * collmass(idfiber) * depdt;
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
  for (int i = 0; i < 3; i++)
    A(i) = a(i)*a(i);

  A(3) = a(0)*a(1); A(4) = a(1)*a(2); A(5) = a(0)*a(2);

  double I4 =  A(0)*C(0) + A(1)*C(1) + A(2)*C(2)
             + 1.*(A(3)*C(3) + A(4)*C(4) + A(5)*C(5)); // I4 = trace(A C)
  I4 = I4 * stretch*stretch;  // account for prestretch and stretch at deposition time

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
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmataniso;
    for (int i = 0; i < 6; ++i) {
      for (int j = 0; j < 6; ++j) {
          cmataniso(i,j) = delta7 * A(i) * A(j);  // delta7 A x A
      }
    }
    *cmat = cmataniso;
  }
}

/*----------------------------------------------------------------------*
 |  EvaluateElastin                               (private)        12/10|
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
 |  EvaluateMuscle                                (private)        12/11|
 *----------------------------------------------------------------------*
 strain energy function

 Psi    = k1^m/(2.0*k2^m)*(exp(k2^m*(I_4 - 1.0)^2)-1.0)

        + S_{max} * (\lambda + 1/3*(\lambda_M-\lambda)^3/(\lambda_M-\lambda_0)^2

 I_4 .. invariant accounting for the fiber direction

*/
void MAT::ConstraintMixture::EvaluateMuscle
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress,
  const int gp,
  double* currmassdens
)
{
  const double k1 = params_->k1muscle_; //14175.0;
  const double k2 = params_->k2muscle_; //8.5;
  const double prestretchmuscle = params_->prestretchmuscle_; //1.2;
  const double massfrac = params_->phimuscle_; //0.2;

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
  LINALG::Matrix<3,1>  a = a1_->at(gp);
  LINALG::Matrix<NUM_STRESS_3D,1>  A;
  for (int i = 0; i < 3; i++)
    A(i) = a(i)*a(i);

  A(3) = a(0)*a(1); A(4) = a(1)*a(2); A(5) = a(0)*a(2);

  double I4 =  A(0)*C(0) + A(1)*C(1) + A(2)*C(2)
             + 1.*(A(3)*C(3) + A(4)*C(4) + A(5)*C(5)); // I4 = trace(A C)

  //++++++++++++++++++++++ passive part ++++++++++++++++++
  //--- determine 2nd Piola Kirchhoff stresses S -----------------------------------------
  double preI4 = I4 * prestretchmuscle*prestretchmuscle;
  LINALG::Matrix<NUM_STRESS_3D,1> Saniso(A); // first compute S = 2 dW/dI4 A
  const double exp1 = exp(k2*(preI4-1.)*(preI4-1.));
  if (isinf(exp1))
    dserror("stretch in fiber direction is too high");
  const double fib1 = 2.*(k1*(preI4-1.)*exp1);  // 2 dW/dI4
  Saniso.Scale(fib1);  //S

  // consider mass fraction and prestretch
  Saniso.Scale(massfrac * prestretchmuscle*prestretchmuscle);

  *stress = Saniso;

  //--- do elasticity matrix -------------------------------------------------------------
  const double delta7 = 4.0*(k1*exp1 + 2.0*k1*k2*(preI4-1.0)*(preI4-1.0)*exp1); // 4 d^2Wf/dI4dI4
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmataniso;
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
        cmataniso(i,j) = delta7 * A(i) * A(j);  // delta7 A x A
    }
  }
  // consider mass fraction and prestretch
  cmataniso.Scale(massfrac * prestretchmuscle*prestretchmuscle*prestretchmuscle*prestretchmuscle);

  *cmat = cmataniso;

  //++++++++++++++++++++++ active part ++++++++++++++++++
  // does not depend on mass fraction
  double lambda = sqrt(I4);
  double Smax = params_->Smax_; //50000;
  double lambda_M = 1.2;
  double lambda_0 = 0.7;
  // perhaps check if lambda is between lambda_M and lambda_0

  Saniso.Update(A);
  double facS = Smax / lambda * (1.0 - ((lambda_M-lambda)*(lambda_M-lambda)/(lambda_M-lambda_0)/(lambda_M-lambda_0)));
  Saniso.Scale(facS);
  *stress += Saniso;

  cmataniso.Scale(0.0);
  double faccmat = - Smax / (lambda*lambda) * ((1.0 - ((lambda_M-lambda)*(lambda_M-lambda)/(lambda_M-lambda_0)/(lambda_M-lambda_0)))/lambda
                 - 2.0 * (lambda_M-lambda)/((lambda_M-lambda_0)*(lambda_M-lambda_0)));
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
        cmataniso(i,j) = faccmat * A(i) * A(j);
    }
  }
  *cmat += cmataniso;

  *currmassdens += massfrac * params_->density_;
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
  for (int i = 0; i < 6; i++)
    Svol(i) = J*p * Cinv(i);

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
driving force depends on S

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
  history_->back().GetStretches(gp, &currentstretch);

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
    (*massprodcomp)(2) = massprodbasal_ * (1.0 + params_->growthfactor_ * (massstress3 / params_->homstress_ - 1.0)); // /4.0 massstress3
  else
    (*massprodcomp)(2) = massprodbasal_ * (1.0 + params_->growthfactor_ * massstress3);
//  (*massprodcomp)(2) = 4. * (*massprodcomp)(2);

  double massstress4 = 0.0;
  for (int i=0; i < 3; i++){
    for (int j=0; j < 3; j++){
      massstress4 += a4_->at(gp)(i) * temp2(i,j) * a4_->at(gp)(j);
    }
  }
  massstress4 = sqrt(massstress4) / currentstretch(3);
  if (params_->homstress_ != 0.0)
    (*massprodcomp)(3) = massprodbasal_ * (1.0 + params_->growthfactor_ * (massstress4 / params_->homstress_ - 1.0)); // /4.0 massstress4
  else
    (*massprodcomp)(3) = massprodbasal_ * (1.0 + params_->growthfactor_ * massstress4);
//  (*massprodcomp)(3) = 4. * (*massprodcomp)(3);

  (*massstress)(0) = massstress1;
  (*massstress)(1) = massstress2;
  (*massstress)(2) = massstress3;
  (*massstress)(3) = massstress4;

}

/*----------------------------------------------------------------------*
 |  MassProductionSingleFiber                     (private)        05/11|
 *----------------------------------------------------------------------*
compute new deposition rate for one fiber family
driving force depends on S^k

m^k = mbasal * (1 + K * ( sqrt(a_i^T*C*S^k*C*S^k*C*a_i/detC)/lambda_k(t) /homeo - 1))

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
  history_->back().GetStretches(gp, &currentstretch);

  (*massstress) = temp(0);
  (*massstress) = sqrt(*massstress) / currentstretch(idfiber);
  double homstress = params_->homstress_;
//  if (idfiber == 2 || idfiber == 3)
//    homstress = 4. * homstress;
  if (params_->homstress_ != 0.0)
    (*massprodcomp) = massprodbasal_ * (1.0 + params_->growthfactor_ * ((*massstress) / homstress - 1.0));
  else
    (*massprodcomp) = massprodbasal_ * (1.0 + params_->growthfactor_ * (*massstress));
//  if (idfiber == 2 || idfiber == 3)
//    (*massprodcomp) = 4. * (*massprodcomp);
}

/*----------------------------------------------------------------------*
 |  Degradation                                   (private)        10/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::Degradation(double t, double* degr)
{
  if (params_->degoption_ == "Lin")  // heaviside step function
  {
    if (t <= params_->lifetime_ + 1.0e-12) {
      *degr = 1.0;
    } else {
      *degr = 0.0;
    }
  }
  else if (params_->degoption_ == "Exp")  // exponential decrease
  {
    *degr = exp(- t / params_->lifetime_);
  }
  else if (params_->degoption_ == "Cos")  // transition zone with cos shape
  {
    if (t < 0.2 * params_->lifetime_ - 1.0e-12) {
      *degr = 1.0;
    } else if (t < params_->lifetime_ - 1.0e-12) {
      *degr = 0.5*(cos(PI*(t-0.2*params_->lifetime_)/(0.8*params_->lifetime_))+1.0);
    } else {
      *degr = 0.0;
    }
  } else dserror("Degradation option not implemented! Valid options are Lin, Exp and Cos !");
}

/*----------------------------------------------------------------------*
 |  EvaluateImplicitAll                           (private)        12/11|
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
  LINALG::Matrix<4,1> massprod,
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
  history_->back().GetStretches(gp, &actcollstretch);

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

  history_->back().SetMass(gp,massprod);
  LINALG::Matrix<NUM_STRESS_3D,1> stressresidual(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmattemp(true);
  EvaluateStress(glstrain, gp, &cmattemp, &stressresidual, firstiter, time);

  // determine residual
  LINALG::Matrix<NUM_STRESS_3D,1> Residual(*stress);
  Residual.Update(-1.0,stressresidual,1.0);


  //--------------------------------------------------------------------------------------
  // local Newton iteration
  int localistep = 0;
  int maxstep = 50;
  while (Residual.Norm2() > params_->abstol_*(*stress).NormInf() && localistep < maxstep)
  {
    localistep += 1;
    //--------------------------------------------------------------------------------------
    // derivative of residual
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> DResidual(true);
    for (unsigned int id = 0; id < NUM_STRESS_3D; id++)
      DResidual(id,id) = 1.0;

    // for all 4 fiber families
    double stretch = prestretchcollagen/actcollstretch(0);
    LINALG::Matrix<NUM_STRESS_3D,1> dstressdmass(true);
    LINALG::Matrix<NUM_STRESS_3D,1> dmassdstress(true);
    GradStressDMass(glstrain, &dstressdmass, Cinv, a1_->at(gp), stretch, J, dt, true);
    GradMassDStress(&dmassdstress, Cmatrix, Smatrix, a1_->at(gp), J, massstress(0), actcollstretch(0));
    DResidual.MultiplyNT(-1.0,dstressdmass,dmassdstress,1.0);

    stretch = prestretchcollagen/actcollstretch(1);
    dstressdmass.Scale(0.0);
    dmassdstress.Scale(0.0);
    GradStressDMass(glstrain, &dstressdmass, Cinv, a2_->at(gp), stretch, J, dt, true);
    GradMassDStress(&dmassdstress, Cmatrix, Smatrix, a2_->at(gp), J, massstress(1), actcollstretch(1));
    DResidual.MultiplyNT(-1.0,dstressdmass,dmassdstress,1.0);

    stretch = prestretchcollagen/actcollstretch(2);
    dstressdmass.Scale(0.0);
    dmassdstress.Scale(0.0);
    GradStressDMass(glstrain, &dstressdmass, Cinv, a3_->at(gp), stretch, J, dt, true);
    GradMassDStress(&dmassdstress, Cmatrix, Smatrix, a3_->at(gp), J, massstress(2), actcollstretch(2));
    DResidual.MultiplyNT(-1.0,dstressdmass,dmassdstress,1.0);

    stretch = prestretchcollagen/actcollstretch(3);
    dstressdmass.Scale(0.0);
    dmassdstress.Scale(0.0);
    GradStressDMass(glstrain, &dstressdmass, Cinv, a4_->at(gp), stretch, J, dt, true);
    GradMassDStress(&dmassdstress, Cmatrix, Smatrix, a4_->at(gp), J, massstress(3), actcollstretch(3));
    DResidual.MultiplyNT(-1.0,dstressdmass,dmassdstress,1.0);

    //----------------------------------------------------
    // solve linear system of equations: gradF * incr = -F
    //----------------------------------------------------
    // F = F*-1.0
    Residual.Scale(-1.0);
    LINALG::Matrix<NUM_STRESS_3D,1> increment(true);
    // solve A.X=B
    LINALG::FixedSizeSerialDenseSolver<NUM_STRESS_3D,NUM_STRESS_3D,1> solver;
    solver.SetMatrix(DResidual);              // set A=DResidual
    solver.SetVectors(increment, Residual);           // set X=increment, B=Residual
    solver.FactorWithEquilibration(true); // "some easy type of preconditioning" (Michael)
    int err2 = solver.Factor();           // ?
    int err = solver.Solve();             // X = A^-1 B
    if ((err!=0) || (err2!=0))
      dserror("solving linear system in Newton-Raphson method for implicit integration failed");

    // damping strategy
    double omega = 2.0;
    LINALG::Matrix<NUM_STRESS_3D,1> stepstress(true);
    LINALG::Matrix<NUM_STRESS_3D,1> Residualtemp(Residual);
    double omegamin = 1.0/64.0;
    while (Residualtemp.Norm2() > (1.0-0.5*omega)*Residual.Norm2() && omega > omegamin)
    {
      // update of stress and mass
      omega = omega/2.0;
      stepstress.Update(1.0,*stress,omega,increment);

      MassProduction(gp, C, stepstress, &massstress, &massprod);
      history_->back().SetMass(gp, massprod);
      stressresidual.Scale(0.0);
      EvaluateStress(glstrain, gp, &cmattemp, &stressresidual, firstiter, time);

      Residualtemp.Update(1.0,stepstress,-1.0,stressresidual);
    }
    if (omega <= omegamin && Residualtemp.Norm2() > (1.0-0.5*omega)*Residual.Norm2())
      dserror("no damping coefficient found");

    *stress = stepstress;
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
  if (localistep == maxstep && Residual.Norm2() > params_->abstol_*(*stress).NormInf())
    dserror("local Newton iteration did not converge %e", Residual.Norm2());

  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatelastic(true);
  EvaluateStress(glstrain, gp, &cmatelastic, stress, firstiter, time);

  //--------------------------------------------------------------------------------------
  // compute cmat
  // right handside of the linear equations
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> RHS(cmatelastic);
  // left matrix of the linear equations
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> LM(true);
  for (unsigned int id = 0; id < NUM_STRESS_3D; id++) LM(id,id) = 1.0;

  // Fiber1
  double stretch = prestretchcollagen/actcollstretch(0);
  LINALG::Matrix<NUM_STRESS_3D,1> dmassdstress(true);
  LINALG::Matrix<NUM_STRESS_3D,1> dmassdstretch(true);
  LINALG::Matrix<NUM_STRESS_3D,1> dstressdmass(true);
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
  LINALG::FixedSizeSerialDenseSolver<NUM_STRESS_3D,NUM_STRESS_3D,NUM_STRESS_3D> solver;
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
 |  EvaluateImplicitSingle                        (private)        11/11|
 *----------------------------------------------------------------------*
evaluate stress and cmat for implicit integration
driving force of massproduction is the fiber stress S^k
Newton loop only for stresses
*/
void MAT::ConstraintMixture::EvaluateImplicitSingle
(
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  const int gp,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress,
  double dt,
  double time
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
  Degradation(0.0, &qdegrad);
  double density = params_->density_;
  // store actual collagen stretches, do not change anymore
  LINALG::Matrix<4,1> actcollstretch(true);
  history_->back().GetStretches(gp, &actcollstretch);
  LINALG::Matrix<4,1> massprod(true);
  history_->back().GetMass(gp, &massprod);
  LINALG::Matrix<4,1> massstress(true);

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

    double stretch = prestretchcollagen/actcollstretch(idfiber);
    LINALG::Matrix<NUM_STRESS_3D,1> stressfiber(true);
    LINALG::Matrix<NUM_STRESS_3D,1> stressresidual(true);
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatfiber(true);
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmattemp(true);
    double currmassdensfiber = 0.0;
    double currmassdenstemp = 0.0;
    double massprodfiber = 0.0;
    double massstressfiber = 0.0;
    EvaluateFiberFamily(glstrain, gp, &cmatfiber, &stressfiber, a, &currmassdensfiber, firstiter, time, idfiber);
    // mass always corresponds to the current stress
    MassProductionSingleFiber(gp, C, stressfiber, &massstressfiber, &massprodfiber, a, idfiber);
    massprod(idfiber) = massprodfiber;
    history_->back().SetMass(gp, massprod);
    massstress(idfiber) = massstressfiber;
    // compute stresses for the computed mass
    EvaluateFiberFamily(glstrain, gp, &cmattemp, &stressresidual, a, &currmassdenstemp, firstiter, time, idfiber);

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
    LINALG::Matrix<NUM_STRESS_3D,1> Residual(stressfiber);
    Residual.Update(-1.0,stressresidual,1.0);

    //--------------------------------------------------------------------------------------
    // local Newton iteration to determine stress
    int localistep = 0;
    int maxstep = 50;
    while (Residual.Norm2() > params_->abstol_*stressfiber.NormInf() && localistep < maxstep)
    {
      localistep += 1;
      //--------------------------------------------------------------------------------------
      // derivative of residual
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> DResidual(true);
      for (unsigned int id = 0; id < NUM_STRESS_3D; id++)
        DResidual(id,id) = 1.0;

      // linearisation of stress formula
      LINALG::Matrix<NUM_STRESS_3D,1> dstressdmass(true);
      GradStressDMass(glstrain, &dstressdmass, Cinv, a, stretch, J, dt, false);
      LINALG::Matrix<NUM_STRESS_3D,1> dmassdstress(true);
      GradMassDStress(&dmassdstress, Cmatrix, Smatrix, a, J, massstress(idfiber), actcollstretch(idfiber));
      DResidual.MultiplyNT(-1.0,dstressdmass,dmassdstress,1.0);

      //----------------------------------------------------
      // solve linear system of equations: gradF * incr = -F
      //----------------------------------------------------
      // F = F*-1.0
      Residual.Scale(-1.0);
      LINALG::Matrix<NUM_STRESS_3D,1> increment(true);
      // solve A.X=B
      LINALG::FixedSizeSerialDenseSolver<NUM_STRESS_3D,NUM_STRESS_3D,1> solver;
      solver.SetMatrix(DResidual);              // set A=DResidual
      solver.SetVectors(increment, Residual);           // set X=increment, B=Residual
      solver.FactorWithEquilibration(true); // "some easy type of preconditioning" (Michael)
      int err2 = solver.Factor();           // ?
      int err = solver.Solve();             // X = A^-1 B
      if ((err!=0) || (err2!=0))
        dserror("solving linear system in Newton-Raphson method for implicit integration failed");

      // damping strategy
      double omega = 2.0;
      LINALG::Matrix<NUM_STRESS_3D,1> stepstress(true);
      LINALG::Matrix<NUM_STRESS_3D,1> Residualtemp(Residual);
      double omegamin = 1.0/64.0;
      while (Residualtemp.Norm2() > (1.0-0.5*omega)*Residual.Norm2() && omega > omegamin)
      {
        // update of stress
        omega = omega/2.0;
        stepstress.Update(1.0,stressfiber,omega,increment);

        // corresponding mass
        MassProductionSingleFiber(gp, C, stepstress, &massstressfiber, &massprodfiber, a, idfiber);
        massprod(idfiber) = massprodfiber;
        history_->back().SetMass(gp, massprod);
        massstress(idfiber) = massstressfiber;
        // compute stresses for the computed mass
        stressresidual.Scale(0.0);
        EvaluateFiberFamily(glstrain, gp, &cmattemp, &stressresidual, a, &currmassdenstemp, firstiter, time, idfiber);

        Residualtemp.Update(1.0,stepstress,-1.0,stressresidual);
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
    if (localistep == maxstep && Residual.Norm2() > params_->abstol_*stressfiber.NormInf())
      dserror("local Newton iteration did not converge %e", Residual.Norm2());

    currmassdensfiber = 0.0;
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatelastic(true);
    stressfiber.Scale(0.0);
    EvaluateFiberFamily(glstrain, gp, &cmatelastic, &stressfiber, a, &currmassdensfiber, firstiter, time, idfiber);

    //--------------------------------------------------------------------------------------
    // compute cmat
    // right handside of the linear equations
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> RHS(cmatelastic);
    // left matrix of the linear equations
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> LM(true);
    for (unsigned int id = 0; id < NUM_STRESS_3D; id++) LM(id,id) = 1.0;

    LINALG::Matrix<NUM_STRESS_3D,1> dmassdstress(true);
    LINALG::Matrix<NUM_STRESS_3D,1> dmassdstretch(true);
    LINALG::Matrix<NUM_STRESS_3D,1> dstressdmass(true);
    GradMassDStretch(&dmassdstretch, Cmatrix, Smatrix, Cinv, a, J, massstress(idfiber), actcollstretch(idfiber), dt);
    GradMassDStress(&dmassdstress, Cmatrix, Smatrix, a, J, massstress(idfiber), actcollstretch(idfiber));
    GradStressDMass(glstrain, &dstressdmass, Cinv, a, stretch, J, dt, false);
    RHS.MultiplyNT(2.0,dstressdmass,dmassdstretch,1.0);
    LM.MultiplyNT(-1.0,dstressdmass,dmassdstress,1.0);

    cmatfiber.Scale(0.0);
    //----------------------------------------------------
    // solve linear system of equations: A.X=B
    //----------------------------------------------------
    LINALG::FixedSizeSerialDenseSolver<NUM_STRESS_3D,NUM_STRESS_3D,NUM_STRESS_3D> solver;
    solver.SetMatrix(LM);              // set A=LM
    solver.SetVectors(cmatfiber, RHS);           // set X=increment, B=RHS
    solver.FactorWithEquilibration(true); // "some easy type of preconditioning" (Michael)
    int err2 = solver.Factor();           // ?
    int err = solver.Solve();             // X = A^-1 B
    if ((err!=0) || (err2!=0))
      dserror("solving linear system for cmat failed");

    (*stress) += stressfiber;
    (*cmat) += cmatfiber;
    currmassdens += currmassdensfiber;

    // volumetric part, that is related to this fiber family
    LINALG::Matrix<NUM_STRESS_3D,1> tempvol(true);
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatvolfiber(true);
    tempvol.MultiplyTN(cmatfiber, dmassdstress);
    tempvol.Update(2.0,dmassdstretch,1.0);
    cmatvolfiber.MultiplyNT(-1.0*dt/density*qdegrad*params_->kappa_*J,Cinv,tempvol);
    //(*cmat).MultiplyNT(-2.0*dt/density*qdegrad*params_->kappa_*J,Cinv,dmassdstretch,1.0);
    (*cmat) += cmatvolfiber;
  }

  // elastin
  double refmassdenselastin = params_->phielastin_ * density;
  double prestretchelastin = params_->prestretchelastin_;
  // account for isotropic prestretch of elastin
  LINALG::Matrix<NUM_STRESS_3D,1> glstrainiso(*glstrain);
  glstrainiso.Scale(prestretchelastin*prestretchelastin);
  LINALG::Matrix<NUM_STRESS_3D,1> Siso(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatiso(true);
  EvaluateElastin(&glstrainiso, &cmatiso, &Siso);
  Siso.Scale(refmassdenselastin/density*prestretchelastin*prestretchelastin);
  (*stress) += Siso;
  cmatiso.Scale(refmassdenselastin/density*prestretchelastin*prestretchelastin*prestretchelastin*prestretchelastin);
  (*cmat) += cmatiso;
  currmassdens += refmassdenselastin;

  // smooth muscle cells
  LINALG::Matrix<NUM_STRESS_3D,1> Smus(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatmus(true);
  EvaluateMuscle(glstrain, &cmatmus, &Smus, gp, &currmassdens);
  (*stress) += Smus;
  (*cmat) += cmatmus;

  // volumetric part
  LINALG::Matrix<NUM_STRESS_3D,1> Svol(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatvol(true);
  EvaluateVolumetric(glstrain, &cmatvol, &Svol, currmassdens, density);
  (*stress) += Svol;
  (*cmat) += cmatvol;

  vismassstress_->at(gp)(0) = massstress(0);
  vismassstress_->at(gp)(1) = massstress(1);
  vismassstress_->at(gp)(2) = massstress(2);
  refmassdens_->at(gp) = currmassdens;
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
  double facS = stretch*stretch * qdegrad / density * dt;
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
  double fac = massprodbasal_ * params_->growthfactor_ / homstressfixed / (J*J)
             / massstress / (actcollstretch*actcollstretch);
  if (massstress == 0.0) fac = 0.0;
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
  temp.Scale(0.0);
  temp.Multiply(Cmatrix,SCa);
  SCSCa.Multiply(Smatrix,temp);

  (*derivative)(0) = a(0)*a(0); (*derivative)(1) = a(1)*a(1);
  (*derivative)(2) = a(2)*a(2); (*derivative)(3) = a(0)*a(1);
  (*derivative)(4) = a(1)*a(2); (*derivative)(5) = a(0)*a(2);
  (*derivative).Update(1.0,Cinv,1.0/(actcollstretch*actcollstretch));
  (*derivative).Scale(- 1.0 * massstress);
  double fac = 1.0 / massstress / (actcollstretch*actcollstretch) / (J*J);
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
// History
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
  AddtoPack(data, type);

  // Pack internal variables
  AddtoPack(data, depositiontime_);
  AddtoPack(data, dt_);
  AddtoPack(data, numgp_);
  for (int gp = 0; gp < numgp_; ++gp)
  {
    AddtoPack(data, collagenstretch1_->at(gp));
    AddtoPack(data, collagenstretch2_->at(gp));
    AddtoPack(data, collagenstretch3_->at(gp));
    AddtoPack(data, collagenstretch4_->at(gp));
    AddtoPack(data, massprod1_->at(gp));
    AddtoPack(data, massprod2_->at(gp));
    AddtoPack(data, massprod3_->at(gp));
    AddtoPack(data, massprod4_->at(gp));
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
  ExtractfromPack(position, data, a);
  depositiontime_ = a;
  ExtractfromPack(position, data, a);
  dt_ = a;
  int b;
  ExtractfromPack(position, data, b);
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
    ExtractfromPack(position, data, a);
    collagenstretch1_->at(gp) = a;
    ExtractfromPack(position, data, a);
    collagenstretch2_->at(gp) = a;
    ExtractfromPack(position, data, a);
    collagenstretch3_->at(gp) = a;
    ExtractfromPack(position, data, a);
    collagenstretch4_->at(gp) = a;
    ExtractfromPack(position, data, a);
    massprod1_->at(gp) = a;
    ExtractfromPack(position, data, a);
    massprod2_->at(gp) = a;
    ExtractfromPack(position, data, a);
    massprod3_->at(gp) = a;
    ExtractfromPack(position, data, a);
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
    massprod3_->at(gp) = massprodbasal; //*4.;
    massprod4_->at(gp) = massprodbasal; //*4.;
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
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

    Teuchos::RefCountPtr<MAT::Material> mat = actele->Material();
    MAT::ConstraintMixture* grow = static_cast <MAT::ConstraintMixture*>(mat.get());

    // material plot at gauss points
    int ngp = grow->Geta1()->size();

    // update element geometry
    const int numnode = actele->NumNode();
    const int numdof = 3;
    Epetra_SerialDenseMatrix xcurr(numnode,3);  // material coord. of element
    for (int i=0; i<numnode; ++i)
    {
      xcurr(i,0) = actele->Nodes()[i]->X()[0]+ mydisp[i*numdof+0];
      xcurr(i,1) = actele->Nodes()[i]->X()[1]+ mydisp[i*numdof+1];
      xcurr(i,2) = actele->Nodes()[i]->X()[2]+ mydisp[i*numdof+2];
    }
    const DRT::Element::DiscretizationType distype = actele->Shape();
    Epetra_SerialDenseVector funct(numnode);

    // define gauss rule
    DRT::UTILS::GaussRule3D gaussrule_ = DRT::UTILS::intrule3D_undefined;
    switch (distype)
    {
    case DRT::Element::hex8:
    {
      gaussrule_ = DRT::UTILS::intrule_hex_8point;
      if (ngp != 8)
        dserror("hex8 has not 8 gauss points: %d", ngp);
      break;
    }
    case DRT::Element::wedge6:
    {
      gaussrule_ = DRT::UTILS::intrule_wedge_6point;
      if (ngp != 6)
        dserror("wedge6 has not 6 gauss points: %d", ngp);
      break;
    }
    case DRT::Element::tet4:
    {
      gaussrule_ = DRT::UTILS::intrule_tet_1point;
      if (ngp != 1)
        dserror("tet4 has not 1 gauss point: %d", ngp);
      break;
    }
    default:
      dserror("unknown element in ConstraintMixtureOutputToGmsh");
    }

    const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule_);

    for (int gp = 0; gp < ngp; ++gp)
    {
      DRT::UTILS::shape_function_3D(funct, intpoints.qxg[gp][0], intpoints.qxg[gp][1], intpoints.qxg[gp][2], distype);
      Epetra_SerialDenseMatrix point(1,3);
      point.Multiply('T','N',1.0,funct,xcurr,0.0);

      // write mandel stress
      LINALG::Matrix<3,1> mandelgp = grow->GetVis(gp);
      gmshfilecontent << "SP(" << scientific << point(0,0) << ",";
      gmshfilecontent << scientific << point(0,1) << ",";
      gmshfilecontent << scientific << point(0,2) << ")";
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

