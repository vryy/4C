/*----------------------------------------------------------------------------*/
/*!
\file constraintmixture.cpp
\brief This file contains routines for constraint mixture growth and remodeling.
example input line
MAT 1 MAT_ConstraintMixture DENS 0.001 MUE 1.0 NUE 0.49 PHIE 0.08
PREELA 1.0 K1 1.0 K2 1.0 NUMHOM 1 PRECOLL 1.06 DAMAGE 1.5 K1M 1.0 K2M 1.0
PHIM 1.0 PREMUS 1.0 SMAX 0.0 KAPPA 1.0E6 LIFETIME 5.0 GROWTHFAC 0.5
HOMSTR 6.75E4 SHEARGROWTHFAC 0.0 HOMRAD 1.0 STARTTIME 5.0
INTEGRATION Explicit TOL 1.0E-4 GROWTHFORCE Single ELASTINDEGRAD None
MASSPROD Lin INITSTRETCH None CURVE 0 DEGOPTION Cos MAXMASSPRODFAC 4.0
STOREHISTORY No

Here an approach for growth and remodeling of an artery is modeled.
For a detailed description see:
- Humphrey, J. & Rajagopal, K.: A constrained mixture model for arterial
  adaptations to a sustained step change in blood flow,
  Biomechanics and Modeling in Mechanobiology, 2003, 2, 109-126

\level 2

\maintainer Fabian Braeu
            braeu@lnm.mw.tum.de

*/
/*----------------------------------------------------------------------------*/

#include "constraintmixture.H"
#include "constraintmixture_history.H"
#include "../drt_lib/drt_globalproblem.H"
//#include "../drt_lib/standardtypes_cpp.H" // for PI, not needed here?
#include "matpar_bundle.H"
#include "../drt_mat/material_service.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_comm/comm_utils.H"
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
  nue_(matdata->GetDouble("NUE")),
  phielastin_(matdata->GetDouble("PHIE")),
  prestretchelastin_(matdata->GetDouble("PREELA")),
  k1_(matdata->GetDouble("K1")),
  k2_(matdata->GetDouble("K2")),
  numhom_(matdata->GetInt("NUMHOM")),
  prestretchcollagen_(*(matdata->Get<std::vector<double> >("PRECOLL"))),
  damagestretch_(matdata->GetDouble("DAMAGE")),
  k1muscle_(matdata->GetDouble("K1M")),
  k2muscle_(matdata->GetDouble("K2M")),
  phimuscle_(matdata->GetDouble("PHIM")),
  prestretchmuscle_(matdata->GetDouble("PREMUS")),
  Smax_(matdata->GetDouble("SMAX")),
  kappa_(matdata->GetDouble("KAPPA")),
  lifetime_(matdata->GetDouble("LIFETIME")),
  //growthfactor_(matdata->GetDouble("GROWTHFAC")),
  homstress_(*(matdata->Get<std::vector<double> >("HOMSTR"))),
  sheargrowthfactor_(matdata->GetDouble("SHEARGROWTHFAC")),
  homradius_(matdata->GetDouble("HOMRAD")),
  starttime_(matdata->GetDouble("STARTTIME")),
  integration_(matdata->Get<std::string>("INTEGRATION")),
  abstol_(matdata->GetDouble("TOL")),
  growthforce_(matdata->Get<std::string>("GROWTHFORCE")),
  elastindegrad_(matdata->Get<std::string>("ELASTINDEGRAD")),
  massprodfunc_(matdata->Get<std::string>("MASSPROD")),
  initstretch_(matdata->Get<std::string>("INITSTRETCH")),
  timecurve_(matdata->GetInt("CURVE")),
  degoption_(*(matdata->Get<std::string>("DEGOPTION"))),
  maxmassprodfac_(matdata->GetDouble("MAXMASSPRODFAC")),
  storehistory_(matdata->GetInt("STOREHISTORY")),
  degtol_(1.0e-6)
{
  Epetra_Map dummy_map(1,1,0,*(DRT::Problem::Instance()->GetNPGroup()->LocalComm()));
  for(int i=first ; i<=last; i++)
  {
    matparams_.push_back(Teuchos::rcp(new Epetra_Vector(dummy_map,true)));
  }
  matparams_.at(growthfactor)->PutScalar(matdata->GetDouble("GROWTHFAC"));
  matparams_.at(elastin_survival)->PutScalar(matdata->GetDouble("ELASTINFAC"));
}


Teuchos::RCP<MAT::Material> MAT::PAR::ConstraintMixture::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ConstraintMixture(this));
}

void MAT::PAR::ConstraintMixture::OptParams(std::map<std::string,int>* pnames)
{
  pnames->insert(std::pair<std::string,int>("GROWTHFAC", growthfactor));
  pnames->insert(std::pair<std::string,int>("ELASTINFAC", elastin_survival));
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
  for (int gp = 0; gp < numgp; gp++)
  {
    AddtoPack(data, a1_->at(gp));
    AddtoPack(data, a2_->at(gp));
    AddtoPack(data, a3_->at(gp));
    AddtoPack(data, a4_->at(gp));
    AddtoPack(data, vismassstress_->at(gp));
    AddtoPack(data, refmassdens_->at(gp));
    AddtoPack(data, visrefmassdens_->at(gp));
    AddtoPack(data, localprestretch_->at(gp));
    AddtoPack(data, localhomstress_->at(gp));
  }
  if (numgp > 0)
  {
    AddtoPack(data,massprodbasal_);
    AddtoPack(data,homradius_);

    int delsize = deletemass_->size();
    AddtoPack(data,delsize);
    for (int iddel = 0; iddel < delsize; iddel++)
      AddtoPack(data,deletemass_->at(iddel));

    // Pack history
    int minindex = minindex_;
    AddtoPack(data, minindex);
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
void MAT::ConstraintMixture::Unpack(const std::vector<char>& data)
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
  a1_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,1> >(numgp));
  a2_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,1> >(numgp));
  a3_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,1> >(numgp));
  a4_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,1> >(numgp));
  vismassstress_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,1> > (numgp));
  refmassdens_ = Teuchos::rcp(new std::vector<double> (numgp));
  visrefmassdens_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,1> > (numgp));
  localprestretch_ = Teuchos::rcp(new std::vector<LINALG::Matrix<4,1> > (numgp));
  localhomstress_ = Teuchos::rcp(new std::vector<LINALG::Matrix<4,1> > (numgp));

  for (int gp = 0; gp < numgp; gp++)
  {
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
    LINALG::Matrix<4,1> pre;
    ExtractfromPack(position, data, pre);
    localprestretch_->at(gp) = pre;
    ExtractfromPack(position, data, pre);
    localhomstress_->at(gp) = pre;
  }
  double basal;
  ExtractfromPack(position,data,basal);
  massprodbasal_ = basal;
  ExtractfromPack(position,data,basal);
  homradius_ = basal;

  int delsize = 0;
  ExtractfromPack(position,data,delsize);
  deletemass_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,1> > (delsize));
  for (int iddel = 0; iddel < delsize; iddel++)
  {
    LINALG::Matrix<3,1> deldata;
    ExtractfromPack(position,data,deldata);
    deletemass_->at(iddel) = deldata;
  }

  // unpack history
  ExtractfromPack(position,data,delsize);
  minindex_ = delsize;
  int sizehistory;
  ExtractfromPack(position, data, sizehistory);
  history_ = Teuchos::rcp(new std::vector<ConstraintMixtureHistory> (sizehistory));
  for (int idpast = 0; idpast < sizehistory; idpast++)
  {
    std::vector<char> datahistory;
    ExtractfromPack(position, data, datahistory);
    history_->at(idpast).Unpack(datahistory);
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  /*
  double oldesttime = 0.0;
  double acttime = 0.0;
  double tempdt = 0.0;
  history_->begin()->GetTime(&oldesttime, &tempdt);
  history_->at(sizehistory-2).GetTime(&acttime, &tempdt);
  double intdegr = 0.0;
  double degrtime = 0.0;
  double degrdt = 0.0;
  for (int idpast = 0; idpast < sizehistory -1; idpast++)
  {
    double degr = 0.0;
    history_->at(idpast).GetTime(&degrtime, &degrdt);
    double timeloc = 0.0;
    double dtloc = 0.0;
    history_->at(idpast+1).GetTime(&timeloc, &dtloc);
    degrdt = dtloc;
    Degradation(acttime-degrtime, degr);
    intdegr += degr * degrdt;
  }

  for (int idpast = 0; idpast < sizehistory-1; idpast++)
  {
    double temptime = 0.0;
    history_->at(idpast).GetTime(&temptime, &tempdt);
    for (int idgauss = 0; idgauss < numgp; idgauss++)
    {
      LINALG::Matrix<4,1> stretchtemp(true);
      LINALG::Matrix<4,1> stretchact(true);
      LINALG::Matrix<4,1> stretchold(true);
      history_->at(sizehistory-2).GetStretches(idgauss, &stretchact);
      history_->begin()->GetStretches(idgauss, &stretchold);
      stretchtemp.Update(stretchact);
      // linear interpolated stretch
      //double scalar = (temptime - acttime) / (oldesttime - acttime);
      //stretchtemp.Update(scalar,stretchold,1.0);
      //stretchtemp.Update(-scalar,stretchact,1.0);
      // modify stretch
      //history_->at(idpast).SetStretches(idgauss,stretchtemp);
      // distribute mass equally
      double massprodbasal = (refmassdens_->at(idgauss) - (params_->phimuscle_ + params_->phielastin_) * params_->density_) / 4.0 / intdegr;
      LINALG::Matrix<4,1> masstemp(true);
      masstemp.PutScalar(massprodbasal);
      history_->at(idpast).SetMass(idgauss,masstemp);
    }
  }
  std::cout << "Unpack called, history of mass/stretch is lost" << std::endl;
  */

  return;
}

/*----------------------------------------------------------------------*
 |  Setup                                         (public)         12/10|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  if (*params_->integration_ != "Implicit" && *params_->integration_ != "Explicit")
    dserror("unknown option for integration");
  if (*params_->growthforce_ != "Single" && *params_->growthforce_ != "All" && *params_->growthforce_ != "ElaCol")
    dserror("unknown driving force for growth");
  if (*params_->elastindegrad_ != "None" && *params_->elastindegrad_ != "Rectangle" && *params_->elastindegrad_ != "Time"
      && *params_->elastindegrad_ != "RectanglePlate" && *params_->elastindegrad_ != "Wedge"
      && *params_->elastindegrad_ != "Circles" && *params_->elastindegrad_ != "InvEla")
    dserror("unknown option for elastin degradation");
  if (*params_->massprodfunc_ != "Lin" && *params_->massprodfunc_ != "CosCos")
    dserror("unknown option for mass production function");

  // visualization
  vismassstress_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,1> > (numgp));
  refmassdens_ = Teuchos::rcp(new std::vector<double> (numgp));
  visrefmassdens_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,1> > (numgp));
  // homeostatic prestretch of collagen fibers
  localprestretch_ = Teuchos::rcp(new std::vector<LINALG::Matrix<4,1> > (numgp));
  localhomstress_ = Teuchos::rcp(new std::vector<LINALG::Matrix<4,1> > (numgp));

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

  deletemass_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,1> > (0));

  // history
  ResetAll(numgp);

  // fiber vectors
  a1_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,1> > (numgp));
  a2_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,1> > (numgp));
  a3_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,1> > (numgp));
  a4_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,1> > (numgp));

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

  const double gamma = (45.0*PI)/180.; //angle for diagonal fibers

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
 |  ResetAll                                      (public)         03/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::ResetAll(const int numgp)
{
  // homeostatic variables
  for (int gp = 0; gp < numgp; gp++)
  {
    if (params_->numhom_ == 1)
    {
      localprestretch_->at(gp).PutScalar(params_->prestretchcollagen_[0]);
      localhomstress_->at(gp).PutScalar(params_->homstress_[0]);
    }
    else if (params_->numhom_ == 3)
    {
      localprestretch_->at(gp)(0) = params_->prestretchcollagen_[0];
      localprestretch_->at(gp)(1) = params_->prestretchcollagen_[1];
      localprestretch_->at(gp)(2) = params_->prestretchcollagen_[2];
      localprestretch_->at(gp)(3) = params_->prestretchcollagen_[2];
      localhomstress_->at(gp)(0) = params_->homstress_[0];
      localhomstress_->at(gp)(1) = params_->homstress_[1];
      localhomstress_->at(gp)(2) = params_->homstress_[2];
      localhomstress_->at(gp)(3) = params_->homstress_[2];
    }
    else
      dserror("wrong number of homeostatic variables");
  }
  homradius_ = params_->homradius_;

  // history
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
  if (params_->degoption_ == "Exp" || params_->degoption_ == "ExpVar")
  {
    double taumax = - log(params_->degtol_) / log(2.0) * params_->lifetime_;
    numpast = static_cast<int>(round(taumax / dt)) + firstiter;
  }
  minindex_ = 0;

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
    Degradation((numpast-1-idpast)*dt, degr);
    intdegr += degr * dt;
  }
  massprodbasal_ = (1.0 - params_->phimuscle_ - params_->phielastin_) * params_->density_ / 4.0 / intdegr;
//  massprodbasal_ = (1.0 - params_->phimuscle_ - params_->phielastin_) * params_->density_ / 10.0 / intdegr;

  // history
  history_ = Teuchos::rcp(new std::vector<ConstraintMixtureHistory> (numpast));
  bool expvar = false;
  if (params_->degoption_ == "ExpVar") expvar = true;
  for (int idpast = 0; idpast < numpast; idpast++)
  {
    history_->at(idpast).Setup(numgp, massprodbasal_, expvar);
    history_->at(idpast).SetTime(dt-(numpast-1-idpast)*dt, dt);
  }
}

/*----------------------------------------------------------------------*
 |  Update internal variables                     (public)         12/10|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::Update()
{
  // set mass of damaged collagen to zero
  // has to be done before history vector is changed
  // do not erase as they are needed for computation of massprodbasal
  int delsize = deletemass_->size();
  for (int iddel = 0; iddel < delsize; iddel++)
  {
    int gpdel = deletemass_->at(iddel)(0);
    int idpastdel = deletemass_->at(iddel)(1);
    int idfiberdel = deletemass_->at(iddel)(2);
    history_->at(idpastdel).SetMass(gpdel, 0.0, idfiberdel);
  }

  // erase oldest collagen
  int numgp = history_->at(0).NumGP();
  int sizehistory = history_->size();
  double deptime = 0.0;
  double depdt = 0.0;
  history_->back().GetTime(&deptime, &depdt);
  bool expvar = false;
  if (params_->degoption_ == "ExpVar") expvar = true;

  // just do update in case of growth
  if (deptime > params_->starttime_ + 1.0e-11)
  {
    if (!params_->storehistory_)
    {
      // delete just the steps that surely won't be needed, especially with a smaller timestep later
      // thus reference time is deposition time of last collagen fibers
      double acttime = deptime;
      int eraseiter = 0;
      history_->at(eraseiter).GetTime(&deptime, &depdt);
      double degrad = 0.0;
      Degradation(acttime-deptime, degrad);
      while (degrad < params_->degtol_ && eraseiter < sizehistory)
      {
        eraseiter +=1;
        history_->at(eraseiter).GetTime(&deptime, &depdt);
        Degradation(acttime-deptime, degrad);
      }
      if (eraseiter > 0)
      {
        history_->erase(history_->begin(),history_->begin()+eraseiter);
        //std::cout << "erased " << eraseiter << " history variables" << std::endl;
      }
    }

    // append new collagen
    ConstraintMixtureHistory newhis;
    newhis.Setup(numgp, massprodbasal_, expvar);
    // it is very important to set time and dt to 0.0
    // this makes it clear that this step was created in update and has no reliable content
    // they are not known here either
    // in EvaluateStress this is important, if called after Update
    newhis.SetTime(0.0, 0.0);
    history_->push_back(newhis);

  }
  else // time < starttime, adapt deposition time, set history if wanted
  {
    // SetHistory
    if (*params_->initstretch_ == "SetConstantHistory")
    {
      for (int igp=0; igp < numgp; igp++)
      {
        // set stretch
        LINALG::Matrix<4,1> actstretch(true);
        history_->back().GetStretches(igp, &actstretch);
        // special treatment for first timestep
        if (deptime == depdt)
        {
          LINALG::Matrix<4,1> ones(true);
          ones.PutScalar(1.0);
          history_->back().SetStretches(igp, ones);
        }
        for (int istep = 0; istep < sizehistory; istep++)
        {
          LINALG::Matrix<4,1> oldstretch(true);
          history_->at(istep).GetStretches(igp, &oldstretch);
          oldstretch.EDivide(actstretch);
          history_->at(istep).SetStretches(igp, oldstretch);
        }

        // set mass
        LINALG::Matrix<4,1> actmass(true);
        history_->back().GetMass(igp, &actmass);
        if (abs(actmass(0)-massprodbasal_) > 1.0e-8 || abs(actmass(1)-massprodbasal_) > 1.0e-8
            || abs(actmass(2)-massprodbasal_) > 1.0e-8 || abs(actmass(3)-massprodbasal_) > 1.0e-8)
        {
          LINALG::Matrix<4,1> ones(true);
          ones.PutScalar(1.0);
          LINALG::Matrix<1,1> summass(true);
          summass.MultiplyTN(ones,actmass);
          LINALG::Matrix<4,1> newmass(true);
          newmass.PutScalar(4.0*massprodbasal_);
          if (summass(0) > 0.0 + 1.0e-12)
          {
            actmass.Scale(1.0/summass(0));
            newmass.EMultiply(actmass);
          }
          for (int istep = 0; istep < sizehistory; istep++)
          {
            history_->at(istep).SetMass(igp, newmass);
          }
          actmass.PutScalar(massprodbasal_);
          history_->back().SetMass(igp, actmass);
        }
      }
    }

    if (*params_->initstretch_ == "SetLinearHistory")
    {
      for (int igp=0; igp < numgp; igp++)
      {
        LINALG::Matrix<4,1> actstretch(true);
        history_->back().GetStretches(igp, &actstretch);
        LINALG::Matrix<4,1> actmass(true);
        history_->back().GetMass(igp, &actmass);
        if (abs(actmass(0)-massprodbasal_) < 1.0e-8 && abs(actmass(1)-massprodbasal_) < 1.0e-8
            && abs(actmass(2)-massprodbasal_) < 1.0e-8 && abs(actmass(3)-massprodbasal_) < 1.0e-8)
        {
          if (abs(actstretch(0)-1.0) > 1.0e-8 || abs(actstretch(1)-1.0) > 1.0e-8
              || abs(actstretch(2)-1.0) > 1.0e-8 || abs(actstretch(3)-1.0) > 1.0e-8)
          {
            // set stretch
            for (int istep = 0; istep < sizehistory-1; istep++)
            {
              LINALG::Matrix<4,1> oldstretch(true);
              for (int idfiber=0;idfiber<4;idfiber++)
                oldstretch(idfiber) = (actstretch(idfiber)-1.0)*(1.0-istep/(sizehistory-1.0))+1.0;
              LINALG::Matrix<4,1> ones(true);
              ones.PutScalar(1.0);
              ones.EDivide(oldstretch);
              history_->at(istep).SetStretches(igp, ones);
            }
            actstretch.PutScalar(1.0);
            history_->back().SetStretches(igp, actstretch);
          }
        }
        else
        {
          // set mass
          LINALG::Matrix<4,1> ones(true);
          ones.PutScalar(1.0);
          LINALG::Matrix<1,1> summass(true);
          summass.MultiplyTN(ones,actmass);
          LINALG::Matrix<4,1> newmass(true);
          newmass.PutScalar(4.0*massprodbasal_);
          if (summass(0) > 0.0 + 1.0e-12)
          {
            actmass.Scale(1.0/summass(0));
            newmass.EMultiply(actmass);
          }
          for (int istep = 0; istep < sizehistory; istep++)
          {
            history_->at(istep).SetMass(igp, newmass);
          }
          actmass.PutScalar(massprodbasal_);
          history_->back().SetMass(igp, actmass);
        }
      }
    }

    // just adopt deposition time, the rest stays the same
    double newtime = 0.0;
    double newdt = 0.0;
    history_->at(0).GetTime(&newtime, &newdt);
    double degrad = 0.0;
    Degradation(deptime-newtime, degrad);
    if (degrad < params_->degtol_)
    {
      for (int iter = 0; iter < sizehistory-1; iter++)
      {
        history_->at(iter+1).GetTime(&newtime, &newdt);
        history_->at(iter).SetTime(newtime, newdt);
      }
      history_->back().SetTime(0.0, 0.0);
    }
    else if (deptime == depdt || (deptime == 2*depdt && *params_->integration_ == "Implicit"))
    {
      // special case of first time step
      ConstraintMixtureHistory newhis;
      newhis.Setup(numgp, massprodbasal_, expvar);
      newhis.SetTime(0.0, 0.0);
      history_->push_back(newhis);
    }
    else
    {
      dserror("You should not change your timestep size in the case time < starttime! %f", deptime);
    }
  } // time < starttime
}

/*----------------------------------------------------------------------*
 |  Reset internal variables                      (public)         01/12|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::ResetStep()
{
  history_->back().SetTime(0.0, 0.0);
  if (params_->degoption_ == "ExpVar")
    dserror("variable degradation not combinable with adaptive time stepping");
}

/*----------------------------------------------------------------------*
 |  Evaluate                                      (public)         12/10|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::Evaluate(const LINALG::Matrix<3,3>* defgrd,
                                      const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
                                      Teuchos::ParameterList& params,
                                      LINALG::Matrix<NUM_STRESS_3D,1>* stress,
                                      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
                                      const int eleGID)
{
  // get gauss point number
  const int gp = params.get<int>("gp",-1);
  if (gp == -1) dserror("no Gauss point number provided in material");
  // get element lID incase we have element specific material parameters
  int eleID = DRT::Problem::Instance()->GetDis("structure")->ElementColMap()->LID(eleGID);
  const double growthfactor =  params_->GetParameter(params_->growthfactor,eleID);

  // get variables from params
  double dt = params.get<double>("delta time",-1.0);
  double time = params.get<double>("total time",-1.0);
  std::string action = params.get<std::string>("action","none");
  bool output = false;
  if (action == "calc_struct_stress") output = true;
  // int eleid = params.get<int>("eleID",0);
  double inner_radius = 0.0;
  if (params_->sheargrowthfactor_ > 0.0)
    inner_radius = params.get<double>("inner radius");
  double eps = 1.0e-11;

  // elastin degradation
  double elastin_survival = 1.0;
  if (time > params_->starttime_ + eps)
  {
    if (*params_->elastindegrad_ == "Rectangle" || *params_->elastindegrad_ == "RectanglePlate"
        || *params_->elastindegrad_ == "Wedge" || *params_->elastindegrad_ == "Circles")
    {
      LINALG::Matrix<1,3> point_refe(true);
      point_refe = params.get<LINALG::Matrix<1,3> >("gprefecoord");
      ElastinDegradation(point_refe, elastin_survival);
    }
    else if (*params_->elastindegrad_ == "Time")
    {
      elastin_survival = exp(-(time - params_->starttime_) * log(2.0) / 14600.0);
    }
  }
  else if (*params_->initstretch_ == "SetHomeo" || *params_->initstretch_ == "SetLinearHistory"
      || *params_->initstretch_ == "SetConstantHistory")
  {
    if (*params_->elastindegrad_ == "Rectangle" || *params_->elastindegrad_ == "RectanglePlate"
        || *params_->elastindegrad_ == "Wedge" || *params_->elastindegrad_ == "Circles")
    {
      LINALG::Matrix<1,3> point_refe(true);
      point_refe = params.get<LINALG::Matrix<1,3> >("gprefecoord");
      ElastinDegradation(point_refe, elastin_survival);
    }
  }
  if (*params_->elastindegrad_ == "InvEla")
    elastin_survival = params_->GetParameter(params_->elastin_survival,eleID);

  // stuff for collagen damage
  deletemass_->resize(0);

  // differentiate between explicit and implicit description
  int firstiter = 0;
  if (*params_->integration_ == "Explicit")
    firstiter = 1;

  // determine minimal index
  if (params_->storehistory_)
  {
    double acttime = 0.0;
    minindex_ = 0;
    int sizehistory = history_->size();
    double temptime = 0.0;
    double tempdt = 0.0;
    history_->at(minindex_).GetTime(&temptime, &tempdt);
    history_->back().GetTime(&acttime, &tempdt);
    if (acttime == 0.0 || time <= acttime)
      acttime = time;
    double degrad = 0.0;
    Degradation(acttime-temptime, degrad);
    while (degrad < params_->degtol_ && minindex_ < sizehistory - 1)
    {
      minindex_ +=1;
      history_->at(minindex_).GetTime(&temptime, &tempdt);
      Degradation(acttime-temptime, degrad);
    }
  }

  // this is not nice, but I see no different way to make inverse analysis of prestretch work
  // only in first time step
  if (abs(time -dt) < eps && *params_->initstretch_ == "UpdatePrestretch")
  {
    if (params_->numhom_ == 1)
    {
      localprestretch_->at(gp).PutScalar(params_->prestretchcollagen_[0]);
      localhomstress_->at(gp).PutScalar(params_->homstress_[0]);
    }
    else
    {
      localprestretch_->at(gp)(0) = params_->prestretchcollagen_[0];
      localprestretch_->at(gp)(1) = params_->prestretchcollagen_[1];
      localprestretch_->at(gp)(2) = params_->prestretchcollagen_[2];
      localprestretch_->at(gp)(3) = params_->prestretchcollagen_[2];
      localhomstress_->at(gp)(0) = params_->homstress_[0];
      localhomstress_->at(gp)(1) = params_->homstress_[1];
      localhomstress_->at(gp)(2) = params_->homstress_[2];
      localhomstress_->at(gp)(3) = params_->homstress_[2];
    }
  }

  if (!output)
  {
    // set actual time as it might have changed after an restart etc. but just once
    double temptime = 0.0;
    double tempdt = 0.0;
    history_->back().GetTime(&temptime, &tempdt);
    if (temptime == 0.0 && tempdt == 0.0)
    {
      int sizehistory = history_->size();
      history_->at(sizehistory-2).GetTime(&temptime, &tempdt);
      // for restart the function ApplyForceInternal calls the material with the old time
      // (i.e. time = temptime) thus make sure not to store it
      if (time > temptime + 1.0e-11)
      {
        history_->back().SetTime(time, dt);
        temptime = time;
        // if you change your time step size the basal mass production rate changes
        // basal mass production rate determined by DENS, PHIE and degradation function
        double intdegr = 0.0;
        double degrtime = 0.0;
        double degrdt = 0.0;
        for (int idpast = minindex_; idpast < sizehistory - firstiter; idpast++)
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
          Degradation(time-degrtime, degr);
          intdegr += degr * degrdt;
        }
        massprodbasal_ = (1.0 - params_->phimuscle_ - params_->phielastin_) * params_->density_ / 4.0 / intdegr;
//        massprodbasal_ = (1.0 - params_->phimuscle_ - params_->phielastin_) * params_->density_ / 10.0 / intdegr;

        // update degradation values (not possible in Update, as dt not known there)
        if (params_->degoption_ == "ExpVar")
        {
          if (*params_->integration_ == "Implicit")
            dserror("ExpVar not implemented for implicit time integration");
          int sizehistory = history_->size();
          // stretch of previous time step
          LINALG::Matrix<4,1> actstretch(true);
          history_->at(sizehistory-2).GetStretches(gp, &actstretch);
          for (int idfiber = 0; idfiber < 4; idfiber ++)
          {
            double homstrain = 0.0;
            double fac_cmat = 0.0;
            EvaluateSingleFiberScalars(localprestretch_->at(gp)(idfiber)*localprestretch_->at(gp)(idfiber), fac_cmat, homstrain);
            homstrain = homstrain * localprestretch_->at(gp)(idfiber);
            for (int idtime = minindex_; idtime < sizehistory-2; idtime++)
            {
              LINALG::Matrix<4,1> depstretch(true);
              history_->at(idtime).GetStretches(gp, &depstretch);
              double strain = 0.0;
              double fac_cmat = 0.0;
              double stretch = localprestretch_->at(gp)(idfiber) * actstretch(idfiber) / depstretch(idfiber);
              double I4 = stretch * stretch;
              EvaluateSingleFiberScalars(I4, fac_cmat, strain);
              strain = strain * stretch;
              double olddegrad = 0.0;
              history_->at(idtime).GetVarDegrad(gp,idfiber,&olddegrad);
              double newdegrad = exp(-(strain/homstrain - 1.0) * (strain/homstrain - 1.0) * dt * log(2.0) / params_->lifetime_);
              newdegrad = newdegrad * olddegrad;
              history_->at(idtime).SetVarDegrad(gp, idfiber, newdegrad);
            }
          }
        } // expvar
      }
    }
    else if (time > temptime + 1.0e-11)
    {
      // in remodeling time might be wrong depending on the time integration used
      // correct this for the computation but do not store it
      time = temptime;
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
    else if (*params_->initstretch_ == "SetConstantHistory" && time > (0.9*params_->starttime_ - dt + 1.0e-12) && time <= (0.9*params_->starttime_ + 1.0e-12))
    {
      LINALG::Matrix<4,1> tempstretch(actstretch);
      for (int i=0; i<4; i++)
      {
        if (tempstretch(i) < 1.0)
          tempstretch(i) = 1.0;
      }
      history_->back().SetStretches(gp, tempstretch);
    }
    else if (*params_->initstretch_ == "SetLinearHistory" && time <= (0.9*params_->starttime_ + 1.0e-12) && time > (0.9*params_->starttime_ - dt + 1.0e-12))
    {
      LINALG::Matrix<4,1> tempstretch(actstretch);
      history_->back().SetStretches(gp, tempstretch);
    }

    // set prestretch according to time curve or adapt prestretch
    if (time < params_->starttime_ - eps && params_->timecurve_ != 0)
    {
      // increase prestretch according to time curve
      int curvenum = params_->timecurve_;
      double curvefac = 1.0;
      // numbering starts from zero here, thus use curvenum-1
      if (curvenum)
        curvefac = DRT::Problem::Instance()->Curve(curvenum-1).f(time);
      if (curvefac > (1.0 + eps) || curvefac < (0.0 - eps))
        dserror("correct your time curve for prestretch, just values in [0,1] are allowed %f", curvefac);
      if (params_->numhom_ == 1)
      {
        double prestretch = 1.0 + (params_->prestretchcollagen_[0] - 1.0)*curvefac;
        localprestretch_->at(gp).PutScalar(prestretch);
      }
      else
      {
        double prestretch = 1.0 + (params_->prestretchcollagen_[0] - 1.0)*curvefac;
        localprestretch_->at(gp)(0) = prestretch;
        prestretch = 1.0 + (params_->prestretchcollagen_[1] - 1.0)*curvefac;
        localprestretch_->at(gp)(1) = prestretch;
        prestretch = 1.0 + (params_->prestretchcollagen_[2] - 1.0)*curvefac;
        localprestretch_->at(gp)(2) = prestretch;
        localprestretch_->at(gp)(3) = prestretch;
      }
    }
    else if (abs(time - params_->starttime_) < eps && *params_->initstretch_ == "UpdatePrestretch")
    {
      // use current stretch as prestretch
      if (params_->numhom_ == 1)
        localprestretch_->at(gp).Update(params_->prestretchcollagen_[0],actstretch);
      else
      {
        localprestretch_->at(gp)(0) = params_->prestretchcollagen_[0] * actstretch(0);
        localprestretch_->at(gp)(1) = params_->prestretchcollagen_[1] * actstretch(1);
        localprestretch_->at(gp)(2) = params_->prestretchcollagen_[2] * actstretch(2);
        localprestretch_->at(gp)(3) = params_->prestretchcollagen_[2] * actstretch(3);
      }
      // adopt deposition stretch
      int numsteps = history_->size();
      for (int i = 0; i < numsteps; i++)
        history_->at(i).SetStretches(gp, actstretch);
      homradius_ = inner_radius;
    }

    // start in every iteration from the original value, this is important for implicit only
    LINALG::Matrix<4,1> massprodstart(true);
    for (int id = 0; id < 4; id++)
      massprodstart(id) = massprodbasal_;
//    massprodstart(2) = massprodstart(2)*4;
//    massprodstart(3) = massprodstart(3)*4;
    history_->back().SetMass(gp, massprodstart);

    EvaluateStress(glstrain, gp, cmat, stress, firstiter, time, elastin_survival);

    //--------------------------------------------------------------------------------------
    // compute new deposition rates
    // either for future use or just for visualization
    LINALG::Matrix<4,1> massstress(true);
    LINALG::Matrix<4,1> massprodcomp(true);
    if (*params_->growthforce_ == "All")
    {
      MassProduction(gp, *defgrd, *stress, &massstress, inner_radius, &massprodcomp, growthfactor);
    }
    else if (*params_->growthforce_ == "ElaCol")
    {
      double masstemp = 0.0;
      LINALG::Matrix<NUM_STRESS_3D,1> stresstemp(true);
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmattemp(true);

      // 1st step: elastin
      //==========================
      LINALG::Matrix<NUM_STRESS_3D,1> Siso(true);
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatiso(true);
      EvaluateElastin(C, &cmatiso, &Siso, time, &masstemp, elastin_survival);
      stresstemp = Siso;
      cmattemp = cmatiso;

      // 2nd step: collagen
      //==========================
      EvaluateFiberFamily(C, gp, &cmattemp, &stresstemp, a1_->at(gp), &masstemp, firstiter, time, 0);
      EvaluateFiberFamily(C, gp, &cmattemp, &stresstemp, a2_->at(gp), &masstemp, firstiter, time, 1);
      EvaluateFiberFamily(C, gp, &cmattemp, &stresstemp, a3_->at(gp), &masstemp, firstiter, time, 2);
      EvaluateFiberFamily(C, gp, &cmattemp, &stresstemp, a4_->at(gp), &masstemp, firstiter, time, 3);

      MassProduction(gp, *defgrd, stresstemp, &massstress, inner_radius, &massprodcomp, growthfactor);
    }
    else
    {
      double massstresstemp = 0.0;
      double massprodtemp = 0.0;
      LINALG::Matrix<NUM_STRESS_3D,1> stresstemp(true);
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmattemp(true);
      double masstemp = 0.0;
      EvaluateFiberFamily(C, gp, &cmattemp, &stresstemp, a1_->at(gp), &masstemp, firstiter, time, 0);
      MassProductionSingleFiber(gp, *defgrd, stresstemp, &massstresstemp, inner_radius, &massprodtemp, a1_->at(gp), 0, growthfactor);
      massstress(0) = massstresstemp;
      massprodcomp(0) = massprodtemp;
      stresstemp.Scale(0.0);
      EvaluateFiberFamily(C, gp, &cmattemp, &stresstemp, a2_->at(gp), &masstemp, firstiter, time, 1);
      MassProductionSingleFiber(gp, *defgrd, stresstemp, &massstresstemp, inner_radius, &massprodtemp, a2_->at(gp), 1, growthfactor);
      massstress(1) = massstresstemp;
      massprodcomp(1) = massprodtemp;
      stresstemp.Scale(0.0);
      EvaluateFiberFamily(C, gp, &cmattemp, &stresstemp, a3_->at(gp), &masstemp, firstiter, time, 2);
      MassProductionSingleFiber(gp, *defgrd, stresstemp, &massstresstemp, inner_radius, &massprodtemp, a3_->at(gp), 2, growthfactor);
      massstress(2) = massstresstemp;
      massprodcomp(2) = massprodtemp;
      stresstemp.Scale(0.0);
      EvaluateFiberFamily(C, gp, &cmattemp, &stresstemp, a4_->at(gp), &masstemp, firstiter, time, 3);
      MassProductionSingleFiber(gp, *defgrd, stresstemp, &massstresstemp, inner_radius, &massprodtemp, a4_->at(gp), 3, growthfactor);
      massstress(3) = massstresstemp;
      massprodcomp(3) = massprodtemp;
    }

    // set new homstress
    if (abs(time - params_->starttime_) < eps && *params_->initstretch_ == "UpdatePrestretch")
    {
      localhomstress_->at(gp).Update(massstress);
      // perhaps add time < params_->starttime_? correction factor actstretch needed?
    }

    if (time > params_->starttime_ + eps && (growthfactor != 0.0 || params_->sheargrowthfactor_ != 0.0))
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
      else
      {
        if (*params_->massprodfunc_ != "Lin" || *params_->initstretch_ == "SetConstantHistory"
            || *params_->initstretch_ == "SetLinearHistory")
          dserror("Your desired option of elastin degradation, mass production function or initstretch\n is not implemented in implicit time integration");
        if (*params_->growthforce_ == "All")
        {
          EvaluateImplicitAll(*defgrd, glstrain, gp, cmat, stress, dt, time, massprodcomp, massstress, elastin_survival, growthfactor);
        }
        else if (*params_->growthforce_ == "Single")
        {
          EvaluateImplicitSingle(*defgrd, glstrain, gp, cmat, stress, dt, time, elastin_survival, growthfactor);
        }
        else if (*params_->growthforce_ == "ElaCol")
        {
          dserror("GROWTHFORCE ElaCol not implemented for implicit integration");
        }
      }

    } else {
      // visualization of massstresss for the other cases
      vismassstress_->at(gp)(0) = massstress(0);
      vismassstress_->at(gp)(1) = massstress(1);
      vismassstress_->at(gp)(2) = massstress(2);
      if ((*params_->initstretch_ == "SetConstantHistory" || *params_->initstretch_ == "SetLinearHistory")
          && time > (0.6*params_->starttime_ + 1.0e-12) && time <= (0.9*params_->starttime_ - dt + 1.0e-12))
        history_->back().SetMass(gp,massprodcomp);
    }
  }
  else
  {
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
        EvaluateStress(glstrain, gp, cmat, stress, firstiter, temptime, elastin_survival);
        dserror("has to be checked, update called before output");
      } else
        dserror("times do not match: %f actual time, %f deposition time of last fiber",time,temptime);
    }
    else
    {
      EvaluateStress(glstrain, gp, cmat, stress, firstiter, time, elastin_survival);
    }
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
  double time,
  double elastin_survival
)
{
  //--------------------------------------------------------------------------------------
  // some variables
  double density = params_->density_;
  double currmassdens = 0.0;

  //--------------------------------------------------------------------------------------
  // build identity tensor I
  LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
  C.Scale(2.0);
  C += Id;

  //--------------------------------------------------------------------------------------
  // calculate stress and elasticity matrix

  // 1st step: elastin
  //==========================
  LINALG::Matrix<NUM_STRESS_3D,1> Siso(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatiso(true);
  EvaluateElastin(C, &cmatiso, &Siso, time, &currmassdens, elastin_survival);
  (*stress) = Siso;
  (*cmat) = cmatiso;

  // 2nd step: collagen
  //==========================
  EvaluateFiberFamily(C, gp, cmat, stress, a1_->at(gp), &currmassdens, firstiter, time, 0);

  EvaluateFiberFamily(C, gp, cmat, stress, a2_->at(gp), &currmassdens, firstiter, time, 1);

  EvaluateFiberFamily(C, gp, cmat, stress, a3_->at(gp), &currmassdens, firstiter, time, 2);

  EvaluateFiberFamily(C, gp, cmat, stress, a4_->at(gp), &currmassdens, firstiter, time, 3);

  // 3rd step: smooth muscle
  //==========================
  LINALG::Matrix<NUM_STRESS_3D,1> Smus(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatmus(true);
  EvaluateMuscle(C, &cmatmus, &Smus, gp, &currmassdens);
  (*stress) += Smus;
  (*cmat) += cmatmus;

  // 4th step: volumetric part
  //==========================
  LINALG::Matrix<NUM_STRESS_3D,1> Svol(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatvol(true);
  EvaluateVolumetric(C, &cmatvol, &Svol, currmassdens, density);
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
  const LINALG::Matrix<NUM_STRESS_3D,1> C,
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
  double prestretchcollagen = localprestretch_->at(gp)(idfiber);
  double density = params_->density_;
  int sizehistory = history_->size();
  double eps = 1.0e-11;
  if (idfiber != 3)
    visrefmassdens_->at(gp)(idfiber) = 0.0;

  //--------------------------------------------------------------------------------------
  // structural tensors in voigt notation
  // A = a x a
  LINALG::Matrix<NUM_STRESS_3D,1>  A;
  for (int i = 0; i < 3; i++)
    A(i) = a(i)*a(i);

  A(3) = a(0)*a(1); A(4) = a(1)*a(2); A(5) = a(0)*a(2);

  double I4 =  A(0)*C(0) + A(1)*C(1) + A(2)*C(2)
             + 1.*(A(3)*C(3) + A(4)*C(4) + A(5)*C(5)); // I4 = trace(A C)

  // variables for stress and cmat
  double fac_cmat = 0.0;
  double fac_stress = 0.0;

  //--------------------------------------------------------------------------------------
  // prestress time
//  const Teuchos::ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
//  INPAR::STR::PreStress pstype = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(pslist,"PRESTRESS");
//  double pstime = -1.0 * params_->lifetime_ - dt;
//  if (pstype == INPAR::STR::prestress_mulf)
//    pstime = pslist.get<double>("PRESTRESSTIME");

  //--------------------------------------------------------------------------------------
  // calculate stress and elasticity matrix
  for (int idpast = minindex_; idpast < sizehistory - firstiter; idpast++)
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
    }
    else if (deptime >= (time - params_->lifetime_ - eps) &&
             (deptime - depdt) < (time - params_->lifetime_ - eps) &&
             (params_->degoption_ != "Exp"  && params_->degoption_ != "ExpVar"))
    {
      depdt = deptime - time + params_->lifetime_;
    }

    LINALG::Matrix<4,1> collstretch(true);
    history_->at(idpast).GetStretches(gp, &collstretch);
    double stretch = prestretchcollagen / collstretch(idfiber);
    // prestretch of collagen fibers is not applied, might be reasonable combined with prestress
    if (*params_->initstretch_ == "experimental" && deptime <= params_->starttime_ + eps)
      stretch = 1.0 / collstretch(idfiber);

    double I4_loc = I4;

    // linear distribution of stretch
    if (*params_->initstretch_ == "SetLinearHistory" && time <= (0.9*params_->starttime_ + 1.0e-12))
    {
      I4_loc = ((sqrt(I4)-1.0)*(1.0-idpast/(sizehistory-1.0))+1.0)*((sqrt(I4)-1.0)*(1.0-idpast/(sizehistory-1.0))+1.0);
      LINALG::Matrix<4,1> tempstretch(true);
      history_->at(idpast).GetStretches(gp,&tempstretch);
      if (abs(tempstretch(idfiber)-1.0)>1.0e-12)
        dserror("linear stretch when stretch history has been modified");
    }

    I4_loc = I4_loc* stretch*stretch;  // account for prestretch and stretch at deposition time
    if (sqrt(I4_loc) > params_->damagestretch_)
    {
      LINALG::Matrix<3,1> deldata(true);
      deldata(0) = gp;
      deldata(1) = idpast;
      deldata(2) = idfiber;
      deletemass_->push_back(deldata);
    }

    double fac_cmat_loc = 0.0;
    double fac_stress_loc = 0.0;
    EvaluateSingleFiberScalars(I4_loc, fac_cmat_loc, fac_stress_loc);

    if (*params_->initstretch_ == "SetLinearHistory" && time <= (0.9*params_->starttime_ + 1.0e-12))
      fac_cmat_loc = fac_cmat_loc * ((sqrt(I4)-1.0)*(1.0-idpast/(sizehistory-1.0))+1.0)*(1.0-idpast/(sizehistory-1.0))/sqrt(I4);

    double qdegrad = 0.0;
    Degradation(time - deptime, qdegrad);
    if (params_->degoption_ == "ExpVar")
    {
      double vardegrad;
      history_->at(idpast).GetVarDegrad(gp,idfiber,&vardegrad);
      qdegrad = qdegrad * vardegrad;
    }
    LINALG::Matrix<4,1> collmass(true);
    history_->at(idpast).GetMass(gp, &collmass);
    fac_stress += fac_stress_loc * stretch*stretch * qdegrad * collmass(idfiber) / density * depdt;

    if (*params_->initstretch_ == "Homeo" || (idpast == sizehistory - 1 && time > params_->starttime_ + 1.0e-11))
    {
      LINALG::Matrix<NUM_STRESS_3D,1> Saniso_loc(A);
      Saniso_loc.Scale(fac_stress_loc * stretch*stretch * qdegrad * collmass(idfiber) / density * depdt);
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatanisoadd(true);
      cmatanisoadd.MultiplyNT(A, Saniso_loc);
      cmatanisoadd.Scale(-2.0/(collstretch(idfiber)*collstretch(idfiber)));
      (*cmat) += cmatanisoadd;
    }
    else
    {
      fac_cmat += fac_cmat_loc * stretch*stretch*stretch*stretch * qdegrad * collmass(idfiber) / density * depdt;
    }

    (*currmassdens) += qdegrad * collmass(idfiber) * depdt;
    if (idfiber != 3)
      visrefmassdens_->at(gp)(idfiber) += qdegrad * collmass(idfiber) * depdt;
  }

  // matrices for stress and cmat
  LINALG::Matrix<NUM_STRESS_3D,1> Saniso(A);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmataniso(true);
  cmataniso.MultiplyNT(A,A);
  Saniso.Scale(fac_stress);
  cmataniso.Scale(fac_cmat);
  (*stress) += Saniso;
  (*cmat) += cmataniso;
}

/*----------------------------------------------------------------------*
 |  EvaluateSingleFiberScalars                    (private)        02/12|
 *----------------------------------------------------------------------*
 strain energy function

 Psi    = k1/(2.0*k2)*(exp(k2*(I_4 - 1.0)^2)-1.0)

 I_4 .. invariant accounting for the fiber direction
 */
void MAT::ConstraintMixture::EvaluateSingleFiberScalars
(
  double I4,
  double& fac_cmat,
  double& fac_stress
)
{
  const double k1 = params_->k1_;
  const double k2 = params_->k2_;

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

  // stress
  const double exp1 = std::exp(k2*(I4-1.)*(I4-1.));
  if (std::isinf(exp1))
    dserror("stretch in fiber direction is too high %e", sqrt(I4));
  fac_stress = 2.*(k1*(I4-1.)*exp1);  // 2 dW/dI4

  // cmat
  fac_cmat = fib1_tension* 4.0*(k1*exp1 + 2.0*k1*k2*(I4-1.0)*(I4-1.0)*exp1); // 4 d^2Wf/dI4dI4
}

/*----------------------------------------------------------------------*
 |  EvaluateElastin                               (private)        12/10|
 *----------------------------------------------------------------------*
 strain energy function

 W    = 1/2 mue (I1-3) + mue / (2 beta) (I3^-beta -1)

*/
void MAT::ConstraintMixture::EvaluateElastin
(
  const LINALG::Matrix<NUM_STRESS_3D,1> C,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress,
  double time,
  double* currmassdens,
  double elastin_survival
)
{
  double density = params_->density_;
  double prestretchelastin = params_->prestretchelastin_;
  double refmassdenselastin = params_->phielastin_ * density;
  if (refmassdenselastin < 0.0 || refmassdenselastin > density) dserror("mass fraction of elastin not in [0;1]");

  if (time < params_->starttime_ - 1.0e-11 && params_->timecurve_ != 0)
  {
    // increase prestretch according to time curve
    int curvenum = params_->timecurve_;
    double curvefac = 1.0;
    // numbering starts from zero here, thus use curvenum-1
    if (curvenum)
      curvefac = DRT::Problem::Instance()->Curve(curvenum-1).f(time);
    if (curvefac > 1.0 || curvefac < 0.0)
      dserror("correct your time curve for prestretch, just values in [0,1] are allowed %f", curvefac);
    prestretchelastin = 1.0 + (params_->prestretchelastin_ - 1.0)*curvefac;
  }
  // account for isotropic prestretch of elastin
  LINALG::Matrix<NUM_STRESS_3D,1> Ciso(C);
  Ciso.Scale(prestretchelastin*prestretchelastin);

  double mue = params_->mue_;
  double nue = params_->nue_;
  double beta = nue / (1.0-2.0*nue);

  // isotropic invariant
  const double I3 = Ciso(0)*Ciso(1)*Ciso(2)
        + 0.25 * Ciso(3)*Ciso(4)*Ciso(5)
        - 0.25 * Ciso(1)*Ciso(5)*Ciso(5)
        - 0.25 * Ciso(2)*Ciso(3)*Ciso(3)
        - 0.25 * Ciso(0)*Ciso(4)*Ciso(4);    // 3rd invariant, determinant

  //--------------------------------------------------------------------------------------
  // invert Ciso
  LINALG::Matrix<NUM_STRESS_3D,1> Cisoinv(true);
  Cisoinv(0) = Ciso(1)*Ciso(2) - 0.25*Ciso(4)*Ciso(4);
  Cisoinv(1) = Ciso(0)*Ciso(2) - 0.25*Ciso(5)*Ciso(5);
  Cisoinv(2) = Ciso(0)*Ciso(1) - 0.25*Ciso(3)*Ciso(3);
  Cisoinv(3) = 0.25*Ciso(5)*Ciso(4) - 0.5*Ciso(3)*Ciso(2);
  Cisoinv(4) = 0.25*Ciso(3)*Ciso(5) - 0.5*Ciso(0)*Ciso(4);
  Cisoinv(5) = 0.25*Ciso(3)*Ciso(4) - 0.5*Ciso(5)*Ciso(1);
  Cisoinv.Scale(1.0/I3);

  LINALG::Matrix<NUM_STRESS_3D,1> Siso(true);
  for (int i = 0; i < 3; i++) Siso(i) = mue;

  double gamma2 = - mue*pow(I3,-beta);
  Siso.Update(gamma2, Cisoinv, 1.0);
  Siso.Scale(refmassdenselastin/density*prestretchelastin*prestretchelastin*elastin_survival);
  *stress = Siso;

  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatiso(true);
  double delta6 = 2.*beta*mue*pow(I3, -beta);
  cmatiso.MultiplyNT(delta6, Cisoinv, Cisoinv);
  double delta7 = 2.*mue*pow(I3, -beta);
  AddtoCmatHolzapfelProduct(cmatiso, Cisoinv, delta7);
  cmatiso.Scale(refmassdenselastin/density*prestretchelastin*prestretchelastin*prestretchelastin*prestretchelastin*elastin_survival);
  *cmat = cmatiso;

  (*currmassdens) += refmassdenselastin;
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
  const LINALG::Matrix<NUM_STRESS_3D,1> C,
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
  // structural tensors in voigt notation
  // A = a x a
  LINALG::Matrix<3,1>  a = a1_->at(gp);
  LINALG::Matrix<NUM_STRESS_3D,1>  A;
  for (int i = 0; i < 3; i++)
    A(i) = a(i)*a(i);

  A(3) = a(0)*a(1); A(4) = a(1)*a(2); A(5) = a(0)*a(2);

  double I4 =  A(0)*C(0) + A(1)*C(1) + A(2)*C(2)
             + 1.*(A(3)*C(3) + A(4)*C(4) + A(5)*C(5)); // I4 = trace(A C)

  if (massfrac > 0.0 + 1.0e-12 && k1 > 0.0 + 1.0e-12)
  {
    //++++++++++++++++++++++ passive part ++++++++++++++++++
    //--- determine 2nd Piola Kirchhoff stresses S -----------------------------------------
    double preI4 = I4 * prestretchmuscle*prestretchmuscle;
    LINALG::Matrix<NUM_STRESS_3D,1> Saniso(A); // first compute S = 2 dW/dI4 A
    const double exp1 = std::exp(k2*(preI4-1.)*(preI4-1.));
    if (std::isinf(exp1))
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
  }

  //++++++++++++++++++++++ active part ++++++++++++++++++
  // does not depend on mass fraction
  double lambda = sqrt(I4);
  double Smax = params_->Smax_; //50000;
  double lambda_M = 1.2; //1.2,1.4;
  double lambda_0 = 0.7; //0.7,0.8;
  // perhaps check if lambda is between lambda_M and lambda_0

  LINALG::Matrix<NUM_STRESS_3D,1> Saniso(A);
  double facS = Smax / lambda * (1.0 - ((lambda_M-lambda)*(lambda_M-lambda)/(lambda_M-lambda_0)/(lambda_M-lambda_0)));
  Saniso.Scale(facS);
  *stress += Saniso;

  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmataniso(true);
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
  const LINALG::Matrix<NUM_STRESS_3D,1> C,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress,
  double currMassDens,
  double refMassDens
)
{
  double kappa = params_->kappa_;

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
  LINALG::Matrix<3,3> defgrd,
  LINALG::Matrix<NUM_STRESS_3D,1> S,
  LINALG::Matrix<4,1>* massstress,
  double inner_radius,
  LINALG::Matrix<4,1>* massprodcomp,
  double growthfactor
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
  // one has to use the defgrd here, as with EAS C != F^T F
  LINALG::Matrix<3,3> Cmatrix(true);
  Cmatrix.MultiplyTN(defgrd,defgrd);
  double detC = Cmatrix.Determinant();

  LINALG::Matrix<3,3> temp1(true);
  LINALG::Matrix<3,3> temp2(true);
  temp1.Multiply(Smatrix,Cmatrix);
  temp2.Multiply(Cmatrix,temp1);
  temp1.Scale(0.0);
  temp1.Multiply(Smatrix,temp2);
  temp2.Scale(0.0);
  temp2.Multiply(1.0/detC,Cmatrix,temp1);

  double sheardiff = 0.0;
  double sheargrowthfactor = params_->sheargrowthfactor_;
  if (sheargrowthfactor > 0.0)
    sheardiff = 1.0 - homradius_*homradius_*homradius_ / (inner_radius*inner_radius*inner_radius);
  double maxmassprodfac = params_->maxmassprodfac_;

  // Fiber1
  LINALG::Matrix<3,1> temp_vector(true);
  LINALG::Matrix<1,1> temp_scalar(true);
  temp_vector.Multiply(Cmatrix,a1_->at(gp));
  temp_scalar.MultiplyTN(a1_->at(gp),temp_vector);
  double currentstretch = sqrt(temp_scalar(0));
  temp_vector.Multiply(temp2,a1_->at(gp));
  temp_scalar.MultiplyTN(a1_->at(gp),temp_vector);
  double massstress1 = sqrt(temp_scalar(0)) / currentstretch;
  if (*params_->massprodfunc_ == "Lin")
  {
    if (localhomstress_->at(gp)(0) != 0.0)
      (*massprodcomp)(0) = massprodbasal_ * (1.0 + growthfactor * (massstress1 / localhomstress_->at(gp)(0) - 1.0) + sheargrowthfactor * sheardiff);
    else
      (*massprodcomp)(0) = massprodbasal_ * (1.0 + growthfactor * massstress1 + sheargrowthfactor * sheardiff);
    if ((*massprodcomp)(0) > (maxmassprodfac * massprodbasal_))
      (*massprodcomp)(0) = maxmassprodfac * massprodbasal_;
  }
  else if (*params_->massprodfunc_ == "CosCos")
  {
    if (localhomstress_->at(gp)(0) != 0.0)
    {
      double facstress = 0.0;
      double deltastress = massstress1 / localhomstress_->at(gp)(0) - 1.0;
      MassFunction(growthfactor, deltastress, maxmassprodfac, facstress);
      double facshear = 0.0;
      MassFunction(sheargrowthfactor, sheardiff, maxmassprodfac, facshear);
      (*massprodcomp)(0) = massprodbasal_ * 0.5 * (facstress + facshear);
    }
    else
    {
      double facstress = 0.0;
      double deltastress = massstress1;
      MassFunction(growthfactor, deltastress, maxmassprodfac, facstress);
      double facshear = 0.0;
      MassFunction(sheargrowthfactor, sheardiff, maxmassprodfac, facshear);
      (*massprodcomp)(0) = massprodbasal_ * 0.5 * (facstress + facshear);
    }
  }
  if ((*massprodcomp)(0) < 0.0)
  {
    (*massprodcomp)(0) = 0.0;
    //std::cout << "warning negative massproduction rate" << std::endl;
  }

  // Fiber2
  temp_vector.Multiply(Cmatrix,a2_->at(gp));
  temp_scalar.MultiplyTN(a2_->at(gp),temp_vector);
  currentstretch = sqrt(temp_scalar(0));
  temp_vector.Multiply(temp2,a2_->at(gp));
  temp_scalar.MultiplyTN(a2_->at(gp),temp_vector);
  double massstress2 = sqrt(temp_scalar(0)) / currentstretch;
  if (*params_->massprodfunc_ == "Lin")
  {
    if (localhomstress_->at(gp)(1) != 0.0)
      (*massprodcomp)(1) = massprodbasal_ * (1.0 + growthfactor * (massstress2 / localhomstress_->at(gp)(1) - 1.0) + sheargrowthfactor * sheardiff);
    else
      (*massprodcomp)(1) = massprodbasal_ * (1.0 + growthfactor * massstress2 + sheargrowthfactor * sheardiff);
    if ((*massprodcomp)(1) > (maxmassprodfac * massprodbasal_))
      (*massprodcomp)(1) = maxmassprodfac * massprodbasal_;
  }
  else if (*params_->massprodfunc_ == "CosCos")
  {
    if (localhomstress_->at(gp)(1) != 0.0)
    {
      double facstress = 0.0;
      double deltastress = massstress2 / localhomstress_->at(gp)(1) - 1.0;
      MassFunction(growthfactor, deltastress, maxmassprodfac, facstress);
      double facshear = 0.0;
      MassFunction(sheargrowthfactor, sheardiff, maxmassprodfac, facshear);
      (*massprodcomp)(1) = massprodbasal_ * 0.5 * (facstress + facshear);
    }
    else
    {
      double facstress = 0.0;
      double deltastress = massstress2;
      MassFunction(growthfactor, deltastress, maxmassprodfac, facstress);
      double facshear = 0.0;
      MassFunction(sheargrowthfactor, sheardiff, maxmassprodfac, facshear);
      (*massprodcomp)(1) = massprodbasal_ * 0.5 * (facstress + facshear);
    }
  }
  if ((*massprodcomp)(1) < 0.0)
  {
    (*massprodcomp)(1) = 0.0;
    //std::cout << "warning negative massproduction rate" << std::endl;
  }

  // Fiber3
  temp_vector.Multiply(Cmatrix,a3_->at(gp));
  temp_scalar.MultiplyTN(a3_->at(gp),temp_vector);
  currentstretch = sqrt(temp_scalar(0));
  temp_vector.Multiply(temp2,a3_->at(gp));
  temp_scalar.MultiplyTN(a3_->at(gp),temp_vector);
  double massstress3 = sqrt(temp_scalar(0)) / currentstretch;
  if (*params_->massprodfunc_ == "Lin")
  {
    if (localhomstress_->at(gp)(2) != 0.0)
      (*massprodcomp)(2) = massprodbasal_ * (1.0 + growthfactor * (massstress3 / localhomstress_->at(gp)(2) - 1.0) + sheargrowthfactor * sheardiff); // /4.0 massstress3
    else
      (*massprodcomp)(2) = massprodbasal_ * (1.0 + growthfactor * massstress3 + sheargrowthfactor * sheardiff);
    if ((*massprodcomp)(2) > (maxmassprodfac * massprodbasal_))
      (*massprodcomp)(2) = maxmassprodfac * massprodbasal_;
  }
  else if (*params_->massprodfunc_ == "CosCos")
  {
    if (localhomstress_->at(gp)(2) != 0.0)
    {
      double facstress = 0.0;
      double deltastress = massstress3 / localhomstress_->at(gp)(2) - 1.0;
      MassFunction(growthfactor, deltastress, maxmassprodfac, facstress);
      double facshear = 0.0;
      MassFunction(sheargrowthfactor, sheardiff, maxmassprodfac, facshear);
      (*massprodcomp)(2) = massprodbasal_ * 0.5 * (facstress + facshear);
    }
    else
    {
      double facstress = 0.0;
      double deltastress = massstress3;
      MassFunction(growthfactor, deltastress, maxmassprodfac, facstress);
      double facshear = 0.0;
      MassFunction(sheargrowthfactor, sheardiff, maxmassprodfac, facshear);
      (*massprodcomp)(2) = massprodbasal_ * 0.5 * (facstress + facshear);
    }
  }
//  (*massprodcomp)(2) = 4. * (*massprodcomp)(2);
  if ((*massprodcomp)(2) < 0.0)
  {
    (*massprodcomp)(2) = 0.0;
    //std::cout << "warning negative massproduction rate" << std::endl;
  }

  // Fiber4
  temp_vector.Multiply(Cmatrix,a4_->at(gp));
  temp_scalar.MultiplyTN(a4_->at(gp),temp_vector);
  currentstretch = sqrt(temp_scalar(0));
  temp_vector.Multiply(temp2,a4_->at(gp));
  temp_scalar.MultiplyTN(a4_->at(gp),temp_vector);
  double massstress4 = sqrt(temp_scalar(0)) / currentstretch;
  if (*params_->massprodfunc_ == "Lin")
  {
    if (localhomstress_->at(gp)(3) != 0.0)
      (*massprodcomp)(3) = massprodbasal_ * (1.0 + growthfactor * (massstress4 / localhomstress_->at(gp)(3) - 1.0) + sheargrowthfactor * sheardiff); // /4.0 massstress4
    else
      (*massprodcomp)(3) = massprodbasal_ * (1.0 + growthfactor * massstress4 + sheargrowthfactor * sheardiff);
    if ((*massprodcomp)(3) > (maxmassprodfac * massprodbasal_))
      (*massprodcomp)(3) = maxmassprodfac * massprodbasal_;
  }
  else if (*params_->massprodfunc_ == "CosCos")
  {
    if (localhomstress_->at(gp)(3) != 0.0)
    {
      double facstress = 0.0;
      double deltastress = massstress4 / localhomstress_->at(gp)(3) - 1.0;
      MassFunction(growthfactor, deltastress, maxmassprodfac, facstress);
      double facshear = 0.0;
      MassFunction(sheargrowthfactor, sheardiff, maxmassprodfac, facshear);
      (*massprodcomp)(3) = massprodbasal_ * 0.5 * (facstress + facshear);
    }
    else
    {
      double facstress = 0.0;
      double deltastress = massstress4;
      MassFunction(growthfactor, deltastress, maxmassprodfac, facstress);
      double facshear = 0.0;
      MassFunction(sheargrowthfactor, sheardiff, maxmassprodfac, facshear);
      (*massprodcomp)(3) = massprodbasal_ * 0.5 * (facstress + facshear);
    }
  }
//  (*massprodcomp)(3) = 4. * (*massprodcomp)(3);
  if ((*massprodcomp)(3) < 0.0)
  {
    (*massprodcomp)(3) = 0.0;
    //std::cout << "warning negative massproduction rate" << std::endl;
  }

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
  LINALG::Matrix<3,3> defgrd,
  LINALG::Matrix<NUM_STRESS_3D,1> S,
  double* massstress,
  double inner_radius,
  double* massprodcomp,
  LINALG::Matrix<3,1> a,
  const int idfiber,
  double growthfactor
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
  // one has to use the defgrd here, as with EAS C != F^T F
  LINALG::Matrix<3,3> Cmatrix(true);
  Cmatrix.MultiplyTN(defgrd,defgrd);
  double detC = Cmatrix.Determinant();

  double sheardiff = 0.0;
  double sheargrowthfactor = params_->sheargrowthfactor_;
  if (sheargrowthfactor > 0.0)
    sheardiff = 1.0 - homradius_*homradius_*homradius_ / (inner_radius*inner_radius*inner_radius);
  double maxmassprodfac = params_->maxmassprodfac_;

  LINALG::Matrix<3,1> CSCa(true);
  LINALG::Matrix<3,1> SCa(true);
  LINALG::Matrix<1,1> temp(true);
  CSCa.Multiply(Cmatrix,a);
  SCa.Multiply(Smatrix,CSCa);
  CSCa.Multiply(Cmatrix,SCa);
  temp.MultiplyTN(1.0/detC,SCa,CSCa);
  (*massstress) = temp(0);
  LINALG::Matrix<3,1> temp_vector(true);
  LINALG::Matrix<1,1> temp_scalar(true);
  temp_vector.Multiply(Cmatrix,a);
  temp_scalar.MultiplyTN(a,temp_vector);
  double currentstretch = sqrt(temp_scalar(0));
  (*massstress) = sqrt(*massstress) / currentstretch;
  double homstress = localhomstress_->at(gp)(idfiber);
//  if (idfiber == 2 || idfiber == 3)
//    homstress = 4. * homstress;
  if (*params_->massprodfunc_ == "Lin")
  {
    if (homstress != 0.0)
      (*massprodcomp) = massprodbasal_ * (1.0 + growthfactor * ((*massstress) / homstress - 1.0) + sheargrowthfactor * sheardiff);
    else
      (*massprodcomp) = massprodbasal_ * (1.0 + growthfactor * (*massstress) + sheargrowthfactor * sheardiff);
    if (*massprodcomp > (maxmassprodfac * massprodbasal_))
      *massprodcomp = maxmassprodfac * massprodbasal_;
  }
  else if (*params_->massprodfunc_ == "CosCos")
  {
    if (homstress != 0.0)
    {
      double facstress = 0.0;
      double deltastress = (*massstress) / homstress - 1.0;
      MassFunction(growthfactor, deltastress, maxmassprodfac, facstress);
      double facshear = 0.0;
      MassFunction(sheargrowthfactor, sheardiff, maxmassprodfac, facshear);
      (*massprodcomp) = massprodbasal_ * 0.5 * (facstress + facshear);
    }
    else
    {
      double facstress = 0.0;
      double deltastress = (*massstress);
      MassFunction(growthfactor, deltastress, maxmassprodfac, facstress);
      double facshear = 0.0;
      MassFunction(sheargrowthfactor, sheardiff, maxmassprodfac, facshear);
      (*massprodcomp) = massprodbasal_ * 0.5 * (facstress + facshear);
    }
  }
//  if (idfiber == 2 || idfiber == 3)
//    (*massprodcomp) = 4. * (*massprodcomp);
  if ((*massprodcomp) < 0.0)
  {
    (*massprodcomp) = 0.0;
    //std::cout << "warning negative massproduction rate" << std::endl;
  }
}

/*----------------------------------------------------------------------*
 |  MassFunction                                  (private)        06/13|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::MassFunction(double growthfac, double delta, double mmax, double& massfac)
{
  if (delta < - PI / (2.0 * growthfac))
  {
    massfac = 0.0;
  }
  else if (delta < 0.0) // && - PI / (2.0 * growthfac) <= delta
  {
    massfac = 0.5 * ( 1.0 + cos(2.0 * growthfac * delta));
  }
  else if (delta < (mmax -1.0) * PI / (2.0*growthfac)) // && 0.0 < delta
  {
    massfac = 0.5 * (mmax + 1.0 - (mmax - 1.0) * cos(2.0 * growthfac * delta / (mmax - 1.0)));
  }
  else  // (mmax -1.0) * PI / (2.0*growthfac) <= delta
  {
    massfac = mmax;
  }

}

/*----------------------------------------------------------------------*
 |  Degradation                                   (private)        10/11|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::Degradation(double t, double& degr)
{
  if (params_->degoption_ == "Lin")  // heaviside step function
  {
    if (t <= params_->lifetime_ + 1.0e-11) {
      degr = 1.0;
    } else {
      degr = 0.0;
    }
  }
  else if (params_->degoption_ == "Exp" || params_->degoption_ == "ExpVar")  // exponential decrease
  {
    degr = exp(- t / params_->lifetime_ * log(2.0));
  }
  else if (params_->degoption_ == "Cos")  // transition zone with cos shape
  {
    if (t < 0.2 * params_->lifetime_ - 1.0e-11) {
      degr = 1.0;
    } else if (t < params_->lifetime_ - 1.0e-11) {
      degr = 0.5*(cos(PI*(t-0.2*params_->lifetime_)/(0.8*params_->lifetime_))+1.0);
    } else {
      degr = 0.0;
    }
  } else dserror("Degradation option not implemented! Valid options are Lin, Exp, ExpVar and Cos !");
}

/*----------------------------------------------------------------------*
 |  ElastinDegradation                            (private)        05/13|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::ElastinDegradation(LINALG::Matrix<1,3> coord, double& degr)
{
  if (*params_->elastindegrad_ == "Rectangle")
  {
    double funcz = 0.0;
    double z1 = -16.0;
    double z2 = -12.0;
    double z3 = 12.0;
    double z4 = 16.0;
    if (z1 < coord(2) && coord(2) < z2) {
      funcz = 0.5* (1.0 - cos((coord(2) - z1) / (z2-z1) * PI));
    } else if (z3 < coord(2) && coord(2) < z4) {
      funcz = 0.5* (1.0 + cos((coord(2) - z3) / (z4-z3) * PI));
    } else if (z2 <= coord(2) && coord(2) <= z3) {
      funcz = 1.0;
    }

    double funcphi = 0.0;
    double phi = atan2(coord(1),coord(0));
    double phi1 = -0.55 * PI;
    double phi2 = -0.5 * PI;
    double phi3 = -0.25 * PI;
    double phi4 = -0.2 * PI;
    if (phi1 < phi && phi < phi2) {
      funcphi = 0.5* (1.0 - cos((phi - phi1) / (phi2-phi1) * PI));
    } else if (phi3 < phi && phi < phi4) {
      funcphi = 0.5* (1.0 + cos((phi - phi3) / (phi4-phi3) * PI));
    } else if (phi2 <= phi && phi <= phi3) {
      funcphi = 1.0;
    }

    degr = 1.0 - funcz * funcphi;
  }
  else if (*params_->elastindegrad_ == "RectanglePlate")
  {
    double funcz = 0.0;
    double z1 = -0.5;
    double z2 = 0.0;
    double z3 = 2.0;
    double z4 = 2.5;
    if (z1 < coord(2) && coord(2) < z2) {
      funcz = 0.5* (1.0 - cos((coord(2) - z1) / (z2-z1) * PI));
    } else if (z3 < coord(2) && coord(2) < z4) {
      funcz = 0.5* (1.0 + cos((coord(2) - z3) / (z4-z3) * PI));
    } else if (z2 <= coord(2) && coord(2) <= z3) {
      funcz = 1.0;
    }

    double funcx = 0.0;
    double x1 = -2.5;
    double x2 = -2.0;
    double x3 = 2.0;
    double x4 = 2.5;
    if (x1 < coord(0) && coord(0) < x2) {
      funcx = 0.5* (1.0 - cos((coord(0) - x1) / (x2-x1) * PI));
    } else if (x3 < coord(0) && coord(0) < x4) {
      funcx = 0.5* (1.0 + cos((coord(0) - x3) / (x4-x3) * PI));
    } else if (x2 <= coord(0) && coord(0) <= x3) {
      funcx = 1.0;
    }

    degr = 1.0 - funcz * funcx;
  }
  else if (*params_->elastindegrad_ == "Wedge")
  {
    double funcz = 0.0;
    double z1 = -16.0; //-22.0; //-16.0;
    double z2 = -12.0; //-19.0; //-12.0;
    double z3 = 12.0; //19.0; //12.0;
    double z4 = 16.0; //22.0; //16.0;
    double phi = atan2(coord(1),coord(0));
    z1 = ((PI - abs(phi)) / PI *0.75 + 0.25)*z1;
    z2 = ((PI - abs(phi)) / PI *0.75 + 0.25)*z2;
    z3 = ((PI - abs(phi)) / PI *0.75 + 0.25)*z3;
    z4 = ((PI - abs(phi)) / PI *0.75 + 0.25)*z4;
    if (z1 < coord(2) && coord(2) < z2) {
      funcz = 0.5* (1.0 - cos((coord(2) - z1) / (z2-z1) * PI));
    } else if (z3 < coord(2) && coord(2) < z4) {
      funcz = 0.5* (1.0 + cos((coord(2) - z3) / (z4-z3) * PI));
    } else if (z2 <= coord(2) && coord(2) <= z3) {
      funcz = 1.0;
    }

    degr = 1.0 - funcz;
  }
  else if (*params_->elastindegrad_ == "Circles")
  {
    double radmin = 10.0;
    double radmax = 15.0;
    double func1 = 0.0;
    LINALG::Matrix<1,3> center1(true);
    center1(0) = 12.0; center1(1) = 0.0; center1(2) = 10.0;
    LINALG::Matrix<1,3> diff(coord);
    diff.Update(-1.0,center1,1.0);
    double rad1 = diff.Norm2();
    if (rad1 < radmin) {
      func1 = 1.0;
    } else if (rad1 < radmax) {
      func1 = 0.5* (1.0 - cos((rad1 - radmax) / (radmax-radmin) * PI));
    }
    double func2 = 0.0;
    LINALG::Matrix<1,3> center2(true);
    center2(0) = -12.0; center2(1) = 0.0; center2(2) = -10.0;
    diff.Update(coord);
    diff.Update(-1.0,center2,1.0);
    double rad2 = diff.Norm2();
    if (rad2 < radmin) {
      func2 = 1.0;
    } else if (rad2 < radmax) {
      func2 = 0.5* (1.0 - cos((rad2 - radmax) / (radmax-radmin) * PI));
    }
    double func = std::max(func1, func2);
    degr = 1.0 - func;
  }
}

/*----------------------------------------------------------------------*
 |  EvaluateImplicitAll                           (private)        12/11|
 *----------------------------------------------------------------------*
 evaluate stress and cmat for implicit integration
 driving force of massproduction is the total stress S
 */
void MAT::ConstraintMixture::EvaluateImplicitAll
(
  LINALG::Matrix<3,3> defgrd,
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  const int gp,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress,
  double dt,
  double time,
  LINALG::Matrix<4,1> massprod,
  LINALG::Matrix<4,1> massstress,
  double elastin_survival,
  double growthfactor
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
  LINALG::Matrix<4,1> prestretchcollagen = localprestretch_->at(gp);
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
  EvaluateStress(glstrain, gp, &cmattemp, &stressresidual, firstiter, time, elastin_survival);

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
    for (int id = 0; id < NUM_STRESS_3D; id++)
      DResidual(id,id) = 1.0;

    // for all 4 fiber families
    double stretch = prestretchcollagen(0)/actcollstretch(0);
    LINALG::Matrix<NUM_STRESS_3D,1> dstressdmass(true);
    LINALG::Matrix<NUM_STRESS_3D,1> dmassdstress(true);
    GradStressDMass(glstrain, &dstressdmass, Cinv, a1_->at(gp), stretch, J, dt, true);
    GradMassDStress(&dmassdstress, defgrd, Smatrix, a1_->at(gp), J, massstress(0), localhomstress_->at(gp)(0), actcollstretch(0), growthfactor);
    DResidual.MultiplyNT(-1.0,dstressdmass,dmassdstress,1.0);

    stretch = prestretchcollagen(1)/actcollstretch(1);
    dstressdmass.Scale(0.0);
    dmassdstress.Scale(0.0);
    GradStressDMass(glstrain, &dstressdmass, Cinv, a2_->at(gp), stretch, J, dt, true);
    GradMassDStress(&dmassdstress, defgrd, Smatrix, a2_->at(gp), J, massstress(1), localhomstress_->at(gp)(1), actcollstretch(1), growthfactor);
    DResidual.MultiplyNT(-1.0,dstressdmass,dmassdstress,1.0);

    stretch = prestretchcollagen(2)/actcollstretch(2);
    dstressdmass.Scale(0.0);
    dmassdstress.Scale(0.0);
    GradStressDMass(glstrain, &dstressdmass, Cinv, a3_->at(gp), stretch, J, dt, true);
    GradMassDStress(&dmassdstress, defgrd, Smatrix, a3_->at(gp), J, massstress(2), localhomstress_->at(gp)(2), actcollstretch(2), growthfactor);
    DResidual.MultiplyNT(-1.0,dstressdmass,dmassdstress,1.0);

    stretch = prestretchcollagen(3)/actcollstretch(3);
    dstressdmass.Scale(0.0);
    dmassdstress.Scale(0.0);
    GradStressDMass(glstrain, &dstressdmass, Cinv, a4_->at(gp), stretch, J, dt, true);
    GradMassDStress(&dmassdstress, defgrd, Smatrix, a4_->at(gp), J, massstress(3), localhomstress_->at(gp)(3), actcollstretch(3), growthfactor);
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

      MassProduction(gp, defgrd, stepstress, &massstress, 0.0, &massprod, growthfactor);
      history_->back().SetMass(gp, massprod);
      stressresidual.Scale(0.0);
      EvaluateStress(glstrain, gp, &cmattemp, &stressresidual, firstiter, time, elastin_survival);

      Residualtemp.Update(1.0,stepstress,-1.0,stressresidual);
    }
    //if (omega <= omegamin && Residualtemp.Norm2() > (1.0-0.5*omega)*Residual.Norm2())
    //  dserror("no damping coefficient found");

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
      std::cout << "1: " << massprod(0) << std::endl;
      std::cout << "2: " << massprod(1) << std::endl;
      std::cout << "3: " << massprod(2) << std::endl;
      std::cout << "4: " << massprod(3) << std::endl;
      dserror("negative mass production computed for at least one collagen fiber family!");
    }

  } //while loop
  if (localistep == maxstep && Residual.Norm2() > params_->abstol_*(*stress).NormInf())
    dserror("local Newton iteration did not converge %e", Residual.Norm2());

  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatelastic(true);
  EvaluateStress(glstrain, gp, &cmatelastic, stress, firstiter, time, elastin_survival);

  //--------------------------------------------------------------------------------------
  // compute cmat
  // right handside of the linear equations
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> RHS(cmatelastic);
  // left matrix of the linear equations
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> LM(true);
  for (int id = 0; id < NUM_STRESS_3D; id++) LM(id,id) = 1.0;

  // Fiber1
  double stretch = prestretchcollagen(0)/actcollstretch(0);
  LINALG::Matrix<NUM_STRESS_3D,1> dmassdstress(true);
  LINALG::Matrix<NUM_STRESS_3D,1> dmassdstretch(true);
  LINALG::Matrix<NUM_STRESS_3D,1> dstressdmass(true);
  GradMassDStretch(&dmassdstretch, defgrd, Smatrix, Cinv, a1_->at(gp), J, massstress(0), localhomstress_->at(gp)(0), actcollstretch(0), dt, growthfactor);
  GradMassDStress(&dmassdstress, defgrd, Smatrix, a1_->at(gp), J, massstress(0), localhomstress_->at(gp)(0), actcollstretch(0), growthfactor);
  GradStressDMass(glstrain, &dstressdmass, Cinv, a1_->at(gp), stretch, J, dt, true);
  RHS.MultiplyNT(2.0,dstressdmass,dmassdstretch,1.0);
  LM.MultiplyNT(-1.0,dstressdmass,dmassdstress,1.0);

  // Fiber2
  stretch = prestretchcollagen(1)/actcollstretch(1);
  GradMassDStretch(&dmassdstretch, defgrd, Smatrix, Cinv, a2_->at(gp), J, massstress(1), localhomstress_->at(gp)(1), actcollstretch(1), dt, growthfactor);
  GradMassDStress(&dmassdstress, defgrd, Smatrix, a2_->at(gp), J, massstress(1), localhomstress_->at(gp)(1), actcollstretch(1), growthfactor);
  GradStressDMass(glstrain, &dstressdmass, Cinv, a2_->at(gp), stretch, J, dt, true);
  RHS.MultiplyNT(2.0,dstressdmass,dmassdstretch,1.0);
  LM.MultiplyNT(-1.0,dstressdmass,dmassdstress,1.0);

  // Fiber3
  stretch = prestretchcollagen(2)/actcollstretch(2);
  GradMassDStretch(&dmassdstretch, defgrd, Smatrix, Cinv, a3_->at(gp), J, massstress(2), localhomstress_->at(gp)(2), actcollstretch(2), dt, growthfactor);
  GradMassDStress(&dmassdstress, defgrd, Smatrix, a3_->at(gp), J, massstress(2), localhomstress_->at(gp)(2), actcollstretch(2), growthfactor);
  GradStressDMass(glstrain, &dstressdmass, Cinv, a3_->at(gp), stretch, J, dt, true);
  RHS.MultiplyNT(2.0,dstressdmass,dmassdstretch,1.0);
  LM.MultiplyNT(-1.0,dstressdmass,dmassdstress,1.0);

  // Fiber4
  stretch = prestretchcollagen(3)/actcollstretch(3);
  GradMassDStretch(&dmassdstretch, defgrd, Smatrix, Cinv, a4_->at(gp), J, massstress(3), localhomstress_->at(gp)(3), actcollstretch(3), dt, growthfactor);
  GradMassDStress(&dmassdstress, defgrd, Smatrix, a4_->at(gp), J, massstress(3), localhomstress_->at(gp)(3), actcollstretch(3), growthfactor);
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
  LINALG::Matrix<3,3> defgrd,
  const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
  const int gp,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
  LINALG::Matrix<NUM_STRESS_3D,1>* stress,
  double dt,
  double time,
  double elastin_survival,
  double growthfactor
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
  double currmassdens = 0.0;
  double qdegrad = 0.0;
  Degradation(0.0, qdegrad);
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

    double homstress = localhomstress_->at(gp)(idfiber);
    double prestretchcollagen = localprestretch_->at(gp)(idfiber);
    double stretch = prestretchcollagen/actcollstretch(idfiber);
    LINALG::Matrix<NUM_STRESS_3D,1> stressfiber(true);
    LINALG::Matrix<NUM_STRESS_3D,1> stressresidual(true);
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatfiber(true);
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmattemp(true);
    double currmassdensfiber = 0.0;
    double currmassdenstemp = 0.0;
    double massprodfiber = 0.0;
    double massstressfiber = 0.0;
    EvaluateFiberFamily(C, gp, &cmatfiber, &stressfiber, a, &currmassdensfiber, firstiter, time, idfiber);
    // mass always corresponds to the current stress
    MassProductionSingleFiber(gp, defgrd, stressfiber, &massstressfiber, 0.0, &massprodfiber, a, idfiber, growthfactor);
    massprod(idfiber) = massprodfiber;
    history_->back().SetMass(gp, massprod);
    massstress(idfiber) = massstressfiber;
    // compute stresses for the computed mass
    EvaluateFiberFamily(C, gp, &cmattemp, &stressresidual, a, &currmassdenstemp, firstiter, time, idfiber);

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
      for (int id = 0; id < NUM_STRESS_3D; id++)
        DResidual(id,id) = 1.0;

      // linearisation of stress formula
      LINALG::Matrix<NUM_STRESS_3D,1> dstressdmass(true);
      GradStressDMass(glstrain, &dstressdmass, Cinv, a, stretch, J, dt, false);
      LINALG::Matrix<NUM_STRESS_3D,1> dmassdstress(true);
      GradMassDStress(&dmassdstress, defgrd, Smatrix, a, J, massstress(idfiber), homstress, actcollstretch(idfiber), growthfactor);
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
        MassProductionSingleFiber(gp, defgrd, stepstress, &massstressfiber, 0.0, &massprodfiber, a, idfiber, growthfactor);
        massprod(idfiber) = massprodfiber;
        history_->back().SetMass(gp, massprod);
        massstress(idfiber) = massstressfiber;
        // compute stresses for the computed mass
        stressresidual.Scale(0.0);
        EvaluateFiberFamily(C, gp, &cmattemp, &stressresidual, a, &currmassdenstemp, firstiter, time, idfiber);

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
        std::cout << idfiber+1 << ": " << massprod(idfiber) << std::endl;
        dserror("negative mass production computed for one collagen fiber family!");
      }

    } //while loop
    if (localistep == maxstep && Residual.Norm2() > params_->abstol_*stressfiber.NormInf())
      dserror("local Newton iteration did not converge %e", Residual.Norm2());

    currmassdensfiber = 0.0;
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatelastic(true);
    stressfiber.Scale(0.0);
    EvaluateFiberFamily(C, gp, &cmatelastic, &stressfiber, a, &currmassdensfiber, firstiter, time, idfiber);

    //--------------------------------------------------------------------------------------
    // compute cmat
    // right handside of the linear equations
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> RHS(cmatelastic);
    // left matrix of the linear equations
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> LM(true);
    for (int id = 0; id < NUM_STRESS_3D; id++) LM(id,id) = 1.0;

    LINALG::Matrix<NUM_STRESS_3D,1> dmassdstress(true);
    LINALG::Matrix<NUM_STRESS_3D,1> dmassdstretch(true);
    LINALG::Matrix<NUM_STRESS_3D,1> dstressdmass(true);
    GradMassDStretch(&dmassdstretch, defgrd, Smatrix, Cinv, a, J, massstress(idfiber), homstress, actcollstretch(idfiber), dt, growthfactor);
    GradMassDStress(&dmassdstress, defgrd, Smatrix, a, J, massstress(idfiber), homstress, actcollstretch(idfiber), growthfactor);
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
  LINALG::Matrix<NUM_STRESS_3D,1> Siso(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatiso(true);
  EvaluateElastin(C, &cmatiso, &Siso, time, &currmassdens, elastin_survival);
  (*stress) += Siso;
  (*cmat) += cmatiso;

  // smooth muscle cells
  LINALG::Matrix<NUM_STRESS_3D,1> Smus(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatmus(true);
  EvaluateMuscle(C, &cmatmus, &Smus, gp, &currmassdens);
  (*stress) += Smus;
  (*cmat) += cmatmus;

  // volumetric part
  LINALG::Matrix<NUM_STRESS_3D,1> Svol(true);
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatvol(true);
  EvaluateVolumetric(C, &cmatvol, &Svol, currmassdens, density);
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
  double qdegrad = 0.0;
  Degradation(0.0, qdegrad);

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

  double facS = 0.0;
  double temp = 0.0;
  LINALG::Matrix<NUM_STRESS_3D,1> Saniso(A);
  I4 = I4 * stretch*stretch;  // account for prestretch and stretch at deposition time
  EvaluateSingleFiberScalars(I4, temp, facS);
  facS = facS * stretch*stretch * qdegrad / density * dt;
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
  LINALG::Matrix<3,3> defgrd,
  LINALG::Matrix<3,3> Smatrix,
  LINALG::Matrix<3,1> a,
  double J,
  double massstress,
  double homstress,
  double actcollstretch,
  double growthfactor
)
{
  // include the case homstress = 0.0!!
  double homstressfixed = 1.0;
  if (homstress != 0.0) homstressfixed = homstress;

  LINALG::Matrix<3,3> Cmatrix(true);
  Cmatrix.MultiplyTN(defgrd,defgrd);

  LINALG::Matrix<3,1> Ca(true);
  LINALG::Matrix<3,1> CSCa(true);
  Ca.Multiply(Cmatrix,a);
  LINALG::Matrix<3,1> temp1(true);
  temp1.Multiply(Smatrix,Ca);
  CSCa.Multiply(Cmatrix,temp1);
  double fac = massprodbasal_ * growthfactor / homstressfixed / (J*J)
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
  LINALG::Matrix<3,3> defgrd,
  LINALG::Matrix<3,3> Smatrix,
  LINALG::Matrix<NUM_STRESS_3D,1> Cinv,
  LINALG::Matrix<3,1> a,
  double J,
  double massstress,
  double homstress,
  double actcollstretch,
  double dt,
  double growthfactor
)
{
  // include the case homstress = 0.0!!
  double homstressfixed = 1.0;
  if (homstress != 0.0) homstressfixed = homstress;

  LINALG::Matrix<3,3> Cmatrix(true);
  Cmatrix.MultiplyTN(defgrd,defgrd);

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
  (*derivative).Scale(0.5*massprodbasal_*growthfactor/homstressfixed);
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

/*----------------------------------------------------------------------*
 |  Return names of visualization data            (public)         03/13|
 *----------------------------------------------------------------------*/
void MAT::ConstraintMixture::VisNames(std::map<std::string,int>& names)
{
  std::string fiber = "MassStress";
  names[fiber] = 3;
  fiber = "Fiber1";
  names[fiber] = 3; // 3-dim vector
  fiber = "Fiber2";
  names[fiber] = 3; // 3-dim vector
  fiber = "referentialMassDensity";
  names[fiber] = 1;
  fiber = "CollagenMassDensity";
  names[fiber] = 3;
  fiber = "Prestretch";
  names[fiber] = 3;
  fiber = "Homstress";
  names[fiber] = 3;
  fiber = "MassProd";
  names[fiber] = 3;
  fiber = "growthfactor";
  names[fiber] = 1;
  fiber = "elastin_survival";
  names[fiber] = 1;
}

/*----------------------------------------------------------------------*
 |  Return visualization data                     (public)         03/13|
 *----------------------------------------------------------------------*/
bool MAT::ConstraintMixture::VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  if (name == "MassStress")
  {
    if ((int)data.size()!=3)
      dserror("size mismatch");
    LINALG::Matrix<3,1> temp(true);
    for (int iter=0; iter<numgp; iter++)
      temp.Update(1.0,vismassstress_->at(iter),1.0);
    data[0] = temp(0)/numgp;
    data[1] = temp(1)/numgp;
    data[2] = temp(2)/numgp;
  }
  else if (name == "Fiber1")
  {
    if ((int)data.size()!=3)
      dserror("size mismatch");
    LINALG::Matrix<3,1> a1 = a1_->at(0); // get a1 of first gp
    data[0] = a1(0);
    data[1] = a1(1);
    data[2] = a1(2);
  }
  else if (name == "Fiber2")
  {
    if ((int)data.size()!=3)
      dserror("size mismatch");
    LINALG::Matrix<3,1> a2 = a2_->at(0); // get a2 of first gp
    data[0] = a2(0);
    data[1] = a2(1);
    data[2] = a2(2);
  }
  else if (name == "referentialMassDensity")
  {
    if ((int)data.size()!=1)
      dserror("size mismatch");
    double temp = 0.0;
    for (int iter=0; iter<numgp; iter++)
      temp += refmassdens_->at(iter);
    data[0] = temp/numgp;
  }
  else if (name == "CollagenMassDensity")
  {
    if ((int)data.size()!=3)
      dserror("size mismatch");
    LINALG::Matrix<3,1> temp(true);
    for (int iter=0; iter<numgp; iter++)
      temp.Update(1.0,visrefmassdens_->at(iter),1.0);
    data[0] = temp(0)/numgp;
    data[1] = temp(1)/numgp;
    data[2] = temp(2)/numgp;
  }
  else if (name == "Prestretch")
  {
    if ((int)data.size()!=3)
      dserror("size mismatch");
    LINALG::Matrix<3,1> temp(true);
    for (int iter=0; iter<numgp; iter++)
      temp.Update(1.0,GetPrestretch(iter),1.0);
    data[0] = temp(0)/numgp;
    data[1] = temp(1)/numgp;
    data[2] = temp(2)/numgp;
  }
  else if (name == "Homstress")
  {
    if ((int)data.size()!=3)
      dserror("size mismatch");
    LINALG::Matrix<3,1> temp(true);
    for (int iter=0; iter<numgp; iter++)
      temp.Update(1.0,GetHomstress(iter),1.0);
    data[0] = temp(0)/numgp;
    data[1] = temp(1)/numgp;
    data[2] = temp(2)/numgp;
  }
  else if (name == "MassProd")
  {
    if ((int)data.size()!=3)
      dserror("size mismatch");
    LINALG::Matrix<4,1> temp(true);
    int sizehistory = history_->size();
    for (int iter=0; iter<numgp; iter++)
    {
      LINALG::Matrix<4,1> temp_loc(true);
      history_->at(sizehistory-2).GetMass(iter,&temp_loc);
      //history_->at(0).GetMass(iter,&temp_loc);
      //history_->at(0).GetStretches(iter,&temp_loc);
      temp.Update(1.0,temp_loc,1.0);
    }
    data[0] = temp(0)/numgp;
    data[1] = temp(1)/numgp;
    data[2] = temp(2)/numgp;
  }
  else if (name == "growthfactor")
  {
    if ((int)data.size()!=1)
      dserror("size mismatch");
    int eleLID = DRT::Problem::Instance()->GetDis("structure")->ElementColMap()->LID(eleID);
    data[0] = params_->GetParameter(params_->growthfactor,eleLID);
  }
  else if (name == "elastin_survival")
  {
    if ((int)data.size()!=1)
      dserror("size mismatch");
    int eleLID = DRT::Problem::Instance()->GetDis("structure")->ElementColMap()->LID(eleID);
    if (*params_->elastindegrad_ == "InvEla")
      data[0] = params_->GetParameter(params_->elastin_survival,eleLID);
    else if (*params_->elastindegrad_ == "Rectangle" || *params_->elastindegrad_ == "RectanglePlate"
        || *params_->elastindegrad_ == "Wedge" || *params_->elastindegrad_ == "Circles")
    {
      DRT::Element* myele = DRT::Problem::Instance()->GetDis("structure")->gElement(eleID);
      DRT::Node** mynodes = myele->Nodes();
      for (int idnodes = 0; idnodes< myele->NumNode();idnodes++)
      {
        DRT::Node* locnode = mynodes[idnodes];
        double elastin_survival = 0.0;
        LINALG::Matrix<1,3> point_refe;
        point_refe(0) = locnode->X()[0]; point_refe(1) = locnode->X()[1]; point_refe(2) = locnode->X()[2];
        ElastinDegradation(point_refe, elastin_survival);
        data[0] += elastin_survival;
      }
      data[0] = data[0] / myele->NumNode();
    }
    else
      data[0] = 1.0;
  }
  else
  {
    return false;
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Debug output to gmsh-file                                      01/11|
 *----------------------------------------------------------------------*
 this needs to be copied to STR::TimInt::OutputStep() to enable debug output
 {
   discret_->SetState("displacement",Dis());
   MAT::ConstraintMixtureOutputToGmsh(discret_, StepOld(), 1);
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
  filename << filebase << "_massdensity" << std::setw(3) << std::setfill('0') << timestep
      << std::setw(2) << std::setfill('0') << iter << ".pos";
  std::ofstream f_system(filename.str().c_str());
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" Time: " << timestep << " Iter: " << iter << " \" {" << std::endl;
  gmshfilecontent.precision(10);

  // file for prestretch
  std::stringstream filename_pre;
  filename_pre << filebase << "_circollagendens" << std::setw(3) << std::setfill('0') << timestep
      << std::setw(2) << std::setfill('0') << iter << ".pos";
  std::ofstream f_system_pre(filename_pre.str().c_str());
  std::stringstream gmshfilecontent_pre;
  gmshfilecontent_pre << "View \" Time: " << timestep << " Iter: " << iter << " \" {" << std::endl;
  gmshfilecontent_pre.precision(10);


  for (int iele=0; iele<dis->NumMyColElements(); ++iele)
  {
    const DRT::Element* actele = dis->lColElement(iele);

    // build current configuration
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    actele->LocationVector(*dis,lm,lmowner,lmstride);
    Teuchos::RCP<const Epetra_Vector> disp = dis->GetState("displacement");
    std::vector<double> mydisp(lm.size(),0);
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

    Teuchos::RCP<MAT::Material> mat = actele->Material();
    MAT::ConstraintMixture* grow = static_cast <MAT::ConstraintMixture*>(mat.get());

    // material plot at gauss points
    int ngp = grow->Geta1()->size();

    // update element geometry
    const int numnode = actele->NumNode();
    const int numdof = 3;
    Epetra_SerialDenseMatrix xcurr(numnode,3);  // material coord. of element
    for (int i=0; i<numnode; ++i)
    {
      xcurr(i,0) = actele->Nodes()[i]->X()[0] + mydisp[i*numdof+0];
      xcurr(i,1) = actele->Nodes()[i]->X()[1] + mydisp[i*numdof+1];
      xcurr(i,2) = actele->Nodes()[i]->X()[2] + mydisp[i*numdof+2];
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
      break;
    }

    const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule_);

    for (int gp = 0; gp < ngp; ++gp)
    {
      DRT::UTILS::shape_function_3D(funct, intpoints.qxg[gp][0], intpoints.qxg[gp][1], intpoints.qxg[gp][2], distype);
      Epetra_SerialDenseMatrix point(1,3);
      point.Multiply('T','N',1.0,funct,xcurr,0.0);

      // write mandel stress
      //LINALG::Matrix<3,1> mandelgp = grow->GetHomstress(gp);
      double mandelgp = grow->GetMassDensity(gp);
      gmshfilecontent << "SP(" << std::scientific << point(0,0) << ",";
      gmshfilecontent << std::scientific << point(0,1) << ",";
      gmshfilecontent << std::scientific << point(0,2) << ")";
      gmshfilecontent << "{" << std::scientific
      << mandelgp
      << "};" << std::endl;

      // write prestretch
      //LINALG::Matrix<3,1> prestretchgp = grow->GetPrestretch(gp);
      LINALG::Matrix<3,1> prestretchgp = grow->GetMassDensityCollagen(gp);
      gmshfilecontent_pre << "SP(" << std::scientific << point(0,0) << ",";
      gmshfilecontent_pre << std::scientific << point(0,1) << ",";
      gmshfilecontent_pre << std::scientific << point(0,2) << ")";
      gmshfilecontent_pre << "{" << std::scientific
      << prestretchgp(0)
      << "};" << std::endl;

    }
  }

  gmshfilecontent << "};" << std::endl;
  f_system << gmshfilecontent.str();
  f_system.close();

  gmshfilecontent_pre << "};" << std::endl;
  f_system_pre << gmshfilecontent_pre.str();
  f_system_pre.close();

  return;
}
