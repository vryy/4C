/*----------------------------------------------------------------------*/
/*!
 \file growth_scd.cpp

 \brief

This file contains routines for an integration point based and scalar dependend volumetric growth law.
It is derived from the growth laws implemented in growth_ip.cpp and additional adds the scalar dependency, e.g. nutrients.

 <pre>
   Maintainer: Moritz Thon
               thon@lnm.mw.tum.de
               http://www.mhpc.mw.tum.de
               089 - 289-10364
 </pre>
 *----------------------------------------------------------------------*/


#include "growth_scd.H"
#include "growth_law.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"   // for function Factory in Unpack
#include "../drt_lib/drt_utils_factory.H"   // for function Factory in Unpack


/*----------------------------------------------------------------------*
 |                                                           vuong 06/11  |
 *----------------------------------------------------------------------*/
MAT::PAR::GrowthScd::GrowthScd(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Growth(matdata),
  rearate_(matdata->GetDouble("REARATE")),
  satcoeff_(matdata->GetDouble("SATCOEFF")),
  growthcoupl_(matdata->Get<std::string>("GROWTHCOUPL"))
{
  if (growthlaw_->MaterialType() == INPAR::MAT::m_growth_linear or growthlaw_->MaterialType() == INPAR::MAT::m_growth_exponential)
  {
    if (rearate_ <= 0)
        dserror("You need to choose a positive reaction rate!");
    if (satcoeff_ < 0)
        dserror("You need to choose a non-negative saturation coefficient!");
  }
}


Teuchos::RCP<MAT::Material> MAT::PAR::GrowthScd::CreateMaterial()
{
  Teuchos::RCP<MAT::Material> mat;

  switch (growthlaw_->MaterialType())
  {
  case INPAR::MAT::m_growth_linear:
  case INPAR::MAT::m_growth_exponential:
    mat = Teuchos::rcp(new MAT::GrowthScd(this));
    break;
  case INPAR::MAT::m_growth_ac:
    mat = Teuchos::rcp(new MAT::GrowthScdAC(this));
    break;
  default:
    dserror("The growth law you have chosen is not valid one!");
    mat = Teuchos::null;
    break;
  }

  return mat;
}


MAT::GrowthScdType MAT::GrowthScdType::instance_;

DRT::ParObject* MAT::GrowthScdType::Create( const std::vector<char> & data )
{
  MAT::GrowthScd* grow = new MAT::GrowthScd();
  grow->Unpack(data);
  return grow;
}

/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  vuong 02/10|
 *----------------------------------------------------------------------*/
MAT::GrowthScd::GrowthScd()
  : GrowthMandel(),
    detFe_(Teuchos::null), //initialized in GrowthScd::Unpack
    dtheta_(Teuchos::null), //initialized in GrowthScd::Unpack
    concentration_(-1.0),
    stressgrowthfunc(-1.0),
    paramsScd_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)    vuong  02/10|
 *----------------------------------------------------------------------*/
MAT::GrowthScd::GrowthScd(MAT::PAR::GrowthScd* params)
  : GrowthMandel(params),
    detFe_(Teuchos::null), //initialized in GrowthScd::Setup
    dtheta_(Teuchos::null), //initialized in GrowthScd::Setup
    concentration_(-1.0),
    stressgrowthfunc(-1.0),
    paramsScd_(params)
{
}

/*----------------------------------------------------------------------*
 |  ResetAll                                      (public)         11/12|
 *----------------------------------------------------------------------*/
void MAT::GrowthScd::ResetAll(const int numgp)
{
  for (int j=0; j<numgp; ++j)
  {
    detFe_->at(j) = 1.0;
    dtheta_->at(j) = 0.0;
  }
  concentration_= -1.0;
  stressgrowthfunc = -1.0;

  MAT::GrowthMandel::ResetAll(numgp);
}

/*----------------------------------------------------------------------*
 |  Pack                                          (public)    vuong  02/10|
 *----------------------------------------------------------------------*/
void MAT::GrowthScd::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matid
  int matid = -1;

  MAT::PAR::GrowthScd* params=Parameter();
  if (params != NULL) matid = params->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);

  int numgp=0;
  if (isinit_)
  {
    numgp = dtheta_->size();   // size is number of gausspoints
  }
  AddtoPack(data,numgp);
  // Pack internal variables
  for (int gp = 0; gp < numgp; ++gp)
  {
    AddtoPack(data,detFe_->at(gp));
    AddtoPack(data,dtheta_->at(gp));
  }

  // Pack base class material
  GrowthMandel::Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)   vuong 02/10|
 *----------------------------------------------------------------------*/
void MAT::GrowthScd::Unpack(const std::vector<char>& data)
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

  paramsScd_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        paramsScd_ = static_cast<MAT::PAR::GrowthScd*>(mat);
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

  // unpack growth internal variables
  detFe_ = Teuchos::rcp(new std::vector<double> (numgp));
  dtheta_ = Teuchos::rcp(new std::vector<double> (numgp));
  for (int gp = 0; gp < numgp; ++gp) {
    double a;
    ExtractfromPack(position,data,a);
    detFe_->at(gp) = a;
    ExtractfromPack(position,data,a);
    dtheta_->at(gp) = a;
  }

  // extract base class material
  std::vector<char> basedata(0);
  GrowthMandel::ExtractfromPack(position,data,basedata);
  GrowthMandel::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 |  Setup                                         (public) vuong  02/10|
 *----------------------------------------------------------------------*/
void MAT::GrowthScd::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  detFe_ = Teuchos::rcp(new std::vector<double> (numgp,1.0));
  dtheta_ = Teuchos::rcp(new std::vector<double> (numgp,0.0));

  //setup base class
  GrowthMandel::Setup(numgp, linedef);
  return;
}


/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)  vuong   02/10|
 *----------------------------------------------------------------------*
 The deformation gradient is decomposed into an elastic and growth part:
     F = Felastic * F_g
 Only the elastic part contributes to the stresses, thus we have to
 compute the elastic Cauchy Green Tensor Cdach and elastic 2PK stress Sdach.
 */
void MAT::GrowthScd::Evaluate
(
    const LINALG::Matrix<3,3>* defgrd,
    const LINALG::Matrix<6,1>* glstrain,
    Teuchos::ParameterList& params,
    LINALG::Matrix<6,1>* stress,
    LINALG::Matrix<6,6>* cmat,
    const int eleGID
)
{
  // get gauss point number
  const int gp = params.get<int>("gp",-1);
  if (gp == -1) dserror("no Gauss point number provided in material");

  //get pointer vector containing the mean scalar values
  Teuchos::RCP<std::vector<double> > scalars = params.get< Teuchos::RCP<std::vector<double> > >("mean_concentrations",Teuchos::null);
  //in this growth law we always assume the first scalar to induce growth!
  if (scalars != Teuchos::null)
    concentration_ = scalars->at(0);

  //if(concentration_==-1.0) dserror("no scalar concentration provided for concentration dependent growth law");

  GrowthMandel::Evaluate(defgrd,glstrain,params,stress,cmat,eleGID);

  // build identity tensor I
  LINALG::Matrix<NUM_STRESS_3D,1> Id(true);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;
  // right Cauchy-Green Tensor  C = 2 * E + I
  LINALG::Matrix<NUM_STRESS_3D,1> C(*glstrain);
  C.Scale(2.0);
  C += Id;

  const double theta = theta_->at(gp);
  // elastic right Cauchy-Green Tensor Cdach = F_g^-T C F_g^-1
  LINALG::Matrix<NUM_STRESS_3D,1> Cdach(C);
  Cdach.Scale(1.0/theta/theta);

  // determinate of F_e (necessary for scatra/nutrientconsumption)
  detFe_->at(gp) =   pow( + Cdach(0)*(Cdach(1)*Cdach(2)-Cdach(4)*Cdach(4))
                          - Cdach(3)*(Cdach(3)*Cdach(2)-Cdach(5)*Cdach(4))
                          + Cdach(5)*(Cdach(3)*Cdach(4)-Cdach(5)*Cdach(1)),
                          1/2);

  //store dtheta
  const double dt = params.get<double>("delta time",-1.0);
  for (unsigned i=0; i<theta_->size(); i++)
    (*dtheta_)[i] = ((*theta_)[i]-(*thetaold_)[i])/dt;
}

/*----------------------------------------------------------------------*
 |  Evaluate growth function                        (protected)  vuong   02/10|
 *----------------------------------------------------------------------*/
void MAT::GrowthScd::EvaluateGrowthFunction
(
    double & growthfunc, // (o)
    double traceM,       // (i)
    double theta         // (i)
)
{
  //call stress based growth law
  GrowthMandel::EvaluateGrowthFunction(growthfunc,traceM,theta);
  stressgrowthfunc = growthfunc;

  MAT::PAR::GrowthScd* params=Parameter();
  // parameters
  const double rearate = params->rearate_;
  const double satcoeff = params->satcoeff_;
  const std::string* growthcoupl =params->growthcoupl_;

  if (*growthcoupl == "ScaleConc")
    // scale with concentration dependent factor
    growthfunc = rearate * concentration_/(satcoeff+concentration_) * growthfunc;

  else if (*growthcoupl == "StressRed")
    // reduce the growth due to scalar transport because of the presence of stresses (biofilm)
   growthfunc = rearate * concentration_/(satcoeff+concentration_) - abs(growthfunc);

  else
    dserror ("The chosen coupling law between stress dependent growth and reaction dependent growth is not implemented");

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate derivative of growth function       (protected)  vuong 02/10|
 *----------------------------------------------------------------------*/
void MAT::GrowthScd::EvaluateGrowthFunctionDerivTheta
(
    double & dgrowthfunctheta,
    double traceM,
    double theta,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& Cdach,
    const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmatelastic
)
{
  //call stress based growth law
  GrowthMandel::EvaluateGrowthFunctionDerivTheta(dgrowthfunctheta,traceM,theta,Cdach,cmatelastic);

  MAT::PAR::GrowthScd* params=Parameter();
  // parameters
  const double rearate = params->rearate_;
  const double satcoeff = params->satcoeff_;
  const std::string* growthcoupl =params->growthcoupl_;

  if (*growthcoupl == "ScaleConc")
    dgrowthfunctheta = rearate * concentration_/(satcoeff+concentration_) * dgrowthfunctheta;

  if (*growthcoupl == "StressRed")
  {
    if (stressgrowthfunc==0.0) dgrowthfunctheta = 0.0;
    else dgrowthfunctheta = - abs(stressgrowthfunc) / stressgrowthfunc * dgrowthfunctheta;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate derivative of growth function       (protected) vuong 02/10|
 *----------------------------------------------------------------------*/
void MAT::GrowthScd::EvaluateGrowthFunctionDerivC
(
    LINALG::Matrix<NUM_STRESS_3D, 1>& dgrowthfuncdC,
    double traceM,
    double theta,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& C,
    const LINALG::Matrix<NUM_STRESS_3D, 1>& S,
    const LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat
)
{
  //call stress based growth law
  GrowthMandel::EvaluateGrowthFunctionDerivC(dgrowthfuncdC,traceM,theta,C,S,cmat);

  // parameters
  MAT::PAR::GrowthScd* params=Parameter();

  const double rearate = params->rearate_;
  const double satcoeff = params->satcoeff_;
  const std::string* growthcoupl =params->growthcoupl_;

  if (*growthcoupl == "ScaleConc")
    // scale with concentration dependent factor
    dgrowthfuncdC.Scale(rearate * concentration_/(satcoeff+concentration_));

  if (*growthcoupl == "StressRed")
  {
    if (stressgrowthfunc==0) dgrowthfuncdC.Scale(0.0);
    else dgrowthfuncdC.Scale(- abs(stressgrowthfunc) / stressgrowthfunc);
  }
  return;
}






MAT::GrowthScdACType MAT::GrowthScdACType::instance_;

DRT::ParObject* MAT::GrowthScdACType::Create( const std::vector<char> & data )
{
  MAT::GrowthScdAC* grow = new MAT::GrowthScdAC();
  grow->Unpack(data);
  return grow;
}

/*----------------------------------------------------------------------*
 |  Constructor                                               Thon 11/14|
 *----------------------------------------------------------------------*/
MAT::GrowthScdAC::GrowthScdAC()
  : GrowthBasic(),
    concentrations_(Teuchos::null),
    paramsScdAC_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                                          Thon 11/14|
 *----------------------------------------------------------------------*/
MAT::GrowthScdAC::GrowthScdAC(MAT::PAR::GrowthScd* params)
  : GrowthBasic(params),
    concentrations_(Teuchos::null),
    paramsScdAC_(params)
{
}

/*----------------------------------------------------------------------*
 |  ResetAll                                                  Thon 11/14|
 *----------------------------------------------------------------------*/
void MAT::GrowthScdAC::ResetAll(const int numgp)
{
  concentrations_ = Teuchos::rcp(new std::vector<double> (10,0.0));

  MAT::Growth::ResetAll(numgp);
}

/*----------------------------------------------------------------------*
 |  Pack                                                      Thon 11/14|
 *----------------------------------------------------------------------*/
void MAT::GrowthScdAC::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matid
  int matid = -1;

  MAT::PAR::GrowthScd* params=Parameter();
  if (params != NULL) matid = params->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);

  int numscal=0;
  if (isinit_)
  {
    numscal = concentrations_->size();   // size is number of gausspoints
  }
  AddtoPack(data,numscal);
  // Pack internal variables
  for (int sc = 0; sc < numscal; ++sc)
  {
    AddtoPack(data,concentrations_->at(sc));
  }

  // Pack base class material
  Growth::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Setup                                                     Thon 11/14|
 *----------------------------------------------------------------------*/
void MAT::GrowthScdAC::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  concentrations_ = Teuchos::rcp(new std::vector<double> (10,0.0)); //this is just a dummy, since we overwrite this pointer in MAT::GrowthScdAC::Evaluate anyway

  //setup base class
  Growth::Setup(numgp, linedef);
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack                                                    Thon 11/14|
 *----------------------------------------------------------------------*/
void MAT::GrowthScdAC::Unpack(const std::vector<char>& data)
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

  paramsScdAC_=NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        paramsScdAC_ = dynamic_cast<MAT::PAR::GrowthScd*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  int numscal;
  ExtractfromPack(position,data,numscal);
  if (numscal == 0){ // no history data to unpack
    isinit_=false;
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d",data.size(),position);
    return;
  }

  // unpack growth internal variables
  concentrations_ = Teuchos::rcp(new std::vector<double> (numscal));
  for (int sc = 0; sc < numscal; ++sc) {
    double a;
    ExtractfromPack(position,data,a);
    concentrations_->at(sc) = a;
  }

  // extract base class material
  std::vector<char> basedata(0);
  Growth::ExtractfromPack(position,data,basedata);
  Growth::Unpack(basedata);


  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate Material                                         Thon 11/14|
 *----------------------------------------------------------------------*/
void MAT::GrowthScdAC::Evaluate
(
    const LINALG::Matrix<3,3>* defgrd,
    const LINALG::Matrix<6,1>* glstrain,
    Teuchos::ParameterList& params,
    LINALG::Matrix<6,1>* stress,
    LINALG::Matrix<6,6>* cmat,
    const int eleGID
)
{
  // get gauss point number
  const int gp = params.get<int>("gp",-1);
  if (gp == -1)
    dserror("no Gauss point number provided in material");

  //get pointer vector containing the mean scalar values
  concentrations_ = params.get< Teuchos::RCP<std::vector<double> > >("mean_concentrations",Teuchos::rcp(new std::vector<double> (10,0.0)) ); //TODO: (Thon) there must be a better way to catch the first call of the structure evaluate

//  if (concentrations_ == Teuchos::null)
//    dserror("getting the mean concentrations failed!");

  GrowthBasic::Evaluate(defgrd,glstrain,params,stress,cmat,eleGID);
}

/*----------------------------------------------------------------------*
 |  Calculate the volumetric growth parameter                 Thon 11/14|
 *----------------------------------------------------------------------*/
double MAT::GrowthScdAC::CalculateTheta( const double J )
{
  double theta = Parameter()->growthlaw_->CalculateTheta( concentrations_ , J );

  return theta;
}

/*--------------------------------------------------------------------------------*
 |  Calculate the volumetric growth derived w.r.t. chauy-green strains  Thon 01/15|
 *--------------------------------------------------------------------------------*/
void MAT::GrowthScdAC::CalculateThetaDerivC( LINALG::Matrix<3,3>& dThetadC , const LINALG::Matrix<3,3>& C , const double J )
{
  Parameter()->growthlaw_->CalculateThetaDerivC( dThetadC , C , concentrations_ , J );
}
