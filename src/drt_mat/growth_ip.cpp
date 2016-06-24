/*----------------------------------------------------------------------*/
/*!
 \file growth_ip.cpp

 \brief This file contains routines for an integration point based isotropic volumetric growth law.

\level 3

 <pre>
   \maintainer Moritz Thon
               thon@lnm.mw.tum.de
               http://www.mhpc.mw.tum.de
               089 - 289-10364
 </pre>
 *----------------------------------------------------------------------*/


#include "growth_ip.H"
#include "growth_law.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_utils_factory.H"


/*----------------------------------------------------------------------------*/
MAT::PAR::Growth::Growth(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  idmatelastic_(matdata->GetInt("IDMATELASTIC")),
  idgrowthlaw_(matdata->GetInt("GROWTHLAW")),
  starttime_(matdata->GetDouble("STARTTIME")),
  endtime_(matdata->GetDouble("ENDTIME"))
{
  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (DRT::Problem::Instance(probinst)->Materials() == Teuchos::null)
    dserror("Sorry dude, cannot work out problem instance.");
  // yet another safety check
  if (DRT::Problem::Instance(probinst)->Materials()->Num() == 0)
    dserror("Sorry dude, no materials defined.");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat = DRT::Problem::Instance(probinst)->Materials()->ById(idgrowthlaw_);

  switch (curmat->Type())
  {
  case INPAR::MAT::m_growth_linear:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::GrowthLawLinear(curmat));
    MAT::PAR::GrowthLawLinear* params = static_cast<MAT::PAR::GrowthLawLinear*>(curmat->Parameter());
    growthlaw_ = params->CreateGrowthLaw();
    break;
  }
  case INPAR::MAT::m_growth_exponential:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::GrowthLawExp(curmat));
    MAT::PAR::GrowthLawExp* params = static_cast<MAT::PAR::GrowthLawExp*>(curmat->Parameter());
    growthlaw_ = params->CreateGrowthLaw();
    break;
  }
  case INPAR::MAT::m_growth_ac:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::GrowthLawAC(curmat));
      MAT::PAR::GrowthLawAC* params = static_cast<MAT::PAR::GrowthLawAC*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
  case INPAR::MAT::m_growth_ac_radial:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::GrowthLawACRadial(curmat));
      MAT::PAR::GrowthLawACRadial* params = static_cast<MAT::PAR::GrowthLawACRadial*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
  case INPAR::MAT::m_growth_ac_radial_refconc:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::GrowthLawACRadialRefConc(curmat));
      MAT::PAR::GrowthLawACRadialRefConc* params = static_cast<MAT::PAR::GrowthLawACRadialRefConc*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
  case INPAR::MAT::m_growth_biofilm:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::GrowthLawBiofilm(curmat));
      MAT::PAR::GrowthLawBiofilm* params = static_cast<MAT::PAR::GrowthLawBiofilm*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
  case INPAR::MAT::m_growth_const:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::GrowthLawConst(curmat));
      MAT::PAR::GrowthLawConst* params = static_cast<MAT::PAR::GrowthLawConst*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
  default:
    dserror("unknown material type %d", curmat->Type());
    break;
  }
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::Growth::CreateMaterial()
{
  Teuchos::RCP<MAT::Material> mat;

  switch (growthlaw_->MaterialType())
  {
  case INPAR::MAT::m_growth_linear:
  case INPAR::MAT::m_growth_exponential:
  case INPAR::MAT::m_growth_const:
    mat = Teuchos::rcp(new MAT::GrowthMandel(this));
    break;
  default:
    dserror("The growth law you have chosen is not valid for the standard volumetric growth material");
    mat = Teuchos::null;
    break;
  }

  return mat;
}

/*----------------------------------------------------------------------------*/
MAT::Growth::Growth()
  : theta_(Teuchos::null),
    isinit_(false),
    params_(NULL),
    matelastic_(Teuchos::null),
    thetaold_(Teuchos::null),
    histdata_()
{
}


/*----------------------------------------------------------------------------*/
MAT::Growth::Growth(MAT::PAR::Growth* params)
  : theta_(Teuchos::null),
    isinit_(false),
    params_(params),
    matelastic_(Teuchos::null),
    thetaold_(Teuchos::null),
    histdata_()
{
}

/*----------------------------------------------------------------------------*/
void MAT::Growth::Pack(DRT::PackBuffer& data) const
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
    numgp = theta_->size();   // size is number of gausspoints
  }
  AddtoPack(data,numgp);
  // Pack internal variables
  for (int gp = 0; gp < numgp; ++gp)
  {
    AddtoPack(data,thetaold_->at(gp));
    AddtoPack(data,theta_->at(gp));
  }

  AddtoPack(data,histdata_);

  // Pack data of elastic material
  if (matelastic_!=Teuchos::null) {
    matelastic_->Pack(data);
  }

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::Growth::Unpack(const std::vector<char>& data)
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
        params_ = dynamic_cast<MAT::PAR::Growth*>(mat);
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
  theta_ = Teuchos::rcp(new std::vector<double> (numgp));
  thetaold_ = Teuchos::rcp(new std::vector<double> (numgp));
  for (int gp = 0; gp < numgp; ++gp)
  {
    double a;
    ExtractfromPack(position,data,a);
    thetaold_->at(gp) = a;
    ExtractfromPack(position,data,a);
    theta_->at(gp) = a;
  }

  std::vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  histdata_.Unpack(tmp);

  // Unpack data of elastic material (these lines are copied from drt_element.cpp)
  std::vector<char> dataelastic;
  ExtractfromPack(position,data,dataelastic);
  if (dataelastic.size()>0)
  {
    DRT::ParObject* o = DRT::UTILS::Factory(dataelastic);  // Unpack is done here
    MAT::So3Material* matel = dynamic_cast<MAT::So3Material*>(o);
    if (matel==NULL)
      dserror("failed to unpack elastic material");
    matelastic_ = Teuchos::rcp(matel);
  }
  else
    matelastic_ = Teuchos::null;

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::Growth::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  if (isinit_)
    dserror("This function should just be called, if the material is jet not initialized.");

  theta_ = Teuchos::rcp(new std::vector<double> (numgp));
  thetaold_ = Teuchos::rcp(new std::vector<double> (numgp));
  for (int j=0; j<numgp; ++j)
  {
    theta_->at(j) = 1.0;
    thetaold_->at(j) = 1.0;
  }

  // Setup of elastic material
  matelastic_ = Teuchos::rcp_dynamic_cast<MAT::So3Material>(MAT::Material::Factory(params_->idmatelastic_));
  matelastic_->Setup(numgp, linedef);

  // Setup of history container with empty data
  std::map<int,std::vector<double> > data;
  histdata_.Add("thetaold",data);

  isinit_ = true;
  return;
}

/*----------------------------------------------------------------------------*/
void MAT::Growth::ResetAll(const int numgp)
{
  for (int j=0; j<numgp; ++j)
  {
    thetaold_->at(j) = 1.0;
    theta_->at(j) = 1.0;
  }

  // overwrite history with empty data
  std::map<int,std::vector<double> > data;
  histdata_.Add("thetaold",data);

  matelastic_->ResetAll(numgp);
}

/*----------------------------------------------------------------------------*/
void MAT::Growth::Update()
{
  const int histsize = theta_->size();
  for (int i=0; i<histsize; i++)
  {
    thetaold_->at(i) = theta_->at(i);
  }

  matelastic_->Update();
}

/*----------------------------------------------------------------------------*/
void MAT::Growth::ResetStep()
{
  dserror("you need a backup of thetaold_ to be able to reset!");
  matelastic_->ResetStep();
}

void MAT::Growth::StoreHistory(int timestep)
{
  std::map<int,std::vector<double> >* access;
  access = histdata_.GetMutable<std::map<int,std::vector<double> > >("thetaold");
  (*access)[timestep] = *thetaold_;
}

void MAT::Growth::SetHistory(int timestep)
{
  const std::map<int,std::vector<double> >* read;
  read = histdata_.Get<std::map<int,std::vector<double> > >("thetaold");

  if (read->find(timestep) != read->end())
    *thetaold_ = read->at(timestep);
  else
    dserror("there is no data to reset for step %d",timestep);
}

/*----------------------------------------------------------------------------*/
void MAT::Growth::EvaluateElastic(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain,
    LINALG::Matrix<6, 1>* stress,
    LINALG::Matrix<6, 6>* cmat,
    Teuchos::ParameterList& params,
    const int eleGID)
{
  Matelastic()->Evaluate(defgrd, glstrain, params, stress, cmat, eleGID);
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthMandel::VisNames(std::map<std::string,int>& names)
{
  std::string fiber = "Theta";
  names[fiber] = 1;
  fiber = "Mandel";
  names[fiber] = 1;
  Matelastic()->VisNames(names);
}

/*----------------------------------------------------------------------------*/
bool MAT::GrowthMandel::VisData(const std::string& name, std::vector<double>& data, int numgp , int eleID)
{
  if (name == "Theta")
  {
    if ((int)data.size()!=1)
      dserror("size mismatch");
    double temp = 0.0;
    for (int iter=0; iter<numgp; iter++)
      temp += theta_->at(iter);
    data[0] = temp/numgp;
  }
  else if (name == "Mandel")
  {
    if ((int)data.size()!=1)
      dserror("size mismatch");
    double temp = 0.0;
    for (int iter=0; iter<numgp; iter++)
      temp += mandel_->at(iter);
    data[0] = temp/numgp;
  }
  else
  {
    return Matelastic()->VisData(name, data, numgp, eleID);
  }
  return true;
}


MAT::GrowthMandelType MAT::GrowthMandelType::instance_;

/*----------------------------------------------------------------------------*/
DRT::ParObject* MAT::GrowthMandelType::Create( const std::vector<char> & data )
{
  MAT::GrowthMandel* grow = new MAT::GrowthMandel();
  grow->Unpack(data);
  return grow;
}

/*----------------------------------------------------------------------------*/
MAT::GrowthMandel::GrowthMandel()
  : Growth(),
    mandel_(Teuchos::null),
    paramsMandel_(NULL)
{
}

/*----------------------------------------------------------------------------*/
MAT::GrowthMandel::GrowthMandel(MAT::PAR::Growth* params)
  : Growth(params),
    mandel_(Teuchos::null),
    paramsMandel_(params)
{
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthMandel::Evaluate( const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain,
    Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>* stress,
    LINALG::Matrix<6, 6>* cmat,
    const int eleGID)
{
  // get gauss point number
  const int gp = params.get<int>("gp", -1);
  if (gp == -1)
    dserror("no Gauss point number provided in material");

  double time = params.get<double>("total time", -1.0);
  if( abs(time+1.0) < 1e-14 ) dserror("no time step or no total time given for growth material!");
  std::string action = params.get<std::string>("action", "none");
  bool output = false;
  if (action == "calc_struct_stress")
    output = true;

  const double eps = 1.0e-14;
  MAT::PAR::Growth* growth_params = Parameter();
  const double endtime = growth_params->endtime_;
  const double starttime = growth_params->starttime_;
  // when stress output is calculated the final parameters already exist
  // we should not do another local Newton iteration, which uses eventually a wrong thetaold
  if (output)
    time = endtime + 1.0;

  if (time > starttime + eps && time <= endtime + eps)
  {
    double theta = ThetaOld()->at(gp);
    //--------------------------------------------------------------------------------------
    // evaluation of the volumetric growth factor and its derivative wrt cauchy-green
    //--------------------------------------------------------------------------------------
    LINALG::Matrix<6,1> dthetadC(true);
    EvaluateGrowth(&theta,&dthetadC,defgrd,glstrain,params,eleGID);

    //--------------------------------------------------------------------------------------
    // call material law with elastic part of defgr and elastic part of glstrain
    //--------------------------------------------------------------------------------------
    // build identity tensor I
    LINALG::Matrix<NUM_STRESS_3D, 1> Id(true);
    for (int i = 0; i < 3; i++)
      Id(i) = 1.0;

    // right Cauchy-Green Tensor  C = 2 * E + I
    LINALG::Matrix<NUM_STRESS_3D, 1> C(*glstrain);
    C.Scale(2.0);
    C += Id;

    // elastic right Cauchy-Green Tensor Cdach = F_g^-T C F_g^-1
    LINALG::Matrix<NUM_STRESS_3D, 1> Cdach(C);
    Cdach.Scale(1.0 / theta / theta);
    LINALG::Matrix<3, 3> defgrddach(*defgrd);
    defgrddach.Scale(1.0 / theta);
    // elastic Green Lagrange strain
    LINALG::Matrix<NUM_STRESS_3D, 1> glstraindach(Cdach);
    glstraindach -= Id;
    glstraindach.Scale(0.5);

    //modify the parameter list to be passed to the elastic material
     Teuchos::ParameterList paramselast(params);
    paramselast.remove("matparderiv",false);

    // elastic 2 PK stress and constitutive matrix
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatelastic(true);
    LINALG::Matrix<NUM_STRESS_3D, 1> Sdach(true);
    EvaluateElastic(&defgrddach,
                    &glstraindach,
                    &Sdach,
                    &cmatelastic,
                    paramselast,
                    eleGID);


    // 2PK stress S = F_g^-1 Sdach F_g^-T
    LINALG::Matrix<NUM_STRESS_3D, 1> S(Sdach);
    S.Scale(1.0 / theta / theta);
    *stress = S;

    // constitutive matrix including growth
    cmatelastic.Scale(1.0 / theta / theta / theta / theta);

    for (int i = 0; i < 6; i++)
    {
      double cmatelasCi =   cmatelastic(i, 0) * C(0) + cmatelastic(i, 1) * C(1)
                                  + cmatelastic(i, 2) * C(2) + cmatelastic(i, 3) * C(3)
                                  + cmatelastic(i, 4) * C(4) + cmatelastic(i, 5) * C(5);

      for (int j = 0; j < 6; j++)
        (*cmat)(i, j) =  cmatelastic(i, j) - 2.0/theta *(2.0 * S(i) + cmatelasCi) * dthetadC(j);
    }

    int deriv = params.get<int>("matparderiv",-1);
    if (deriv != -1)
    {
      LINALG::Matrix<NUM_STRESS_3D, 1> cmatelasC(true);
      for (int i = 0; i < NUM_STRESS_3D; i++)
       {
         cmatelasC(i,0) = cmatelastic(i, 0) * C(0) + cmatelastic(i, 1) * C(1)
                       + cmatelastic(i, 2) * C(2) + cmatelastic(i, 3) * C(3)
                       + cmatelastic(i, 4) * C(4) + cmatelastic(i, 5) * C(5);
       }

      stress->Update(1.0,cmatelasC,2.0);
      double dthetadp;
      Parameter()->growthlaw_->EvaluatePDeriv(&dthetadp,ThetaOld()->at(gp),Matelastic(),defgrd, glstrain,params,eleGID);
      stress->Scale(-dthetadp/theta);

    }

    // store theta
    theta_->at(gp) = theta;

    // trace of elastic Mandel stress Mdach = Cdach Sdach
    double mandel =   Cdach(0) * Sdach(0) + Cdach(1) * Sdach(1)
                    + Cdach(2) * Sdach(2) + Cdach(3) * Sdach(3) + Cdach(4) * Sdach(4)
                    + Cdach(5) * Sdach(5);
    mandel_->at(gp) = mandel;

  }
  else if (time > endtime + eps)
  { // turn off growth or calculate stresses for output
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatelastic(true);
    LINALG::Matrix<NUM_STRESS_3D, 1> Sdach(true);
    double theta = theta_->at(gp);

    //--------------------------------------------------------------------------------------
    // build identity tensor I
    LINALG::Matrix<NUM_STRESS_3D, 1> Id(true);
    for (int i = 0; i < 3; i++)
      Id(i) = 1.0;

    // right Cauchy-Green Tensor  C = 2 * E + I
    LINALG::Matrix<NUM_STRESS_3D, 1> C(*glstrain);
    C.Scale(2.0);
    C += Id;

    // elastic right Cauchy-Green Tensor Cdach = F_g^-T C F_g^-1
    LINALG::Matrix<NUM_STRESS_3D, 1> Cdach(C);
    Cdach.Scale(1.0 / theta / theta);
    LINALG::Matrix<3, 3> defgrddach(*defgrd);
    defgrddach.Scale(1.0 / theta);
    // elastic Green Lagrange strain
    LINALG::Matrix<NUM_STRESS_3D, 1> glstraindach(Cdach);
    glstraindach -= Id;
    glstraindach.Scale(0.5);
    // elastic 2 PK stress and constitutive matrix
    EvaluateElastic(&defgrddach,
                    &glstraindach,
                    &Sdach,
                    &cmatelastic,
                    params,
                    eleGID);

    // 2PK stress S = F_g^-1 Sdach F_g^-T
    LINALG::Matrix<NUM_STRESS_3D, 1> S(Sdach);
    S.Scale(1.0 / theta / theta);
    *stress = S;

    // constitutive matrix including growth
    cmatelastic.Scale(1.0 / theta / theta / theta / theta);
    *cmat = cmatelastic;

    // trace of elastic Mandel stress Mdach = Cdach Sdach
    double mandel =   Cdach(0) * Sdach(0) + Cdach(1) * Sdach(1)
                    + Cdach(2) * Sdach(2) + Cdach(3) * Sdach(3) + Cdach(4) * Sdach(4)
                    + Cdach(5) * Sdach(5);
    mandel_->at(gp) = mandel;

  }
  else
  {
    EvaluateElastic(defgrd, glstrain, stress, cmat, params, eleGID);
    // build identity tensor I
    LINALG::Matrix<NUM_STRESS_3D, 1> Id(true);
    for (int i = 0; i < 3; i++)
      Id(i) = 1.0;
    // right Cauchy-Green Tensor  C = 2 * E + I
    LINALG::Matrix<NUM_STRESS_3D, 1> C(*glstrain);
    C.Scale(2.0);
    C += Id;
    LINALG::Matrix<NUM_STRESS_3D, 1> S(true);
    S = *stress;
    mandel_->at(gp) =   C(0) * S(0) + C(1) * S(1) + C(2) * S(2) + C(3) * S(3)
                      + C(4) * S(4) + C(5) * S(5);
  }
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthMandel::EvaluateGrowth(double* theta,
                      LINALG::Matrix<6,1>* dthetadC,
                      const LINALG::Matrix<3,3>* defgrd,
                      const LINALG::Matrix<6,1>* glstrain,
                      Teuchos::ParameterList& params,
                      const int eleGID)
{
  // get gauss point number
  const int gp = params.get<int>("gp", -1);
  if (gp == -1)
    dserror("no Gauss point number provided in material");

  double thetaold = ThetaOld()->at(gp);

  MAT::Growth* matgrowth = this;
  Parameter()->growthlaw_->Evaluate(theta,thetaold,dthetadC,*matgrowth,defgrd,glstrain,params,eleGID);
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthMandel::ResetAll(const int numgp)
{
  for (int j=0; j<numgp; ++j)
    mandel_->at(j) = 0.0;

  MAT::Growth::ResetAll(numgp);
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthMandel::EvaluateNonLinMass( const LINALG::Matrix<3, 3>* defgrd,
                            const LINALG::Matrix<6, 1>* glstrain,
                            Teuchos::ParameterList& params,
                            LINALG::Matrix<NUM_STRESS_3D,1>* linmass_disp,
                            LINALG::Matrix<NUM_STRESS_3D,1>* linmass_vel,
                            const int eleGID)
{
  double eps = 1.0e-12;
  double starttime = Parameter()->starttime_;
  double endtime = Parameter()->endtime_;
  double time = params.get<double>("total time", -1.0);

  if (time > starttime + eps and time <= endtime + eps)
  {
    // get gauss point number
    const int gp = params.get<int>("gp", -1);
    if (gp == -1)
      dserror("no Gauss point number provided in material");

    double theta = theta_->at(gp);
    double thetaold = ThetaOld()->at(gp);

    MAT::Growth* matgrowth = this;

    Parameter()->growthlaw_->Evaluate(&theta,thetaold,linmass_disp,*matgrowth,defgrd,glstrain,params,eleGID);
    linmass_disp->Scale(3.0*theta*theta*Matelastic()->Density());

    linmass_vel->Clear();
  }
  else
  {
    //no growth. set to zero
    linmass_disp->Clear();
    linmass_vel->Clear();
  }
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthMandel::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matid
  int matid = -1;

  MAT::PAR::Growth* params=Parameter();
  if (params != NULL) matid = params->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);

  int numgp=0;
  if (isinit_)
  {
    numgp = mandel_->size();   // size is number of gausspoints
  }
  AddtoPack(data,numgp);
  // Pack internal variables
  for (int gp = 0; gp < numgp; ++gp)
  {
    AddtoPack(data,mandel_->at(gp));
  }

  // Pack base class material
  Growth::Pack(data);

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthMandel::Unpack(const std::vector<char>& data)
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

  paramsMandel_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        paramsMandel_ = dynamic_cast<MAT::PAR::Growth*>(mat);
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
  mandel_ = Teuchos::rcp(new std::vector<double> (numgp));
  for (int gp = 0; gp < numgp; ++gp)
  {
    double a;
    ExtractfromPack(position,data,a);
    mandel_->at(gp) = a;
  }

  // extract base class material
  std::vector<char> basedata(0);
  Growth::ExtractfromPack(position,data,basedata);
  Growth::Unpack(basedata);


  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthMandel::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  mandel_ = Teuchos::rcp(new std::vector<double> (numgp));
  for (int j=0; j<numgp; ++j)
    mandel_->at(j) = 0.0;

  //setup base class
  Growth::Setup(numgp, linedef);
  return;
}
