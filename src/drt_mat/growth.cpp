/*----------------------------------------------------------------------*/
/*!
\brief This file contains routines for integration point based isotropic and anisotropic volumetric
growth laws.

\level 2

<pre>
   \maintainer Amadeus Gebauer
</pre>
 *----------------------------------------------------------------------*/


#include "growth.H"

#include "growth_law.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_utils_factory.H"


/*----------------------------------------------------------------------------*/
MAT::PAR::Growth::Growth(Teuchos::RCP<MAT::PAR::Material> matdata)
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
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(idgrowthlaw_);

  switch (curmat->Type())
  {
    case INPAR::MAT::m_growth_aniso_strain:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::GrowthLawAnisoStrain(curmat));
      MAT::PAR::GrowthLawAnisoStrain* params =
          static_cast<MAT::PAR::GrowthLawAnisoStrain*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case INPAR::MAT::m_growth_aniso_stress:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::GrowthLawAnisoStress(curmat));
      MAT::PAR::GrowthLawAnisoStress* params =
          static_cast<MAT::PAR::GrowthLawAnisoStress*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case INPAR::MAT::m_growth_aniso_strain_const_trig:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::GrowthLawAnisoStrainConstTrig(curmat));
      MAT::PAR::GrowthLawAnisoStrainConstTrig* params =
          static_cast<MAT::PAR::GrowthLawAnisoStrainConstTrig*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case INPAR::MAT::m_growth_aniso_stress_const_trig:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::GrowthLawAnisoStressConstTrig(curmat));
      MAT::PAR::GrowthLawAnisoStressConstTrig* params =
          static_cast<MAT::PAR::GrowthLawAnisoStressConstTrig*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case INPAR::MAT::m_growth_iso_stress:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::GrowthLawIsoStress(curmat));
      MAT::PAR::GrowthLawIsoStress* params =
          static_cast<MAT::PAR::GrowthLawIsoStress*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case INPAR::MAT::m_growth_ac:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::GrowthLawAC(curmat));
      MAT::PAR::GrowthLawAC* params = static_cast<MAT::PAR::GrowthLawAC*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case INPAR::MAT::m_growth_ac_radial:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::GrowthLawACRadial(curmat));
      MAT::PAR::GrowthLawACRadial* params =
          static_cast<MAT::PAR::GrowthLawACRadial*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case INPAR::MAT::m_growth_ac_radial_refconc:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::GrowthLawACRadialRefConc(curmat));
      MAT::PAR::GrowthLawACRadialRefConc* params =
          static_cast<MAT::PAR::GrowthLawACRadialRefConc*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case INPAR::MAT::m_growth_const:
    {
      if (curmat->Parameter() == NULL) curmat->SetParameter(new MAT::PAR::GrowthLawConst(curmat));
      MAT::PAR::GrowthLawConst* params =
          static_cast<MAT::PAR::GrowthLawConst*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    default:
      dserror("unknown material type %d", curmat->Type());
      break;
  }

  if (starttime_ > endtime_)
    dserror("WTF! It is not reasonable to have a starttime that is larger than the endtime!");
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::Growth::CreateMaterial()
{
  Teuchos::RCP<MAT::Material> mat;

  switch (growthlaw_->MaterialType())
  {
    case INPAR::MAT::m_growth_aniso_strain:
    case INPAR::MAT::m_growth_aniso_stress:
    case INPAR::MAT::m_growth_aniso_strain_const_trig:
    case INPAR::MAT::m_growth_aniso_stress_const_trig:
    case INPAR::MAT::m_growth_iso_stress:
    case INPAR::MAT::m_growth_ac:
    case INPAR::MAT::m_growth_ac_radial:
    case INPAR::MAT::m_growth_ac_radial_refconc:
    case INPAR::MAT::m_growth_const:
      mat = Teuchos::rcp(new MAT::GrowthVolumetric(this));
      break;
    default:
      dserror(
          "The growth law you have chosen is not valid for the standard volumetric growth "
          "material");
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
  DRT::PackBuffer::SizeMarker sm(data);
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
    numgp = 0;  // not initialized -> nothing to pack
  }
  else
  {
    numgp = theta_->size();  // size is number of Gauss points
  }
  AddtoPack(data, numgp);
  // Pack internal variables
  for (int gp = 0; gp < numgp; ++gp)
  {
    AddtoPack(data, thetaold_->at(gp));
    AddtoPack(data, theta_->at(gp));
  }

  AddtoPack(data, histdata_);

  // Pack data of elastic material
  if (matelastic_ != Teuchos::null)
  {
    matelastic_->Pack(data);
  }

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::Growth::Unpack(const std::vector<char>& data)
{
  isinit_ = true;
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);

  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = dynamic_cast<MAT::PAR::Growth*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  int numgp;
  ExtractfromPack(position, data, numgp);
  if (numgp == 0)
  {  // no history data to unpack
    isinit_ = false;
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d", data.size(), position);
    return;
  }

  // unpack growth internal variables
  theta_ = Teuchos::rcp(new std::vector<double>(numgp));
  thetaold_ = Teuchos::rcp(new std::vector<double>(numgp));
  for (int gp = 0; gp < numgp; ++gp)
  {
    double a;
    ExtractfromPack(position, data, a);
    thetaold_->at(gp) = a;
    ExtractfromPack(position, data, a);
    theta_->at(gp) = a;
  }

  std::vector<char> tmp(0);
  ExtractfromPack(position, data, tmp);
  histdata_.Unpack(tmp);

  // Unpack data of elastic material (these lines are copied from drt_element.cpp)
  std::vector<char> dataelastic;
  ExtractfromPack(position, data, dataelastic);
  if (dataelastic.size() > 0)
  {
    DRT::ParObject* o = DRT::UTILS::Factory(dataelastic);  // Unpack is done here
    MAT::So3Material* matel = dynamic_cast<MAT::So3Material*>(o);
    if (matel == NULL) dserror("failed to unpack elastic material");
    matelastic_ = Teuchos::rcp(matel);
  }
  else
    matelastic_ = Teuchos::null;

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::Growth::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  if (isinit_)
    dserror("This function should just be called if the material is not yet initialized.");

  theta_ = Teuchos::rcp(new std::vector<double>(numgp));
  thetaold_ = Teuchos::rcp(new std::vector<double>(numgp));
  for (int j = 0; j < numgp; ++j)
  {
    theta_->at(j) = 1.0;
    thetaold_->at(j) = 1.0;
  }

  // Setup of elastic material
  matelastic_ =
      Teuchos::rcp_dynamic_cast<MAT::So3Material>(MAT::Material::Factory(params_->idmatelastic_));
  matelastic_->Setup(numgp, linedef);

  // Setup of history container with empty data
  std::map<int, std::vector<double>> data;
  histdata_.Add("thetaold", data);

  isinit_ = true;
  return;
}

/*----------------------------------------------------------------------------*/
void MAT::Growth::ResetAll(const int numgp)
{
  for (int j = 0; j < numgp; ++j)
  {
    thetaold_->at(j) = 1.0;
    theta_->at(j) = 1.0;
  }

  // overwrite history with empty data
  std::map<int, std::vector<double>> data;
  histdata_.Add("thetaold", data);

  matelastic_->ResetAll(numgp);
}

/*----------------------------------------------------------------------------*/
void MAT::Growth::Update()
{
  const int numgp = theta_->size();

  for (int i = 0; i < numgp; i++)
  {
    thetaold_->at(i) = theta_->at(i);
  }

  matelastic_->Update();
}

/*----------------------------------------------------------------------------*/
void MAT::Growth::ResetStep()
{
  const int numgp = theta_->size();

  for (int i = 0; i < numgp; i++)
  {
    theta_->at(i) = thetaold_->at(i);
  }

  matelastic_->ResetStep();
}

void MAT::Growth::StoreHistory(int timestep)
{
  std::map<int, std::vector<double>>* access;
  access = histdata_.GetMutable<std::map<int, std::vector<double>>>("thetaold");
  (*access)[timestep] = *thetaold_;
}

void MAT::Growth::SetHistory(int timestep)
{
  const std::map<int, std::vector<double>>* read;
  read = histdata_.Get<std::map<int, std::vector<double>>>("thetaold");

  if (read->find(timestep) != read->end())
    *thetaold_ = read->at(timestep);
  else
    dserror("there is no data to reset for step %d", timestep);
}

/*----------------------------------------------------------------------------*/
void MAT::Growth::EvaluateElastic(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
    Teuchos::ParameterList& params, const int eleGID)
{
  Matelastic()->Evaluate(defgrd, glstrain, params, stress, cmat, eleGID);
}


double MAT::Growth::Density(int gp) const
{
  const double density_elast = matelastic_->Density();
  const double theta_gp = (*theta_)[gp];

  const double density_scale = Parameter()->growthlaw_->DensityScale(theta_gp);

  return density_elast * density_scale;
}

/*----------------------------------------------------------------------*
 | returns whether material density is constant (public)  schmidt 11/17 |
 *----------------------------------------------------------------------*/
bool MAT::Growth::VaryingDensity() const { return Parameter()->growthlaw_->VaryingDensity(); }


/*----------------------------------------------------------------------------*/
void MAT::GrowthVolumetric::VisNames(std::map<std::string, int>& names)
{
  std::string name = "theta";
  names[name] = 1;

  switch (Parameter()->growthlaw_->MaterialType())
  {
    case INPAR::MAT::m_growth_aniso_stress:
    case INPAR::MAT::m_growth_aniso_stress_const_trig:
    case INPAR::MAT::m_growth_iso_stress:
    {
      name = "tr_mandel_e";
      names[name] = 1;
    }
    break;
    case INPAR::MAT::m_growth_aniso_strain:
    case INPAR::MAT::m_growth_aniso_strain_const_trig:
    {
      name = "lambda_fib_e";
      names[name] = 1;
    }
    break;
    default:
      break;
  }

  Matelastic()->VisNames(names);
}

/*----------------------------------------------------------------------------*/
bool MAT::GrowthVolumetric::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  if (name == "theta")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double temp = 0.0;
    for (int gp = 0; gp < numgp; gp++) temp += theta_->at(gp);
    data[0] = temp / numgp;
  }
  else if (name == "tr_mandel_e")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double temp = 0.0;
    for (int gp = 0; gp < numgp; gp++) temp += tr_mandel_e_->at(gp);
    data[0] = temp / numgp;
  }
  else if (name == "lambda_fib_e")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double temp = 0.0;
    for (int gp = 0; gp < numgp; gp++) temp += lambda_fib_e_->at(gp);
    data[0] = temp / numgp;
  }
  else
  {
    return Matelastic()->VisData(name, data, numgp, eleID);
  }
  return true;
}


MAT::GrowthVolumetricType MAT::GrowthVolumetricType::instance_;

/*----------------------------------------------------------------------------*/
DRT::ParObject* MAT::GrowthVolumetricType::Create(const std::vector<char>& data)
{
  MAT::GrowthVolumetric* grow = new MAT::GrowthVolumetric();
  grow->Unpack(data);
  return grow;
}

/*----------------------------------------------------------------------------*/
MAT::GrowthVolumetric::GrowthVolumetric()
    : Growth(),
      tr_mandel_e_(Teuchos::null),
      lambda_fib_e_(Teuchos::null),
      growthtrig_const_(0.0),
      paramsVolumetric_(NULL),
      refdir_(true),
      curdir_(Teuchos::null),
      curdir_for_update_(Teuchos::null),
      F_g_hist_(Teuchos::null)
{
}

/*----------------------------------------------------------------------------*/
MAT::GrowthVolumetric::GrowthVolumetric(MAT::PAR::Growth* params)
    : Growth(params),
      tr_mandel_e_(Teuchos::null),
      lambda_fib_e_(Teuchos::null),
      growthtrig_const_(0.0),
      paramsVolumetric_(params),
      refdir_(true),
      curdir_(Teuchos::null),
      curdir_for_update_(Teuchos::null),
      F_g_hist_(Teuchos::null)
{
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthVolumetric::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, const int eleGID)
{
  // get gauss point number
  const int gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("no Gauss point number provided in material");

  double time = params.get<double>("total time", -1.0);
  if (abs(time + 1.0) < 1e-14) dserror("no time step or no total time given for growth material!");
  std::string action = params.get<std::string>("action", "none");
  bool output = false;
  if (action == "calc_struct_stress") output = true;

  const double eps = 1.0e-14;
  MAT::PAR::Growth* growth_params = Parameter();
  const double endtime = growth_params->endtime_;
  const double starttime = growth_params->starttime_;

  // growth is allowed in here
  if ((((time > starttime + eps) and (time <= endtime + eps)) or
          ((starttime < 0.0) and (endtime < 0.0))) and
      !output)
  {
    // not nice but currently the only way we may do that....
    switch (Parameter()->growthlaw_->MaterialType())
    {
      case INPAR::MAT::m_growth_ac_radial:
      case INPAR::MAT::m_growth_ac_radial_refconc:
      {
        // directional stuff......
        // push-forward of refdir
        LINALG::Matrix<3, 3> defgrdinv(true);
        defgrdinv.Invert(*defgrd);
        LINALG::Matrix<3, 1> curdir_for_update(true);
        curdir_for_update.MultiplyTN(defgrd->Determinant(), defgrdinv, refdir_);
        // scale n to length of one
        curdir_for_update.Scale(1.0 / curdir_for_update.Norm2());
        // save for time update
        curdir_for_update_.at(gp) = curdir_for_update;
        //  curdir_.at(gp)=curdir_for_update;
      }
      break;
      default:
        break;
    }

    double theta = ThetaOld()->at(gp);
    //--------------------------------------------------------------------------------------
    // evaluation of the volumetric growth factor and its derivative wrt cauchy-green
    //--------------------------------------------------------------------------------------
    LINALG::Matrix<6, 1> dthetadCvec(true);
    EvaluateGrowth(&theta, &dthetadCvec, defgrd, glstrain, params, eleGID);

    // modify the parameter list to be passed to the elastic material
    Teuchos::ParameterList paramselast(params);
    paramselast.remove("matparderiv", false);

    LINALG::Matrix<6, 1> S(true);
    LINALG::Matrix<6, 6> cmatdach(true);

    GetSAndCmatdach(theta, defgrd, &S, &cmatdach, paramselast, eleGID);

    *stress = S;

    // calculate growth part F_g of the deformation gradient F
    LINALG::Matrix<3, 3> F_g(true);

    Parameter()->growthlaw_->CalcFg(
        theta, ThetaOld()->at(gp), gp, defgrd, refdir_, curdir_, F_g_hist_, F_g);

    // calculate F_g^(-1)
    LINALG::Matrix<3, 3> F_ginv(true);
    F_ginv.Invert(F_g);

    // constitutive matrix including growth cmat = F_g^-1 F_g^-1 cmatdach F_g^-T F_g^-T
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatelast(true);
    cmatelast = PullBack4Tensor(F_ginv, cmatdach);

    //--------------------------------------------------------------------------------------
    // call material law with elastic part of defgr and elastic part of glstrain
    //--------------------------------------------------------------------------------------
    // build identity tensor I
    LINALG::Matrix<NUM_STRESS_3D, 1> Id(true);
    for (int i = 0; i < 3; i++) Id(i) = 1.0;

    // right Cauchy-Green Tensor  C = 2 * E + I
    LINALG::Matrix<NUM_STRESS_3D, 1> C(*glstrain);
    C.Scale(2.0);
    C += Id;

    // NOTE: we do this by a FD approximation, which is really cheap here due to the fact
    // that theta is a scalar value (hence only one more material evaluation is needed!)

    const double espilon = 1.0e-8;

    LINALG::Matrix<6, 1> SEps(true);
    LINALG::Matrix<6, 6> cmatdachEps(true);

    GetSAndCmatdach(theta + espilon, defgrd, &SEps, &cmatdachEps, params, eleGID);

    //--------------------------------------------------------------------------------------
    // calculate \frac{d S}{d C} = \frac{\partial S}{\partial C} +
    //    + 2* \left( \frac{\partial S}{\partial \theta} \otimes \frac{\partial \theta}{\partial C}
    //    \right)
    //--------------------------------------------------------------------------------------
    for (int i = 0; i < 6; i++)
    {
      for (int j = 0; j < 6; j++)
      {
        (*cmat)(i, j) = cmatelast(i, j) + 2.0 * (SEps(i) - S(i)) / espilon * dthetadCvec(j);
      }
    }


    // stuff needed for const growth law in combination with parameter estimation
    int deriv = params.get<int>("matparderiv", -1);
    if (deriv != -1)
    {
      LINALG::Matrix<NUM_STRESS_3D, 1> cmatelasC(true);
      for (int i = 0; i < NUM_STRESS_3D; i++)
      {
        cmatelasC(i, 0) = cmatelast(i, 0) * C(0) + cmatelast(i, 1) * C(1) + cmatelast(i, 2) * C(2) +
                          cmatelast(i, 3) * C(3) + cmatelast(i, 4) * C(4) + cmatelast(i, 5) * C(5);
      }

      stress->Update(1.0, cmatelasC, 2.0);
      double dthetadp;
      Parameter()->growthlaw_->EvaluatePDeriv(
          &dthetadp, ThetaOld_atgp(gp), Matelastic(), defgrd, glstrain, params, eleGID);
      stress->Scale(-dthetadp / theta);

    }  // END stuff needed for const growth law in combination with parameter estimation

    // store theta
    theta_->at(gp) = theta;
  }

  // when stress output is calculated the final parameters already exist
  // we should not do another local Newton iteration, which uses eventually a wrong thetaold
  else if ((time > endtime + eps) or output)
  {  // turn off growth or calculate stresses for output
    double theta = theta_->at(gp);

    LINALG::Matrix<6, 1> Svec(true);
    LINALG::Matrix<6, 6> cmatdach(true);

    GetSAndCmatdach(theta, defgrd, &Svec, &cmatdach, params, eleGID);

    *stress = Svec;

    // calculate growth part F_g of the deformation gradient F
    LINALG::Matrix<3, 3> F_g(true);

    Parameter()->growthlaw_->CalcFg(
        theta, ThetaOld_atgp(gp), gp, defgrd, refdir_, curdir_, F_g_hist_, F_g);

    // calculate F_g^(-1)
    LINALG::Matrix<3, 3> F_ginv(true);
    F_ginv.Invert(F_g);

    // constitutive matrix including growth cmat = F_g^-1 F_g^-1 cmatdach F_g^-T F_g^-T
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatelast(true);
    cmatelast = PullBack4Tensor(F_ginv, cmatdach);

    *cmat = cmatelast;
  }
  else
  {
    EvaluateElastic(defgrd, glstrain, stress, cmat, params, eleGID);
    // build identity tensor I
    LINALG::Matrix<NUM_STRESS_3D, 1> Id(true);
    for (int i = 0; i < 3; i++) Id(i) = 1.0;
    // right Cauchy-Green Tensor  C = 2 * E + I
    LINALG::Matrix<NUM_STRESS_3D, 1> Cvec(*glstrain);
    Cvec.Scale(2.0);
    Cvec += Id;
    LINALG::Matrix<NUM_STRESS_3D, 1> Svec(true);
    Svec = *stress;

    tr_mandel_e_->at(gp) = Cvec(0) * Svec(0) + Cvec(1) * Svec(1) + Cvec(2) * Svec(2) +
                           Cvec(3) * Svec(3) + Cvec(4) * Svec(4) + Cvec(5) * Svec(5);

    // elastic fiber stretch
    LINALG::Matrix<3, 3> C(true);
    VectorToMatrix(C, Cvec, MAT::voigt_strain);

    LINALG::Matrix<3, 1> CDir(true);
    CDir.MultiplyNN(1.0, C, refdir_);
    lambda_fib_e_->at(gp) =
        sqrt(CDir(0) * refdir_(0) + CDir(1) * refdir_(1) + CDir(2) * refdir_(2));
  }
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthVolumetric::EvaluateGrowth(double* theta, LINALG::Matrix<6, 1>* dthetadC,
    const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
    Teuchos::ParameterList& params, const int eleGID)
{
  // get gauss point number
  const int gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("No Gauss point number provided in material.");

  double thetaold = ThetaOld_atgp(gp);

  MAT::Growth* matgrowth = this;
  Parameter()->growthlaw_->Evaluate(theta, thetaold, dthetadC, *matgrowth, defgrd, glstrain,
      refdir_, curdir_, F_g_hist_, growthtrig_const_, params, eleGID);
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthVolumetric::ResetAll(const int numgp)
{
  for (int j = 0; j < numgp; ++j)
  {
    tr_mandel_e_->at(j) = 0.0;
    lambda_fib_e_->at(j) = 0.0;
  }

  MAT::Growth::ResetAll(numgp);
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthVolumetric::EvaluateNonLinMass(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<NUM_STRESS_3D, 1>* linmass_disp, LINALG::Matrix<NUM_STRESS_3D, 1>* linmass_vel,
    const int eleGID)
{
  const double eps = 1.0e-14;
  const double starttime = Parameter()->starttime_;
  const double endtime = Parameter()->endtime_;
  const double time = params.get<double>("total time", -1.0);

  if (((time > starttime + eps) and (time <= endtime + eps)) or
      ((starttime < 0.0) and (endtime < 0.0)))
  {
    // get gauss point number
    const int gp = params.get<int>("gp", -1);
    if (gp == -1) dserror("no Gauss point number provided in material");

    double theta = theta_->at(gp);
    double thetaold = ThetaOld()->at(gp);

    MAT::Growth* matgrowth = this;
    Parameter()->growthlaw_->Evaluate(&theta, thetaold, linmass_disp, *matgrowth, defgrd, glstrain,
        refdir_, curdir_, F_g_hist_, growthtrig_const_, params, eleGID);

    const double density_deriv_scale = Parameter()->growthlaw_->DensityDerivScale(theta);
    linmass_disp->Scale(density_deriv_scale * Matelastic()->Density());

    linmass_vel->Clear();
  }
  else
  {
    // no growth, set to zero
    linmass_disp->Clear();
    linmass_vel->Clear();
  }
}



///*----------------------------------------------------------------------*
// | calculate stresses and elastic material tangent                      |
// | (both in Voigt notation)                                   Thon 01/16|
// *----------------------------------------------------------------------*/
void MAT::GrowthVolumetric::GetSAndCmatdach(const double theta, const LINALG::Matrix<3, 3>* defgrd,
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmatdach, Teuchos::ParameterList& params,
    const int eleGID)
{
  // get gauss point number
  const int gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("no Gauss point number provided in material");

  // calculate growth part F_g of the deformation gradient F
  LINALG::Matrix<3, 3> F_g(true);
  Parameter()->growthlaw_->CalcFg(
      theta, ThetaOld_atgp(gp), gp, defgrd, refdir_, curdir_, F_g_hist_, F_g);

  // calculate F_g^(-1)
  LINALG::Matrix<3, 3> F_ginv(true);
  F_ginv.Invert(F_g);

  // elastic deformation gradient F_e = F * F_g^(-1)
  LINALG::Matrix<3, 3> defgrddach(true);   //*defgrd);
  defgrddach.MultiplyNN(*defgrd, F_ginv);  // Scale(1.0 / theta);

  // elastic right Cauchy-Green Tensor Cdach = F_e^T * F_e (= F_g^-T C F_g^-1)
  LINALG::Matrix<3, 3> Cdach(true);
  Cdach.MultiplyTN(defgrddach, defgrddach);

  // transform Cdach into a vector
  LINALG::Matrix<6, 1> Cdachvec(true);
  MatrixToVector(Cdach, Cdachvec, MAT::voigt_strain);

  //--------------------------------------------------------------------------------------
  // call material law with elastic part of defgr and elastic part of glstrain
  //--------------------------------------------------------------------------------------
  // build identity tensor I
  LINALG::Matrix<6, 1> Id(true);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;

  LINALG::Matrix<6, 1> glstraindachvec(Cdachvec);
  glstraindachvec -= Id;
  glstraindachvec.Scale(0.5);

  LINALG::Matrix<6, 1> Sdachvec(true);
  // elastic 2 PK stress and constitutive matrix
  Matelastic()->Evaluate(&defgrddach, &glstraindachvec, params, &Sdachvec, cmatdach, eleGID);

  // calculate stress
  // 2PK stress S = F_g^-1 Sdach F_g^-T
  LINALG::Matrix<3, 3> Sdach(true);
  VectorToMatrix(Sdach, Sdachvec, MAT::voigt_stress);

  LINALG::Matrix<3, 3> tmp(true);
  tmp.MultiplyNT(Sdach, F_ginv);
  LINALG::Matrix<3, 3> S(true);
  S.MultiplyNN(F_ginv, tmp);

  MatrixToVector(S, *stress, MAT::voigt_stress);

  // trace of elastic Mandel stress Mdach = Cdach Sdach
  tr_mandel_e_->at(gp) = Cdachvec(0) * Sdachvec(0) + Cdachvec(1) * Sdachvec(1) +
                         Cdachvec(2) * Sdachvec(2) + Cdachvec(3) * Sdachvec(3) +
                         Cdachvec(4) * Sdachvec(4) + Cdachvec(5) * Sdachvec(5);
  ;

  // elastic fiber stretch lambda = \sqrt(f_0 \cdot Cdach f_0)
  LINALG::Matrix<3, 1> CdachDir(true);
  CdachDir.MultiplyNN(1.0, Cdach, refdir_);
  lambda_fib_e_->at(gp) =
      sqrt(CdachDir(0) * refdir_(0) + CdachDir(1) * refdir_(1) + CdachDir(2) * refdir_(2));
}



/*----------------------------------------------------------------------*
 | transform vector in voigt notation into symmetric                    |
 | 2Tensor (in matrix notation)                              Thon 01/15 |
 *----------------------------------------------------------------------*/
void MAT::GrowthVolumetric::VectorToMatrix(
    LINALG::Matrix<3, 3>& Matrix, const LINALG::Matrix<6, 1>& Vector, const MAT::VoigtType Type)
{
  double alpha;

  switch (Type)
  {
    case MAT::voigt_stress:
      alpha = 1.0;
      break;
    case MAT::voigt_strain:
      alpha = 0.5;
      break;
    default:
      dserror("No supported VoigtType!");
      alpha = 1.0;
      break;
  }

  Matrix(0, 0) = Vector(0);
  Matrix(0, 1) = alpha * Vector(3);
  Matrix(0, 2) = alpha * Vector(5);
  Matrix(1, 0) = Matrix(0, 1);  // alpha*Vector(3);
  Matrix(1, 1) = Vector(1);
  Matrix(1, 2) = alpha * Vector(4);
  Matrix(2, 0) = Matrix(0, 2);  // alpha*Vector(5);
  Matrix(2, 1) = Matrix(1, 2);  // alpha*Vector(4);
  Matrix(2, 2) = Vector(2);
}

/*----------------------------------------------------------------------*
 | ///transform symmetric 2Tensor(in matrix notation) into voigt        |
 | notation ( e.g. vector notation)                          Thon 01/15 |
 *----------------------------------------------------------------------*/
void MAT::GrowthVolumetric::MatrixToVector(
    const LINALG::Matrix<3, 3>& Matrix, LINALG::Matrix<6, 1>& Vector, const MAT::VoigtType Type)
{
  double alpha;
  switch (Type)
  {
    case MAT::voigt_stress:
      alpha = 1.0;
      break;
    case MAT::voigt_strain:
      alpha = 2.0;
      break;
    default:
      dserror("No supported VoigtType!");
      alpha = 1.0;
      break;
  }
  Vector(0) = Matrix(0, 0);
  Vector(1) = Matrix(1, 1);
  Vector(2) = Matrix(2, 2);
  Vector(3) = alpha * Matrix(0, 1);
  Vector(4) = alpha * Matrix(1, 2);
  Vector(5) = alpha * Matrix(0, 2);
}

/*----------------------------------------------------------------------*
 | pull back of a symmertic elastic 4th order tensor (in matrix/voigt   |
 | notation) via the 2th order deformation gradient (also in matrix     |
 | notation)                                                  thon 01/15|
 *----------------------------------------------------------------------*/
LINALG::Matrix<6, 6> MAT::GrowthVolumetric::PullBack4Tensor(
    const LINALG::Matrix<3, 3>& defgr, const LINALG::Matrix<6, 6>& Cmat)
{
  double CMAT[3][3][3][3] = {{{{0.0}}}};
  Setup4Tensor(CMAT, Cmat);
  //  PrintFourTensor(CMAT);

  // This would be the long way....
  //
  //  double tmp[3][3][3][3] = {{{{0.0}}}};
  //  for(int i=0;i<3;++i)
  //    for(int j=0;j<3;++j)
  //      for(int k=0;k<3;++k)
  //        for(int l=0;l<3;++l)
  //          for(int A=0;A<3;++A)
  //            for(int B=0;B<3;++B)
  //              for(int C=0;C<3;++C)
  //                for(int D=0;D<3;++D)
  //                  tmp[i][j][k][l] +=
  //                  defgr(i,A)*defgr(j,B)*defgr(k,C)*defgr(l,D)*CMAT[A][B][C][D];
  //        }
  //  PrintFourTensor(tmp);
  //  LINALG::Matrix<6,6> CResult1(true);
  //  Setup6x6VoigtMatrix(CResult1,tmp1);

  // But we can use the fact that CResult(i,j,k,l)=CResult(k,l,i,j) iff we have a hyperelatic
  // material
  LINALG::Matrix<6, 6> CResult(true);

  CResult(0, 0) = PullBack4Tensorijkl(defgr, CMAT, 0, 0, 0, 0);
  CResult(0, 1) = PullBack4Tensorijkl(defgr, CMAT, 0, 0, 1, 1);
  CResult(0, 2) = PullBack4Tensorijkl(defgr, CMAT, 0, 0, 2, 2);
  CResult(0, 3) = PullBack4Tensorijkl(defgr, CMAT, 0, 0, 0, 1);
  CResult(0, 4) = PullBack4Tensorijkl(defgr, CMAT, 0, 0, 1, 2);
  CResult(0, 5) = PullBack4Tensorijkl(defgr, CMAT, 0, 0, 0, 2);
  CResult(1, 0) = CResult(0, 1);
  CResult(1, 1) = PullBack4Tensorijkl(defgr, CMAT, 1, 1, 1, 1);
  CResult(1, 2) = PullBack4Tensorijkl(defgr, CMAT, 1, 1, 2, 2);
  CResult(1, 3) = PullBack4Tensorijkl(defgr, CMAT, 1, 1, 0, 1);
  CResult(1, 4) = PullBack4Tensorijkl(defgr, CMAT, 1, 1, 1, 2);
  CResult(1, 5) = PullBack4Tensorijkl(defgr, CMAT, 1, 1, 0, 2);
  CResult(2, 0) = CResult(0, 2);
  CResult(2, 1) = CResult(1, 2);
  CResult(2, 2) = PullBack4Tensorijkl(defgr, CMAT, 2, 2, 2, 2);
  CResult(2, 3) = PullBack4Tensorijkl(defgr, CMAT, 2, 2, 0, 1);
  CResult(2, 4) = PullBack4Tensorijkl(defgr, CMAT, 2, 2, 1, 2);
  CResult(2, 5) = PullBack4Tensorijkl(defgr, CMAT, 2, 2, 0, 2);
  CResult(3, 0) = CResult(0, 3);
  CResult(3, 1) = CResult(1, 3);
  CResult(3, 2) = CResult(2, 3);
  CResult(3, 3) = PullBack4Tensorijkl(defgr, CMAT, 0, 1, 0, 1);
  CResult(3, 4) = PullBack4Tensorijkl(defgr, CMAT, 0, 1, 1, 2);
  CResult(3, 5) = PullBack4Tensorijkl(defgr, CMAT, 0, 1, 0, 2);
  CResult(4, 0) = CResult(0, 4);
  CResult(4, 1) = CResult(1, 4);
  CResult(4, 2) = CResult(2, 4);
  CResult(4, 3) = CResult(3, 4);
  CResult(4, 4) = PullBack4Tensorijkl(defgr, CMAT, 1, 2, 1, 2);
  CResult(4, 5) = PullBack4Tensorijkl(defgr, CMAT, 1, 2, 0, 2);
  CResult(5, 0) = CResult(0, 5);
  CResult(5, 1) = CResult(1, 5);
  CResult(5, 2) = CResult(2, 5);
  CResult(5, 3) = CResult(3, 5);
  CResult(5, 4) = CResult(4, 5);
  CResult(5, 5) = PullBack4Tensorijkl(defgr, CMAT, 0, 2, 0, 2);

  return CResult;
}

/*-------------------------------------------------------------------------------------*
 | pull back the ijkl-th entry of a symmetric elastic 4th order                        |
 | tensor (in matrix/voigt notation) /// via the 2th order deformation                 |
 | gradient (also in matrix notation)                                      Thon  01/15 |
 *-------------------------------------------------------------------------------------*/
double MAT::GrowthVolumetric::PullBack4Tensorijkl(const LINALG::Matrix<3, 3>& defgr,
    const double (&FourTensor)[3][3][3][3], const double& i, const double& j, const double& k,
    const double& l)
{
  double CResult_ijkl = 0;

  for (int A = 0; A < 3; ++A)
    for (int B = 0; B < 3; ++B)
      for (int C = 0; C < 3; ++C)
        for (int D = 0; D < 3; ++D)
          CResult_ijkl +=
              defgr(i, A) * defgr(j, B) * defgr(k, C) * defgr(l, D) * FourTensor[A][B][C][D];

  return CResult_ijkl;
}

/*-------------------------------------------------------------------------------------*
 |  Setup 4-Tensor from 6x6 Voigt notation                                 thon  01/15 |
 *-------------------------------------------------------------------------------------*/
void MAT::GrowthVolumetric::Setup4Tensor(
    double (&FourTensor)[3][3][3][3], const LINALG::Matrix<6, 6>& VoigtMatrix)
{
  // Clear4Tensor(FourTensor);
  // Setup 4-Tensor from 6x6 Voigt matrix (which has to be the representative of a 4 tensor with at
  // least minor symmetries)
  FourTensor[0][0][0][0] = VoigtMatrix(0, 0);  // C1111
  FourTensor[0][0][1][1] = VoigtMatrix(0, 1);  // C1122
  FourTensor[0][0][2][2] = VoigtMatrix(0, 2);  // C1133
  FourTensor[0][0][0][1] = VoigtMatrix(0, 3);  // C1112
  FourTensor[0][0][1][0] = VoigtMatrix(0, 3);  // C1121
  FourTensor[0][0][1][2] = VoigtMatrix(0, 4);  // C1123
  FourTensor[0][0][2][1] = VoigtMatrix(0, 4);  // C1132
  FourTensor[0][0][0][2] = VoigtMatrix(0, 5);  // C1113
  FourTensor[0][0][2][0] = VoigtMatrix(0, 5);  // C1131

  FourTensor[1][1][0][0] = VoigtMatrix(1, 0);  // C2211
  FourTensor[1][1][1][1] = VoigtMatrix(1, 1);  // C2222
  FourTensor[1][1][2][2] = VoigtMatrix(1, 2);  // C2233
  FourTensor[1][1][0][1] = VoigtMatrix(1, 3);  // C2212
  FourTensor[1][1][1][0] = VoigtMatrix(1, 3);  // C2221
  FourTensor[1][1][1][2] = VoigtMatrix(1, 4);  // C2223
  FourTensor[1][1][2][1] = VoigtMatrix(1, 4);  // C2232
  FourTensor[1][1][0][2] = VoigtMatrix(1, 5);  // C2213
  FourTensor[1][1][2][0] = VoigtMatrix(1, 5);  // C2231

  FourTensor[2][2][0][0] = VoigtMatrix(2, 0);  // C3311
  FourTensor[2][2][1][1] = VoigtMatrix(2, 1);  // C3322
  FourTensor[2][2][2][2] = VoigtMatrix(2, 2);  // C3333
  FourTensor[2][2][0][1] = VoigtMatrix(2, 3);  // C3312
  FourTensor[2][2][1][0] = VoigtMatrix(2, 3);  // C3321
  FourTensor[2][2][1][2] = VoigtMatrix(2, 4);  // C3323
  FourTensor[2][2][2][1] = VoigtMatrix(2, 4);  // C3332
  FourTensor[2][2][0][2] = VoigtMatrix(2, 5);  // C3313
  FourTensor[2][2][2][0] = VoigtMatrix(2, 5);  // C3331

  FourTensor[0][1][0][0] = VoigtMatrix(3, 0);
  FourTensor[1][0][0][0] = VoigtMatrix(3, 0);  // C1211 = C2111
  FourTensor[0][1][1][1] = VoigtMatrix(3, 1);
  FourTensor[1][0][1][1] = VoigtMatrix(3, 1);  // C1222 = C2122
  FourTensor[0][1][2][2] = VoigtMatrix(3, 2);
  FourTensor[1][0][2][2] = VoigtMatrix(3, 2);  // C1233 = C2133
  FourTensor[0][1][0][1] = VoigtMatrix(3, 3);
  FourTensor[1][0][0][1] = VoigtMatrix(3, 3);  // C1212 = C2112
  FourTensor[0][1][1][0] = VoigtMatrix(3, 3);
  FourTensor[1][0][1][0] = VoigtMatrix(3, 3);  // C1221 = C2121
  FourTensor[0][1][1][2] = VoigtMatrix(3, 4);
  FourTensor[1][0][1][2] = VoigtMatrix(3, 4);  // C1223 = C2123
  FourTensor[0][1][2][1] = VoigtMatrix(3, 4);
  FourTensor[1][0][2][1] = VoigtMatrix(3, 4);  // C1232 = C2132
  FourTensor[0][1][0][2] = VoigtMatrix(3, 5);
  FourTensor[1][0][0][2] = VoigtMatrix(3, 5);  // C1213 = C2113
  FourTensor[0][1][2][0] = VoigtMatrix(3, 5);
  FourTensor[1][0][2][0] = VoigtMatrix(3, 5);  // C1231 = C2131

  FourTensor[1][2][0][0] = VoigtMatrix(4, 0);
  FourTensor[2][1][0][0] = VoigtMatrix(4, 0);  // C2311 = C3211
  FourTensor[1][2][1][1] = VoigtMatrix(4, 1);
  FourTensor[2][1][1][1] = VoigtMatrix(4, 1);  // C2322 = C3222
  FourTensor[1][2][2][2] = VoigtMatrix(4, 2);
  FourTensor[2][1][2][2] = VoigtMatrix(4, 2);  // C2333 = C3233
  FourTensor[1][2][0][1] = VoigtMatrix(4, 3);
  FourTensor[2][1][0][1] = VoigtMatrix(4, 3);  // C2312 = C3212
  FourTensor[1][2][1][0] = VoigtMatrix(4, 3);
  FourTensor[2][1][1][0] = VoigtMatrix(4, 3);  // C2321 = C3221
  FourTensor[1][2][1][2] = VoigtMatrix(4, 4);
  FourTensor[2][1][1][2] = VoigtMatrix(4, 4);  // C2323 = C3223
  FourTensor[1][2][2][1] = VoigtMatrix(4, 4);
  FourTensor[2][1][2][1] = VoigtMatrix(4, 4);  // C2332 = C3232
  FourTensor[1][2][0][2] = VoigtMatrix(4, 5);
  FourTensor[2][1][0][2] = VoigtMatrix(4, 5);  // C2313 = C3213
  FourTensor[1][2][2][0] = VoigtMatrix(4, 5);
  FourTensor[2][1][2][0] = VoigtMatrix(4, 5);  // C2331 = C3231

  FourTensor[0][2][0][0] = VoigtMatrix(5, 0);
  FourTensor[2][0][0][0] = VoigtMatrix(5, 0);  // C1311 = C3111
  FourTensor[0][2][1][1] = VoigtMatrix(5, 1);
  FourTensor[2][0][1][1] = VoigtMatrix(5, 1);  // C1322 = C3122
  FourTensor[0][2][2][2] = VoigtMatrix(5, 2);
  FourTensor[2][0][2][2] = VoigtMatrix(5, 2);  // C1333 = C3133
  FourTensor[0][2][0][1] = VoigtMatrix(5, 3);
  FourTensor[2][0][0][1] = VoigtMatrix(5, 3);  // C1312 = C3112
  FourTensor[0][2][1][0] = VoigtMatrix(5, 3);
  FourTensor[2][0][1][0] = VoigtMatrix(5, 3);  // C1321 = C3121
  FourTensor[0][2][1][2] = VoigtMatrix(5, 4);
  FourTensor[2][0][1][2] = VoigtMatrix(5, 4);  // C1323 = C3123
  FourTensor[0][2][2][1] = VoigtMatrix(5, 4);
  FourTensor[2][0][2][1] = VoigtMatrix(5, 4);  // C1332 = C3132
  FourTensor[0][2][0][2] = VoigtMatrix(5, 5);
  FourTensor[2][0][0][2] = VoigtMatrix(5, 5);  // C1313 = C3113
  FourTensor[0][2][2][0] = VoigtMatrix(5, 5);
  FourTensor[2][0][2][0] = VoigtMatrix(5, 5);  // C1331 = C3131

}  // Setup4Tensor()


///*------------------------------------------------------------------------------------------*
// |  Setup 6x6 matrix in Voigt notation from 4-Tensor                            thon  01/15 |
// *------------------------------------------------------------------------------------------*/
// void MAT::GrowthScdACRadial::Setup6x6VoigtMatrix(
//    LINALG::Matrix<6,6>& VoigtMatrix,
//    const double (&FourTensor)[3][3][3][3]
//)
//{
/////*  [      C1111                 C1122                C1133                0.5*(C1112+C1121)
/// 0.5*(C1123+C1132)                0.5*(C1113+C1131)      ] /    [      C2211 C2222 C2233
/// 0.5*(C2212+C2221)               0.5*(C2223+C2232) 0.5*(C2213+C2231)      ] /    [      C3311
/// C3322                C3333                0.5*(C3312+C3321) 0.5*(C3323+C3332) 0.5*(C3313+C3331)
/// ] /    [0.5*(C1211+C2111) 0.5*(C1222+C2122)    0.5*(C1233+C2133) 0.5*(C1212+C2112+C1221+C2121)
/// 0.5*(C1223+C2123+C1232+C2132)    0.5*(C1213+C2113+C1231+C2131)] /    [0.5*(C2311+C3211)
/// 0.5*(C2322+C3222)    0.5*(C2333+C3233)    0.5*(C2312+C3212+C2321+C3221)
/// 0.5*(C2323+C3223+C2332+C3232)    0.5*(C2313+C3213+C2331+C3231)] /    [0.5*(C1322+C3122)
/// 0.5*(C1322+C3122)    0.5*(C1333+C3133)    0.5*(C1312+C3112+C1321+C3121)
/// 0.5*(C1323+C3123+C1332+C3132)    0.5*(C1313+C3113+C1331+C3131)] */
//
//  // Setup 4-Tensor from 6x6 Voigt matrix
//  VoigtMatrix(0,0) = FourTensor[0][0][0][0]; //C1111
//  VoigtMatrix(0,1) = FourTensor[0][0][1][1]; //C1122
//  VoigtMatrix(0,2) = FourTensor[0][0][2][2]; //C1133
//  VoigtMatrix(0,3) = 0.5 * (FourTensor[0][0][0][1] + FourTensor[0][0][1][0]); //0.5*(C1112+C1121)
//  VoigtMatrix(0,4) = 0.5 * (FourTensor[0][0][1][2] + FourTensor[0][0][2][1]); //0.5*(C1123+C1132)
//  VoigtMatrix(0,5) = 0.5 * (FourTensor[0][0][0][2] + FourTensor[0][0][2][0]); //0.5*(C1113+C1131)
//
//  VoigtMatrix(1,0) = FourTensor[1][1][0][0]; //C2211
//  VoigtMatrix(1,1) = FourTensor[1][1][1][1]; //C2222
//  VoigtMatrix(1,2) = FourTensor[1][1][2][2]; //C2233
//  VoigtMatrix(1,3) = 0.5 * (FourTensor[1][1][0][1] + FourTensor[1][1][1][0]); //0.5*(C2212+C2221)
//  VoigtMatrix(1,4) = 0.5 * (FourTensor[1][1][1][2] + FourTensor[1][1][2][1]); //0.5*(C2223+C2232)
//  VoigtMatrix(1,5) = 0.5 * (FourTensor[1][1][0][2] + FourTensor[1][1][2][0]); //0.5*(C2213+C2231)
//
//  VoigtMatrix(2,0) = FourTensor[2][2][0][0]; //C3311
//  VoigtMatrix(2,1) = FourTensor[2][2][1][1]; //C3322
//  VoigtMatrix(2,2) = FourTensor[2][2][2][2]; //C3333
//  VoigtMatrix(2,3) = 0.5 * (FourTensor[2][2][0][1] + FourTensor[2][2][1][0]); //0.5*(C3312+C3321)
//  VoigtMatrix(2,4) = 0.5 * (FourTensor[2][2][1][2] + FourTensor[2][2][2][1]); //0.5*(C3323+C3332)
//  VoigtMatrix(2,5) = 0.5 * (FourTensor[2][2][0][2] + FourTensor[2][2][2][0]); //0.5*(C3313+C3331)
//
//  VoigtMatrix(3,0) = 0.5 * (FourTensor[0][1][0][0] + FourTensor[1][0][0][0]); //0.5*(C1211+C2111)
//  VoigtMatrix(3,1) = 0.5 * (FourTensor[0][1][1][1] + FourTensor[1][0][1][1]); //0.5*(C1222+C2122)
//  VoigtMatrix(3,2) = 0.5 * (FourTensor[0][1][2][2] + FourTensor[1][0][2][2]); //0.5*(C1233+C2133)
//  VoigtMatrix(3,3) = 0.25 * (FourTensor[0][1][0][1] + FourTensor[1][0][0][1] +
//  FourTensor[0][1][1][0] + FourTensor[1][0][1][0]); //0.5*(C1212+C2112+C1221+C2121)
//  VoigtMatrix(3,4) = 0.25 * (FourTensor[0][1][1][2] + FourTensor[1][0][1][2] +
//  FourTensor[0][1][2][1] + FourTensor[1][0][2][1]); //0.5*(C1223+C2123+C1232+C2132)
//  VoigtMatrix(3,5) = 0.25 * (FourTensor[0][1][0][2] + FourTensor[1][0][0][2] +
//  FourTensor[0][1][2][0] + FourTensor[1][0][2][0]); //0.5*(C1213+C2113+C1231+C2131)
//
//  VoigtMatrix(4,0) = 0.5 * (FourTensor[1][2][0][0] + FourTensor[2][1][0][0]); //0.5*(C2311+C3211)
//  VoigtMatrix(4,1) = 0.5 * (FourTensor[1][2][1][1] + FourTensor[2][1][1][1]); //0.5*(C2322+C3222)
//  VoigtMatrix(4,2) = 0.5 * (FourTensor[1][2][2][2] + FourTensor[2][1][2][2]); //0.5*(C2333+C3233)
//  VoigtMatrix(4,3) = 0.25 * (FourTensor[1][2][0][1] + FourTensor[2][1][0][1] +
//  FourTensor[1][2][1][0] + FourTensor[2][1][1][0]); //0.5*(C2312+C3212+C2321+C3221)
//  VoigtMatrix(4,4) = 0.25 * (FourTensor[1][2][1][2] + FourTensor[2][1][1][2] +
//  FourTensor[1][2][2][1] + FourTensor[2][1][2][1]); //0.5*(C2323+C3223+C2332+C3232)
//  VoigtMatrix(4,5) = 0.25 * (FourTensor[1][2][0][2] + FourTensor[2][1][0][2] +
//  FourTensor[1][2][2][0] + FourTensor[2][1][2][0]); //0.5*(C2313+C3213+C2331+C3231)
//
//  VoigtMatrix(5,0) = 0.5 * (FourTensor[0][2][0][0] + FourTensor[2][0][0][0]); //0.5*(C1311+C3111)
//  VoigtMatrix(5,1) = 0.5 * (FourTensor[0][2][1][1] + FourTensor[2][0][1][1]); //0.5*(C1322+C3122)
//  VoigtMatrix(5,2) = 0.5 * (FourTensor[0][2][2][2] + FourTensor[2][0][2][2]); //0.5*(C1333+C3133)
//  VoigtMatrix(5,3) = 0.25 * (FourTensor[0][2][0][1] + FourTensor[2][0][0][1] +
//  FourTensor[0][2][1][0] + FourTensor[2][0][1][0]); //0.5*(C1312+C3112+C1321+C3121)
//  VoigtMatrix(5,4) = 0.25 * (FourTensor[0][2][1][2] + FourTensor[2][0][1][2] +
//  FourTensor[0][2][2][1] + FourTensor[2][0][2][1]); //0.5*(C1323+C3123+C1332+C3132)
//  VoigtMatrix(5,5) = 0.25 * (FourTensor[0][2][0][2] + FourTensor[2][0][0][2] +
//  FourTensor[0][2][2][0] + FourTensor[2][0][2][0]); //0.5*(C1313+C3113+C1331+C3131)
//
//}  // Setup6x6VoigtMatrix()

/*------------------------------------------------------------------------------------------*
 |  Print Four Tensor                                                           thon  01/15 |
 *------------------------------------------------------------------------------------------*/
void MAT::GrowthVolumetric::PrintFourTensor(double (&FourTensor)[3][3][3][3])
{
  std::cout << "-----------------Print Four Tensor--------------" << std::endl;

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l)
          std::cout << "ELEMENT " << i << j << k << l << " : " << FourTensor[i][j][k][l]
                    << std::endl;

  std::cout << "------------------------------------------------" << std::endl;
  return;
}


/*----------------------------------------------------------------------------*/
void MAT::GrowthVolumetric::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  MAT::PAR::Growth* params = Parameter();
  if (params != NULL) matid = params->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  int numgp = 0;
  if (isinit_)
  {
    numgp = theta_->size();  // size is number of Gauss points
  }
  AddtoPack(data, numgp);

  // Pack internal variables
  for (int gp = 0; gp < numgp; ++gp)
  {
    AddtoPack(data, tr_mandel_e_->at(gp));
  }
  for (int gp = 0; gp < numgp; ++gp)
  {
    AddtoPack(data, lambda_fib_e_->at(gp));
  }

  // Pack prescribed const growth trigger
  AddtoPack(data, growthtrig_const_);

  // Pack internal variables
  for (int i = 0; i < 3; ++i)
  {
    AddtoPack(data, refdir_(i, 0));
  }

  // int numgp=0;
  if (isinit_)
  {
    numgp = theta_->size();
    ;  // size is number of Fauss points
  }
  AddtoPack(data, numgp);

  for (int gp = 0; gp < numgp; gp++)
  {
    LINALG::Matrix<3, 3> F_g_hist = F_g_hist_.at(gp);
    LINALG::Matrix<3, 1> curdir = curdir_.at(gp);
    LINALG::Matrix<3, 1> curdir_for_update = curdir_for_update_.at(gp);

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        AddtoPack(data, F_g_hist(i, j));
      }
      AddtoPack(data, curdir(i, 0));
      AddtoPack(data, curdir_for_update(i, 0));
    }
  }

  // Pack base class material
  Growth::Pack(data);

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthVolumetric::Unpack(const std::vector<char>& data)
{
  isinit_ = true;
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);

  paramsVolumetric_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        paramsVolumetric_ = dynamic_cast<MAT::PAR::Growth*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  int numgp;
  ExtractfromPack(position, data, numgp);
  if (numgp == 0)
  {  // no history data to unpack
    isinit_ = false;
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d", data.size(), position);
    return;
  }

  // unpack growth internal variables
  tr_mandel_e_ = Teuchos::rcp(new std::vector<double>(numgp));
  for (int gp = 0; gp < numgp; ++gp)
  {
    double a;
    ExtractfromPack(position, data, a);
    tr_mandel_e_->at(gp) = a;
  }

  lambda_fib_e_ = Teuchos::rcp(new std::vector<double>(numgp));
  for (int gp = 0; gp < numgp; ++gp)
  {
    double a;
    ExtractfromPack(position, data, a);
    lambda_fib_e_->at(gp) = a;
  }

  double growthtrigconst;
  ExtractfromPack(position, data, growthtrigconst);
  growthtrig_const_ = growthtrigconst;

  // unpack growth internal variables
  // Pack internal variables
  for (int i = 0; i < 3; ++i)
  {
    double refdir_i;
    ExtractfromPack(position, data, refdir_i);
    refdir_(i, 0) = refdir_i;
  }

  // int numgp;
  ExtractfromPack(position, data, numgp);
  if (not(numgp == 0))
  {
    F_g_hist_ = std::vector<LINALG::Matrix<3, 3>>(numgp, LINALG::Matrix<3, 3>(true));
    curdir_ = std::vector<LINALG::Matrix<3, 1>>(numgp, LINALG::Matrix<3, 1>(true));
    curdir_for_update_ = std::vector<LINALG::Matrix<3, 1>>(numgp, LINALG::Matrix<3, 1>(true));

    for (int gp = 0; gp < numgp; gp++)
    {
      LINALG::Matrix<3, 3> F_g_hist(true);
      LINALG::Matrix<3, 1> curdir(true);
      LINALG::Matrix<3, 1> curdir_for_update(true);

      for (int i = 0; i < 3; ++i)
      {
        for (int j = 0; j < 3; ++j)
        {
          double F_g_hist_ij;
          ExtractfromPack(position, data, F_g_hist_ij);
          F_g_hist(i, j) = F_g_hist_ij;
        }
        double curdir_i;
        ExtractfromPack(position, data, curdir_i);
        curdir(i, 0) = curdir_i;

        double curdir_for_update_i;
        ExtractfromPack(position, data, curdir_for_update_i);
        curdir_for_update(i, 0) = curdir_for_update_i;
      }
      F_g_hist_.at(gp) = F_g_hist;
      curdir_.at(gp) = curdir;
      curdir_for_update_.at(gp) = curdir_for_update;
    }
  }

  // extract base class material
  std::vector<char> basedata(0);
  Growth::ExtractfromPack(position, data, basedata);
  Growth::Unpack(basedata);

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);

  return;
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthVolumetric::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  tr_mandel_e_ = Teuchos::rcp(new std::vector<double>(numgp));
  lambda_fib_e_ = Teuchos::rcp(new std::vector<double>(numgp));
  for (int j = 0; j < numgp; ++j)
  {
    tr_mandel_e_->at(j) = 0.0;
    lambda_fib_e_->at(j) = 0.0;
  }

  // setup specific anisotropic growth laws
  switch (Parameter()->growthlaw_->MaterialType())
  {
    case INPAR::MAT::m_growth_ac_radial:
    case INPAR::MAT::m_growth_ac_radial_refconc:
    {
      // CIR-AXI-RAD nomenclature
      if (not(linedef->HaveNamed("RAD")))
        dserror(
            "If you want growth into the radial direction you need to specify RAD in your input "
            "file!");

      ReadFiber(linedef, "RAD", refdir_);
      curdir_ = std::vector<LINALG::Matrix<3, 1>>(numgp, refdir_);
      curdir_for_update_ = std::vector<LINALG::Matrix<3, 1>>(numgp, refdir_);

      LINALG::Matrix<3, 3> Id(true);
      for (int i = 0; i < 3; i++) Id(i, i) = 1.0;
      F_g_hist_ = std::vector<LINALG::Matrix<3, 3>>(numgp, Id);
    }
    break;
    case INPAR::MAT::m_growth_aniso_strain:
    case INPAR::MAT::m_growth_aniso_stress:
    {
      // FIBER1 nomenclature
      if (not(linedef->HaveNamed("FIBER1")))
        dserror(
            "If you want growth in fiber direction you need to specify FIBER1 in your input file!");

      ReadFiber(linedef, "FIBER1", refdir_);

      // only refdir is used - rest remains unused...
      curdir_ = std::vector<LINALG::Matrix<3, 1>>(numgp, refdir_);
      curdir_for_update_ = std::vector<LINALG::Matrix<3, 1>>(numgp, refdir_);
      LINALG::Matrix<3, 3> Id(true);
      for (int i = 0; i < 3; i++) Id(i, i) = 1.0;
      F_g_hist_ = std::vector<LINALG::Matrix<3, 3>>(numgp, Id);
    }
    break;
    case INPAR::MAT::m_growth_aniso_strain_const_trig:
    case INPAR::MAT::m_growth_aniso_stress_const_trig:
    {
      // FIBER1 nomenclature
      if (not(linedef->HaveNamed("FIBER1")))
        dserror(
            "If you want growth in fiber direction you need to specify FIBER1 in your input file!");

      ReadFiber(linedef, "FIBER1", refdir_);

      linedef->ExtractDouble("GROWTHTRIG", growthtrig_const_);
      if (not(linedef->HaveNamed("GROWTHTRIG")))
        dserror("You need to specify GROWTHTRIG in your input file!");

      // only refdir is used - rest remains unused...
      curdir_ = std::vector<LINALG::Matrix<3, 1>>(numgp, refdir_);
      curdir_for_update_ = std::vector<LINALG::Matrix<3, 1>>(numgp, refdir_);
      LINALG::Matrix<3, 3> Id(true);
      for (int i = 0; i < 3; i++) Id(i, i) = 1.0;
      F_g_hist_ = std::vector<LINALG::Matrix<3, 3>>(numgp, Id);
    }
    break;
    default:
    {
      // directions are unused
      refdir_(true);
      curdir_ = std::vector<LINALG::Matrix<3, 1>>(numgp, refdir_);
      curdir_for_update_ = std::vector<LINALG::Matrix<3, 1>>(numgp, refdir_);
      LINALG::Matrix<3, 3> Id(true);
      for (int i = 0; i < 3; i++) Id(i, i) = 1.0;
      F_g_hist_ = std::vector<LINALG::Matrix<3, 3>>(numgp, Id);
    }
    break;
  }

  // setup base class
  Growth::Setup(numgp, linedef);

  return;
}



/*----------------------------------------------------------------------------*/
void MAT::GrowthVolumetric::Update()
{
  const int numgp = theta_->size();

  // setup anisotropic growth laws
  switch (Parameter()->growthlaw_->MaterialType())
  {
    case INPAR::MAT::m_growth_ac_radial:
    case INPAR::MAT::m_growth_ac_radial_refconc:
    {
      const LINALG::Matrix<3, 3> dummydefgrad(true);

      for (int gp = 0; gp < numgp; gp++)
      {
        LINALG::Matrix<3, 3> F_g_hist_new(true);

        Parameter()->growthlaw_->CalcFg(theta_->at(gp), ThetaOld()->at(gp), gp, &dummydefgrad,
            refdir_, curdir_, F_g_hist_, F_g_hist_new);

        F_g_hist_.at(gp) = F_g_hist_new;

        curdir_.at(gp) = curdir_for_update_.at(gp);
      }
    }
    break;
    default:
      break;
  }

  // update base class
  Growth::Update();

  return;
}



/*----------------------------------------------------------------------*
 | Function which reads in the given fiber value             Thon 01/15 |
 *----------------------------------------------------------------------*/
void MAT::GrowthVolumetric::ReadFiber(
    DRT::INPUT::LineDefinition* linedef, std::string specifier, LINALG::Matrix<3, 1>& fiber_vector)
{
  std::vector<double> fiber1;
  linedef->ExtractDoubleVector(specifier, fiber1);
  double f1norm = 0.;
  // normalization
  for (int i = 0; i < 3; ++i)
  {
    f1norm += fiber1[i] * fiber1[i];
  }
  f1norm = sqrt(f1norm);

  // fill final fiber vector
  for (int i = 0; i < 3; ++i) fiber_vector(i) = fiber1[i] / f1norm;
}
