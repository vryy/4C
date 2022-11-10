/*----------------------------------------------------------------------*/
/*! \file
\brief This file contains routines for integration point based isotropic and anisotropic volumetric
growth laws.

\level 2

 *----------------------------------------------------------------------*/


#include <utility>

#include "growth.H"

#include "growth_law.H"
#include "material_service.H"

#include "drt_globalproblem.H"
#include "matpar_bundle.H"
#include "drt_utils_factory.H"
#include "voigt_notation.H"

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
    dserror("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (DRT::Problem::Instance(probinst)->Materials()->Num() == 0)
    dserror("List of materials in the global problem instance is empty.");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(idgrowthlaw_);

  switch (curmat->Type())
  {
    case INPAR::MAT::m_growth_aniso_strain:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::GrowthLawAnisoStrain(curmat));
      auto* params = static_cast<MAT::PAR::GrowthLawAnisoStrain*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case INPAR::MAT::m_growth_aniso_stress:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::GrowthLawAnisoStress(curmat));
      auto* params = static_cast<MAT::PAR::GrowthLawAnisoStress*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case INPAR::MAT::m_growth_aniso_strain_const_trig:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::GrowthLawAnisoStrainConstTrig(curmat));
      auto* params = static_cast<MAT::PAR::GrowthLawAnisoStrainConstTrig*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case INPAR::MAT::m_growth_aniso_stress_const_trig:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::GrowthLawAnisoStressConstTrig(curmat));
      auto* params = static_cast<MAT::PAR::GrowthLawAnisoStressConstTrig*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case INPAR::MAT::m_growth_iso_stress:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::GrowthLawIsoStress(curmat));
      auto* params = static_cast<MAT::PAR::GrowthLawIsoStress*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case INPAR::MAT::m_growth_ac:
    {
      if (curmat->Parameter() == nullptr) curmat->SetParameter(new MAT::PAR::GrowthLawAC(curmat));
      auto* params = static_cast<MAT::PAR::GrowthLawAC*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case INPAR::MAT::m_growth_ac_radial:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::GrowthLawACRadial(curmat));
      auto* params = static_cast<MAT::PAR::GrowthLawACRadial*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case INPAR::MAT::m_growth_ac_radial_refconc:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::GrowthLawACRadialRefConc(curmat));
      auto* params = static_cast<MAT::PAR::GrowthLawACRadialRefConc*>(curmat->Parameter());
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case INPAR::MAT::m_growth_const:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::GrowthLawConst(curmat));
      auto* params = static_cast<MAT::PAR::GrowthLawConst*>(curmat->Parameter());
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
      params_(nullptr),
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
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
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

  params_ = nullptr;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
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
    auto* matel = dynamic_cast<MAT::So3Material*>(o);
    if (matel == nullptr) dserror("failed to unpack elastic material");
    matelastic_ = Teuchos::rcp(matel);
  }
  else
    matelastic_ = Teuchos::null;

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
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
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  Matelastic()->Evaluate(defgrd, glstrain, params, stress, cmat, gp, eleGID);
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
  auto* grow = new MAT::GrowthVolumetric();
  grow->Unpack(data);
  return grow;
}

/*----------------------------------------------------------------------------*/
MAT::GrowthVolumetric::GrowthVolumetric()
    : Growth(),
      tr_mandel_e_(Teuchos::null),
      lambda_fib_e_(Teuchos::null),
      growthtrig_const_(0.0),
      paramsVolumetric_(nullptr),
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
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, const int gp, const int eleGID)
{
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
    EvaluateGrowth(&theta, &dthetadCvec, defgrd, glstrain, params, gp, eleGID);

    // modify the parameter list to be passed to the elastic material
    Teuchos::ParameterList paramselast(params);
    paramselast.remove("matparderiv", false);

    LINALG::Matrix<6, 1> S(true);
    LINALG::Matrix<6, 6> cmatdach(true);

    GetSAndCmatdach(theta, defgrd, &S, &cmatdach, paramselast, gp, eleGID);

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
    cmatelast = MAT::PullBackFourTensor<3>(F_ginv, cmatdach);

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

    GetSAndCmatdach(theta + espilon, defgrd, &SEps, &cmatdachEps, params, gp, eleGID);

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
          &dthetadp, ThetaOldAtGp(gp), Matelastic(), defgrd, glstrain, params, eleGID);
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

    GetSAndCmatdach(theta, defgrd, &Svec, &cmatdach, params, gp, eleGID);

    *stress = Svec;

    // calculate growth part F_g of the deformation gradient F
    LINALG::Matrix<3, 3> F_g(true);

    Parameter()->growthlaw_->CalcFg(
        theta, ThetaOldAtGp(gp), gp, defgrd, refdir_, curdir_, F_g_hist_, F_g);

    // calculate F_g^(-1)
    LINALG::Matrix<3, 3> F_ginv(true);
    F_ginv.Invert(F_g);

    // constitutive matrix including growth cmat = F_g^-1 F_g^-1 cmatdach F_g^-T F_g^-T
    LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatelast(true);
    cmatelast = MAT::PullBackFourTensor<3>(F_ginv, cmatdach);

    *cmat = cmatelast;
  }
  else
  {
    EvaluateElastic(defgrd, glstrain, stress, cmat, params, gp, eleGID);
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
    UTILS::VOIGT::Strains::VectorToMatrix(Cvec, C);

    LINALG::Matrix<3, 1> CDir(true);
    CDir.MultiplyNN(1.0, C, refdir_);
    lambda_fib_e_->at(gp) =
        sqrt(CDir(0) * refdir_(0) + CDir(1) * refdir_(1) + CDir(2) * refdir_(2));
  }
}

/*----------------------------------------------------------------------------*/
void MAT::GrowthVolumetric::EvaluateGrowth(double* theta, LINALG::Matrix<6, 1>* dthetadC,
    const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  // get gauss point number
  if (gp == -1) dserror("No Gauss point number provided in material.");

  double thetaold = ThetaOldAtGp(gp);

  MAT::Growth* matgrowth = this;
  Parameter()->growthlaw_->Evaluate(theta, thetaold, dthetadC, *matgrowth, defgrd, glstrain,
      refdir_, curdir_, F_g_hist_, growthtrig_const_, params, gp, eleGID);
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
    const int gp, const int eleGID)
{
  const double eps = 1.0e-14;
  const double starttime = Parameter()->starttime_;
  const double endtime = Parameter()->endtime_;
  const double time = params.get<double>("total time", -1.0);

  if (((time > starttime + eps) and (time <= endtime + eps)) or
      ((starttime < 0.0) and (endtime < 0.0)))
  {
    double theta = theta_->at(gp);
    double thetaold = ThetaOld()->at(gp);

    MAT::Growth* matgrowth = this;
    Parameter()->growthlaw_->Evaluate(&theta, thetaold, linmass_disp, *matgrowth, defgrd, glstrain,
        refdir_, curdir_, F_g_hist_, growthtrig_const_, params, gp, eleGID);

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
    const int gp, const int eleGID)
{
  // calculate growth part F_g of the deformation gradient F
  LINALG::Matrix<3, 3> F_g(true);
  Parameter()->growthlaw_->CalcFg(
      theta, ThetaOldAtGp(gp), gp, defgrd, refdir_, curdir_, F_g_hist_, F_g);

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
  UTILS::VOIGT::Strains::MatrixToVector(Cdach, Cdachvec);

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
  Matelastic()->Evaluate(&defgrddach, &glstraindachvec, params, &Sdachvec, cmatdach, gp, eleGID);

  // calculate stress
  // 2PK stress S = F_g^-1 Sdach F_g^-T
  LINALG::Matrix<3, 3> Sdach(true);
  UTILS::VOIGT::Stresses::VectorToMatrix(Sdachvec, Sdach);

  LINALG::Matrix<3, 3> tmp(true);
  tmp.MultiplyNT(Sdach, F_ginv);
  LINALG::Matrix<3, 3> S(true);
  S.MultiplyNN(F_ginv, tmp);

  UTILS::VOIGT::Stresses::MatrixToVector(S, *stress);

  // trace of elastic Mandel stress Mdach = Cdach Sdach
  tr_mandel_e_->at(gp) = Cdachvec(0) * Sdachvec(0) + Cdachvec(1) * Sdachvec(1) +
                         Cdachvec(2) * Sdachvec(2) + Cdachvec(3) * Sdachvec(3) +
                         Cdachvec(4) * Sdachvec(4) + Cdachvec(5) * Sdachvec(5);

  // elastic fiber stretch lambda = \sqrt(f_0 \cdot Cdach f_0)
  LINALG::Matrix<3, 1> CdachDir(true);
  CdachDir.MultiplyNN(1.0, Cdach, refdir_);
  lambda_fib_e_->at(gp) =
      sqrt(CdachDir(0) * refdir_(0) + CdachDir(1) * refdir_(1) + CdachDir(2) * refdir_(2));
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
  if (params != nullptr) matid = params->Id();  // in case we are in post-process mode
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
    numgp = theta_->size();  // size is number of Fauss points
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

  paramsVolumetric_ = nullptr;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
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
  if (numgp != 0)
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
      {
        dserror(
            "If you want growth into the radial direction you need to specify RAD in your input "
            "file!");
      }

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
}



/*----------------------------------------------------------------------*
 | Function which reads in the given fiber value             Thon 01/15 |
 *----------------------------------------------------------------------*/
void MAT::GrowthVolumetric::ReadFiber(
    DRT::INPUT::LineDefinition* linedef, std::string specifier, LINALG::Matrix<3, 1>& fiber_vector)
{
  std::vector<double> fiber1;
  linedef->ExtractDoubleVector(std::move(specifier), fiber1);
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
