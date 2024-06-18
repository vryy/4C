/*----------------------------------------------------------------------*/
/*! \file
\brief This file contains routines for integration point based isotropic and anisotropic volumetric
growth laws.

\level 2

 *----------------------------------------------------------------------*/


#include "4C_mat_growth.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_mat_growth_law.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
Mat::PAR::Growth::Growth(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      idmatelastic_(matdata.parameters.get<int>("IDMATELASTIC")),
      idgrowthlaw_(matdata.parameters.get<int>("GROWTHLAW")),
      starttime_(matdata.parameters.get<double>("STARTTIME")),
      endtime_(matdata.parameters.get<double>("ENDTIME"))
{
  // retrieve problem instance to read from
  const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (Global::Problem::Instance(probinst)->Materials() == Teuchos::null)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (Global::Problem::Instance(probinst)->Materials()->Num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  auto* curmat = Global::Problem::Instance(probinst)->Materials()->ParameterById(idgrowthlaw_);

  switch (curmat->Type())
  {
    case Core::Materials::m_growth_aniso_strain:
    {
      auto* params = static_cast<Mat::PAR::GrowthLawAnisoStrain*>(curmat);
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case Core::Materials::m_growth_aniso_stress:
    {
      auto* params = static_cast<Mat::PAR::GrowthLawAnisoStress*>(curmat);
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case Core::Materials::m_growth_aniso_strain_const_trig:
    {
      auto* params = static_cast<Mat::PAR::GrowthLawAnisoStrainConstTrig*>(curmat);
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case Core::Materials::m_growth_aniso_stress_const_trig:
    {
      auto* params = static_cast<Mat::PAR::GrowthLawAnisoStressConstTrig*>(curmat);
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case Core::Materials::m_growth_iso_stress:
    {
      auto* params = static_cast<Mat::PAR::GrowthLawIsoStress*>(curmat);
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case Core::Materials::m_growth_ac:
    {
      auto* params = static_cast<Mat::PAR::GrowthLawAC*>(curmat);
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case Core::Materials::m_growth_ac_radial:
    {
      auto* params = static_cast<Mat::PAR::GrowthLawACRadial*>(curmat);
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case Core::Materials::m_growth_ac_radial_refconc:
    {
      auto* params = static_cast<Mat::PAR::GrowthLawACRadialRefConc*>(curmat);
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    case Core::Materials::m_growth_const:
    {
      auto* params = static_cast<Mat::PAR::GrowthLawConst*>(curmat);
      growthlaw_ = params->CreateGrowthLaw();
      break;
    }
    default:
      FOUR_C_THROW("unknown material type %d", curmat->Type());
      break;
  }

  if (starttime_ > endtime_)
    FOUR_C_THROW("WTF! It is not reasonable to have a starttime that is larger than the endtime!");
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::Growth::create_material()
{
  Teuchos::RCP<Core::Mat::Material> mat;

  switch (growthlaw_->MaterialType())
  {
    case Core::Materials::m_growth_aniso_strain:
    case Core::Materials::m_growth_aniso_stress:
    case Core::Materials::m_growth_aniso_strain_const_trig:
    case Core::Materials::m_growth_aniso_stress_const_trig:
    case Core::Materials::m_growth_iso_stress:
    case Core::Materials::m_growth_ac:
    case Core::Materials::m_growth_ac_radial:
    case Core::Materials::m_growth_ac_radial_refconc:
    case Core::Materials::m_growth_const:
      mat = Teuchos::rcp(new Mat::GrowthVolumetric(this));
      break;
    default:
      FOUR_C_THROW(
          "The growth law you have chosen is not valid for the standard volumetric growth "
          "material");
      mat = Teuchos::null;
      break;
  }

  return mat;
}

/*----------------------------------------------------------------------------*/
Mat::Growth::Growth()
    : theta_(Teuchos::null),
      isinit_(false),
      params_(nullptr),
      matelastic_(Teuchos::null),
      thetaold_(Teuchos::null),
      histdata_()
{
}

/*----------------------------------------------------------------------------*/
Mat::Growth::Growth(Mat::PAR::Growth* params)
    : theta_(Teuchos::null),
      isinit_(false),
      params_(params),
      matelastic_(Teuchos::null),
      thetaold_(Teuchos::null),
      histdata_()
{
}

/*----------------------------------------------------------------------------*/
void Mat::Growth::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  int numgp;
  if (!isinit_)
  {
    numgp = 0;  // not initialized -> nothing to pack
  }
  else
  {
    numgp = theta_->size();  // size is number of Gauss points
  }
  add_to_pack(data, numgp);
  // Pack internal variables
  for (int gp = 0; gp < numgp; ++gp)
  {
    add_to_pack(data, thetaold_->at(gp));
    add_to_pack(data, theta_->at(gp));
  }

  add_to_pack(data, histdata_);

  // Pack data of elastic material
  if (matelastic_ != Teuchos::null)
  {
    matelastic_->pack(data);
  }
}

/*----------------------------------------------------------------------------*/
void Mat::Growth::unpack(const std::vector<char>& data)
{
  isinit_ = true;
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  extract_from_pack(position, data, matid);

  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = dynamic_cast<Mat::PAR::Growth*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  int numgp;
  extract_from_pack(position, data, numgp);
  if (numgp == 0)
  {  // no history data to unpack
    isinit_ = false;
    if (position != data.size())
      FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
    return;
  }

  // unpack growth internal variables
  theta_ = Teuchos::rcp(new std::vector<double>(numgp));
  thetaold_ = Teuchos::rcp(new std::vector<double>(numgp));
  for (int gp = 0; gp < numgp; ++gp)
  {
    double a;
    extract_from_pack(position, data, a);
    thetaold_->at(gp) = a;
    extract_from_pack(position, data, a);
    theta_->at(gp) = a;
  }

  extract_from_pack(position, data, histdata_);

  // Unpack data of elastic material (these lines are copied from element.cpp)
  std::vector<char> dataelastic;
  extract_from_pack(position, data, dataelastic);
  if (dataelastic.size() > 0)
  {
    Core::Communication::ParObject* o =
        Core::Communication::Factory(dataelastic);  // Unpack is done here
    auto* matel = dynamic_cast<Mat::So3Material*>(o);
    if (matel == nullptr) FOUR_C_THROW("failed to unpack elastic material");
    matelastic_ = Teuchos::rcp(matel);
  }
  else
    matelastic_ = Teuchos::null;

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------------*/
void Mat::Growth::setup(int numgp, Input::LineDefinition* linedef)
{
  if (isinit_)
    FOUR_C_THROW("This function should just be called if the material is not yet initialized.");

  theta_ = Teuchos::rcp(new std::vector<double>(numgp));
  thetaold_ = Teuchos::rcp(new std::vector<double>(numgp));
  for (int j = 0; j < numgp; ++j)
  {
    theta_->at(j) = 1.0;
    thetaold_->at(j) = 1.0;
  }

  // Setup of elastic material
  matelastic_ = Teuchos::rcp_dynamic_cast<Mat::So3Material>(Mat::Factory(params_->idmatelastic_));
  matelastic_->setup(numgp, linedef);

  isinit_ = true;
}

/*----------------------------------------------------------------------------*/
void Mat::Growth::update()
{
  const int numgp = theta_->size();

  for (int i = 0; i < numgp; i++)
  {
    thetaold_->at(i) = theta_->at(i);
  }

  matelastic_->update();
}

/*----------------------------------------------------------------------------*/
void Mat::Growth::reset_step()
{
  const int numgp = theta_->size();

  for (int i = 0; i < numgp; i++)
  {
    theta_->at(i) = thetaold_->at(i);
  }

  matelastic_->reset_step();
}

void Mat::Growth::StoreHistory(int timestep) { histdata_[timestep] = *thetaold_; }

void Mat::Growth::SetHistory(int timestep) { *thetaold_ = histdata_.at(timestep); }

/*----------------------------------------------------------------------------*/
void Mat::Growth::EvaluateElastic(const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<6, 1>* glstrain, Core::LinAlg::Matrix<6, 1>* stress,
    Core::LinAlg::Matrix<6, 6>* cmat, Teuchos::ParameterList& params, const int gp,
    const int eleGID)
{
  Matelastic()->evaluate(defgrd, glstrain, params, stress, cmat, gp, eleGID);
}


double Mat::Growth::Density(int gp) const
{
  const double density_elast = matelastic_->Density();
  const double theta_gp = (*theta_)[gp];

  const double density_scale = Parameter()->growthlaw_->DensityScale(theta_gp);

  return density_elast * density_scale;
}

/*----------------------------------------------------------------------*
 | returns whether material density is constant (public)  schmidt 11/17 |
 *----------------------------------------------------------------------*/
bool Mat::Growth::VaryingDensity() const { return Parameter()->growthlaw_->VaryingDensity(); }


/*----------------------------------------------------------------------------*/
void Mat::GrowthVolumetric::VisNames(std::map<std::string, int>& names)
{
  std::string name = "theta";
  names[name] = 1;

  switch (Parameter()->growthlaw_->MaterialType())
  {
    case Core::Materials::m_growth_aniso_stress:
    case Core::Materials::m_growth_aniso_stress_const_trig:
    case Core::Materials::m_growth_iso_stress:
    {
      name = "tr_mandel_e";
      names[name] = 1;
    }
    break;
    case Core::Materials::m_growth_aniso_strain:
    case Core::Materials::m_growth_aniso_strain_const_trig:
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
bool Mat::GrowthVolumetric::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  if (name == "theta")
  {
    if ((int)data.size() != 1) FOUR_C_THROW("size mismatch");
    double temp = 0.0;
    for (int gp = 0; gp < numgp; gp++) temp += theta_->at(gp);
    data[0] = temp / numgp;
  }
  else if (name == "tr_mandel_e")
  {
    if ((int)data.size() != 1) FOUR_C_THROW("size mismatch");
    double temp = 0.0;
    for (int gp = 0; gp < numgp; gp++) temp += tr_mandel_e_->at(gp);
    data[0] = temp / numgp;
  }
  else if (name == "lambda_fib_e")
  {
    if ((int)data.size() != 1) FOUR_C_THROW("size mismatch");
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


Mat::GrowthVolumetricType Mat::GrowthVolumetricType::instance_;

/*----------------------------------------------------------------------------*/
Core::Communication::ParObject* Mat::GrowthVolumetricType::Create(const std::vector<char>& data)
{
  auto* grow = new Mat::GrowthVolumetric();
  grow->unpack(data);
  return grow;
}

/*----------------------------------------------------------------------------*/
Mat::GrowthVolumetric::GrowthVolumetric()
    : Growth(),
      tr_mandel_e_(Teuchos::null),
      lambda_fib_e_(Teuchos::null),
      growthtrig_const_(0.0),
      params_volumetric_(nullptr),
      refdir_(true),
      curdir_(Teuchos::null),
      curdir_for_update_(Teuchos::null),
      f_g_hist_(Teuchos::null)
{
}

/*----------------------------------------------------------------------------*/
Mat::GrowthVolumetric::GrowthVolumetric(Mat::PAR::Growth* params)
    : Growth(params),
      tr_mandel_e_(Teuchos::null),
      lambda_fib_e_(Teuchos::null),
      growthtrig_const_(0.0),
      params_volumetric_(params),
      refdir_(true),
      curdir_(Teuchos::null),
      curdir_for_update_(Teuchos::null),
      f_g_hist_(Teuchos::null)
{
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthVolumetric::evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 1>* stress, Core::LinAlg::Matrix<6, 6>* cmat, const int gp,
    const int eleGID)
{
  double time = params.get<double>("total time", -1.0);
  if (abs(time + 1.0) < 1e-14)
    FOUR_C_THROW("no time step or no total time given for growth material!");
  std::string action = params.get<std::string>("action", "none");
  bool output = false;
  if (action == "calc_struct_stress") output = true;

  const double eps = 1.0e-14;
  Mat::PAR::Growth* growth_params = Parameter();
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
      case Core::Materials::m_growth_ac_radial:
      case Core::Materials::m_growth_ac_radial_refconc:
      {
        // directional stuff......
        // push-forward of refdir
        Core::LinAlg::Matrix<3, 3> defgrdinv(true);
        defgrdinv.Invert(*defgrd);
        Core::LinAlg::Matrix<3, 1> curdir_for_update(true);
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
    Core::LinAlg::Matrix<6, 1> dthetadCvec(true);
    EvaluateGrowth(&theta, &dthetadCvec, defgrd, glstrain, params, gp, eleGID);

    // modify the parameter list to be passed to the elastic material
    Teuchos::ParameterList paramselast(params);
    paramselast.remove("matparderiv", false);

    Core::LinAlg::Matrix<6, 1> S(true);
    Core::LinAlg::Matrix<6, 6> cmatdach(true);

    GetSAndCmatdach(theta, defgrd, &S, &cmatdach, paramselast, gp, eleGID);

    *stress = S;

    // calculate growth part F_g of the deformation gradient F
    Core::LinAlg::Matrix<3, 3> F_g(true);

    Parameter()->growthlaw_->CalcFg(
        theta, ThetaOld()->at(gp), gp, defgrd, refdir_, curdir_, f_g_hist_, F_g);

    // calculate F_g^(-1)
    Core::LinAlg::Matrix<3, 3> F_ginv(true);
    F_ginv.Invert(F_g);

    // constitutive matrix including growth cmat = F_g^-1 F_g^-1 cmatdach F_g^-T F_g^-T
    Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatelast(true);
    cmatelast = Mat::pull_back_four_tensor<3>(F_ginv, cmatdach);

    //--------------------------------------------------------------------------------------
    // call material law with elastic part of defgr and elastic part of glstrain
    //--------------------------------------------------------------------------------------
    // build identity tensor I
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> Id(true);
    for (int i = 0; i < 3; i++) Id(i) = 1.0;

    // right Cauchy-Green Tensor  C = 2 * E + I
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> C(*glstrain);
    C.Scale(2.0);
    C += Id;

    // NOTE: we do this by a FD approximation, which is really cheap here due to the fact
    // that theta is a scalar value (hence only one more material evaluation is needed!)

    const double espilon = 1.0e-8;

    Core::LinAlg::Matrix<6, 1> SEps(true);
    Core::LinAlg::Matrix<6, 6> cmatdachEps(true);

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
      Core::LinAlg::Matrix<NUM_STRESS_3D, 1> cmatelasC(true);
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

    Core::LinAlg::Matrix<6, 1> Svec(true);
    Core::LinAlg::Matrix<6, 6> cmatdach(true);

    GetSAndCmatdach(theta, defgrd, &Svec, &cmatdach, params, gp, eleGID);

    *stress = Svec;

    // calculate growth part F_g of the deformation gradient F
    Core::LinAlg::Matrix<3, 3> F_g(true);

    Parameter()->growthlaw_->CalcFg(
        theta, ThetaOldAtGp(gp), gp, defgrd, refdir_, curdir_, f_g_hist_, F_g);

    // calculate F_g^(-1)
    Core::LinAlg::Matrix<3, 3> F_ginv(true);
    F_ginv.Invert(F_g);

    // constitutive matrix including growth cmat = F_g^-1 F_g^-1 cmatdach F_g^-T F_g^-T
    Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatelast(true);
    cmatelast = Mat::pull_back_four_tensor<3>(F_ginv, cmatdach);

    *cmat = cmatelast;
  }
  else
  {
    EvaluateElastic(defgrd, glstrain, stress, cmat, params, gp, eleGID);
    // build identity tensor I
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> Id(true);
    for (int i = 0; i < 3; i++) Id(i) = 1.0;
    // right Cauchy-Green Tensor  C = 2 * E + I
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> Cvec(*glstrain);
    Cvec.Scale(2.0);
    Cvec += Id;
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> Svec(true);
    Svec = *stress;

    tr_mandel_e_->at(gp) = Cvec(0) * Svec(0) + Cvec(1) * Svec(1) + Cvec(2) * Svec(2) +
                           Cvec(3) * Svec(3) + Cvec(4) * Svec(4) + Cvec(5) * Svec(5);

    // elastic fiber stretch
    Core::LinAlg::Matrix<3, 3> C(true);
    Core::LinAlg::Voigt::Strains::vector_to_matrix(Cvec, C);

    Core::LinAlg::Matrix<3, 1> CDir(true);
    CDir.MultiplyNN(1.0, C, refdir_);
    lambda_fib_e_->at(gp) =
        sqrt(CDir(0) * refdir_(0) + CDir(1) * refdir_(1) + CDir(2) * refdir_(2));
  }
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthVolumetric::EvaluateGrowth(double* theta, Core::LinAlg::Matrix<6, 1>* dthetadC,
    const Core::LinAlg::Matrix<3, 3>* defgrd, const Core::LinAlg::Matrix<6, 1>* glstrain,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  // get gauss point number
  if (gp == -1) FOUR_C_THROW("No Gauss point number provided in material.");

  double thetaold = ThetaOldAtGp(gp);

  Mat::Growth* matgrowth = this;
  Parameter()->growthlaw_->evaluate(theta, thetaold, dthetadC, *matgrowth, defgrd, glstrain,
      refdir_, curdir_, f_g_hist_, growthtrig_const_, params, gp, eleGID);
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthVolumetric::EvaluateNonLinMass(const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1>* linmass_disp,
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1>* linmass_vel, const int gp, const int eleGID)
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

    Mat::Growth* matgrowth = this;
    Parameter()->growthlaw_->evaluate(&theta, thetaold, linmass_disp, *matgrowth, defgrd, glstrain,
        refdir_, curdir_, f_g_hist_, growthtrig_const_, params, gp, eleGID);

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
void Mat::GrowthVolumetric::GetSAndCmatdach(const double theta,
    const Core::LinAlg::Matrix<3, 3>* defgrd, Core::LinAlg::Matrix<6, 1>* stress,
    Core::LinAlg::Matrix<6, 6>* cmatdach, Teuchos::ParameterList& params, const int gp,
    const int eleGID)
{
  // calculate growth part F_g of the deformation gradient F
  Core::LinAlg::Matrix<3, 3> F_g(true);
  Parameter()->growthlaw_->CalcFg(
      theta, ThetaOldAtGp(gp), gp, defgrd, refdir_, curdir_, f_g_hist_, F_g);

  // calculate F_g^(-1)
  Core::LinAlg::Matrix<3, 3> F_ginv(true);
  F_ginv.Invert(F_g);

  // elastic deformation gradient F_e = F * F_g^(-1)
  Core::LinAlg::Matrix<3, 3> defgrddach(true);  //*defgrd);
  defgrddach.MultiplyNN(*defgrd, F_ginv);       // Scale(1.0 / theta);

  // elastic right Cauchy-Green Tensor Cdach = F_e^T * F_e (= F_g^-T C F_g^-1)
  Core::LinAlg::Matrix<3, 3> Cdach(true);
  Cdach.MultiplyTN(defgrddach, defgrddach);

  // transform Cdach into a vector
  Core::LinAlg::Matrix<6, 1> Cdachvec(true);
  Core::LinAlg::Voigt::Strains::matrix_to_vector(Cdach, Cdachvec);

  //--------------------------------------------------------------------------------------
  // call material law with elastic part of defgr and elastic part of glstrain
  //--------------------------------------------------------------------------------------
  // build identity tensor I
  Core::LinAlg::Matrix<6, 1> Id(true);
  for (int i = 0; i < 3; i++) Id(i) = 1.0;

  Core::LinAlg::Matrix<6, 1> glstraindachvec(Cdachvec);
  glstraindachvec -= Id;
  glstraindachvec.Scale(0.5);

  Core::LinAlg::Matrix<6, 1> Sdachvec(true);
  // elastic 2 PK stress and constitutive matrix
  Matelastic()->evaluate(&defgrddach, &glstraindachvec, params, &Sdachvec, cmatdach, gp, eleGID);

  // calculate stress
  // 2PK stress S = F_g^-1 Sdach F_g^-T
  Core::LinAlg::Matrix<3, 3> Sdach(true);
  Core::LinAlg::Voigt::Stresses::vector_to_matrix(Sdachvec, Sdach);

  Core::LinAlg::Matrix<3, 3> tmp(true);
  tmp.MultiplyNT(Sdach, F_ginv);
  Core::LinAlg::Matrix<3, 3> S(true);
  S.MultiplyNN(F_ginv, tmp);

  Core::LinAlg::Voigt::Stresses::matrix_to_vector(S, *stress);

  // trace of elastic Mandel stress Mdach = Cdach Sdach
  tr_mandel_e_->at(gp) = Cdachvec(0) * Sdachvec(0) + Cdachvec(1) * Sdachvec(1) +
                         Cdachvec(2) * Sdachvec(2) + Cdachvec(3) * Sdachvec(3) +
                         Cdachvec(4) * Sdachvec(4) + Cdachvec(5) * Sdachvec(5);

  // elastic fiber stretch lambda = \sqrt(f_0 \cdot Cdach f_0)
  Core::LinAlg::Matrix<3, 1> CdachDir(true);
  CdachDir.MultiplyNN(1.0, Cdach, refdir_);
  lambda_fib_e_->at(gp) =
      sqrt(CdachDir(0) * refdir_(0) + CdachDir(1) * refdir_(1) + CdachDir(2) * refdir_(2));
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthVolumetric::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  Mat::PAR::Growth* params = Parameter();
  if (params != nullptr) matid = params->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  int numgp = 0;
  if (isinit_)
  {
    numgp = theta_->size();  // size is number of Gauss points
  }
  add_to_pack(data, numgp);

  // Pack internal variables
  for (int gp = 0; gp < numgp; ++gp)
  {
    add_to_pack(data, tr_mandel_e_->at(gp));
  }
  for (int gp = 0; gp < numgp; ++gp)
  {
    add_to_pack(data, lambda_fib_e_->at(gp));
  }

  // Pack prescribed const growth trigger
  add_to_pack(data, growthtrig_const_);

  // Pack internal variables
  for (int i = 0; i < 3; ++i)
  {
    add_to_pack(data, refdir_(i, 0));
  }

  // int numgp=0;
  if (isinit_)
  {
    numgp = theta_->size();  // size is number of Fauss points
  }
  add_to_pack(data, numgp);

  for (int gp = 0; gp < numgp; gp++)
  {
    Core::LinAlg::Matrix<3, 3> F_g_hist = f_g_hist_.at(gp);
    Core::LinAlg::Matrix<3, 1> curdir = curdir_.at(gp);
    Core::LinAlg::Matrix<3, 1> curdir_for_update = curdir_for_update_.at(gp);

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        add_to_pack(data, F_g_hist(i, j));
      }
      add_to_pack(data, curdir(i, 0));
      add_to_pack(data, curdir_for_update(i, 0));
    }
  }

  // Pack base class material
  Growth::pack(data);
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthVolumetric::unpack(const std::vector<char>& data)
{
  isinit_ = true;
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  extract_from_pack(position, data, matid);

  params_volumetric_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_volumetric_ = dynamic_cast<Mat::PAR::Growth*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  int numgp;
  extract_from_pack(position, data, numgp);
  if (numgp == 0)
  {  // no history data to unpack
    isinit_ = false;
    if (position != data.size())
      FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
    return;
  }

  // unpack growth internal variables
  tr_mandel_e_ = Teuchos::rcp(new std::vector<double>(numgp));
  for (int gp = 0; gp < numgp; ++gp)
  {
    double a;
    extract_from_pack(position, data, a);
    tr_mandel_e_->at(gp) = a;
  }

  lambda_fib_e_ = Teuchos::rcp(new std::vector<double>(numgp));
  for (int gp = 0; gp < numgp; ++gp)
  {
    double a;
    extract_from_pack(position, data, a);
    lambda_fib_e_->at(gp) = a;
  }

  double growthtrigconst;
  extract_from_pack(position, data, growthtrigconst);
  growthtrig_const_ = growthtrigconst;

  // unpack growth internal variables
  // Pack internal variables
  for (int i = 0; i < 3; ++i)
  {
    double refdir_i;
    extract_from_pack(position, data, refdir_i);
    refdir_(i, 0) = refdir_i;
  }

  // int numgp;
  extract_from_pack(position, data, numgp);
  if (numgp != 0)
  {
    f_g_hist_ = std::vector<Core::LinAlg::Matrix<3, 3>>(numgp, Core::LinAlg::Matrix<3, 3>(true));
    curdir_ = std::vector<Core::LinAlg::Matrix<3, 1>>(numgp, Core::LinAlg::Matrix<3, 1>(true));
    curdir_for_update_ =
        std::vector<Core::LinAlg::Matrix<3, 1>>(numgp, Core::LinAlg::Matrix<3, 1>(true));

    for (int gp = 0; gp < numgp; gp++)
    {
      Core::LinAlg::Matrix<3, 3> F_g_hist(true);
      Core::LinAlg::Matrix<3, 1> curdir(true);
      Core::LinAlg::Matrix<3, 1> curdir_for_update(true);

      for (int i = 0; i < 3; ++i)
      {
        for (int j = 0; j < 3; ++j)
        {
          double F_g_hist_ij;
          extract_from_pack(position, data, F_g_hist_ij);
          F_g_hist(i, j) = F_g_hist_ij;
        }
        double curdir_i;
        extract_from_pack(position, data, curdir_i);
        curdir(i, 0) = curdir_i;

        double curdir_for_update_i;
        extract_from_pack(position, data, curdir_for_update_i);
        curdir_for_update(i, 0) = curdir_for_update_i;
      }
      f_g_hist_.at(gp) = F_g_hist;
      curdir_.at(gp) = curdir;
      curdir_for_update_.at(gp) = curdir_for_update;
    }
  }

  // extract base class material
  std::vector<char> basedata(0);
  Growth::extract_from_pack(position, data, basedata);
  Growth::unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------------*/
void Mat::GrowthVolumetric::setup(int numgp, Input::LineDefinition* linedef)
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
    case Core::Materials::m_growth_ac_radial:
    case Core::Materials::m_growth_ac_radial_refconc:
    {
      // CIR-AXI-RAD nomenclature
      if (not(linedef->has_named("RAD")))
      {
        FOUR_C_THROW(
            "If you want growth into the radial direction you need to specify RAD in your input "
            "file!");
      }

      ReadFiber(linedef, "RAD", refdir_);
      curdir_ = std::vector<Core::LinAlg::Matrix<3, 1>>(numgp, refdir_);
      curdir_for_update_ = std::vector<Core::LinAlg::Matrix<3, 1>>(numgp, refdir_);

      Core::LinAlg::Matrix<3, 3> Id(true);
      for (int i = 0; i < 3; i++) Id(i, i) = 1.0;
      f_g_hist_ = std::vector<Core::LinAlg::Matrix<3, 3>>(numgp, Id);
    }
    break;
    case Core::Materials::m_growth_aniso_strain:
    case Core::Materials::m_growth_aniso_stress:
    {
      // FIBER1 nomenclature
      if (not(linedef->has_named("FIBER1")))
        FOUR_C_THROW(
            "If you want growth in fiber direction you need to specify FIBER1 in your input file!");

      ReadFiber(linedef, "FIBER1", refdir_);

      // only refdir is used - rest remains unused...
      curdir_ = std::vector<Core::LinAlg::Matrix<3, 1>>(numgp, refdir_);
      curdir_for_update_ = std::vector<Core::LinAlg::Matrix<3, 1>>(numgp, refdir_);
      Core::LinAlg::Matrix<3, 3> Id(true);
      for (int i = 0; i < 3; i++) Id(i, i) = 1.0;
      f_g_hist_ = std::vector<Core::LinAlg::Matrix<3, 3>>(numgp, Id);
    }
    break;
    case Core::Materials::m_growth_aniso_strain_const_trig:
    case Core::Materials::m_growth_aniso_stress_const_trig:
    {
      // FIBER1 nomenclature
      if (not(linedef->has_named("FIBER1")))
        FOUR_C_THROW(
            "If you want growth in fiber direction you need to specify FIBER1 in your input file!");

      ReadFiber(linedef, "FIBER1", refdir_);

      linedef->extract_double("GROWTHTRIG", growthtrig_const_);
      if (not(linedef->has_named("GROWTHTRIG")))
        FOUR_C_THROW("You need to specify GROWTHTRIG in your input file!");

      // only refdir is used - rest remains unused...
      curdir_ = std::vector<Core::LinAlg::Matrix<3, 1>>(numgp, refdir_);
      curdir_for_update_ = std::vector<Core::LinAlg::Matrix<3, 1>>(numgp, refdir_);
      Core::LinAlg::Matrix<3, 3> Id(true);
      for (int i = 0; i < 3; i++) Id(i, i) = 1.0;
      f_g_hist_ = std::vector<Core::LinAlg::Matrix<3, 3>>(numgp, Id);
    }
    break;
    default:
    {
      // directions are unused
      refdir_(true);
      curdir_ = std::vector<Core::LinAlg::Matrix<3, 1>>(numgp, refdir_);
      curdir_for_update_ = std::vector<Core::LinAlg::Matrix<3, 1>>(numgp, refdir_);
      Core::LinAlg::Matrix<3, 3> Id(true);
      for (int i = 0; i < 3; i++) Id(i, i) = 1.0;
      f_g_hist_ = std::vector<Core::LinAlg::Matrix<3, 3>>(numgp, Id);
    }
    break;
  }

  // setup base class
  Growth::setup(numgp, linedef);
}



/*----------------------------------------------------------------------------*/
void Mat::GrowthVolumetric::update()
{
  const int numgp = theta_->size();

  // setup anisotropic growth laws
  switch (Parameter()->growthlaw_->MaterialType())
  {
    case Core::Materials::m_growth_ac_radial:
    case Core::Materials::m_growth_ac_radial_refconc:
    {
      const Core::LinAlg::Matrix<3, 3> dummydefgrad(true);

      for (int gp = 0; gp < numgp; gp++)
      {
        Core::LinAlg::Matrix<3, 3> F_g_hist_new(true);

        Parameter()->growthlaw_->CalcFg(theta_->at(gp), ThetaOld()->at(gp), gp, &dummydefgrad,
            refdir_, curdir_, f_g_hist_, F_g_hist_new);

        f_g_hist_.at(gp) = F_g_hist_new;

        curdir_.at(gp) = curdir_for_update_.at(gp);
      }
    }
    break;
    default:
      break;
  }

  // update base class
  Growth::update();
}



/*----------------------------------------------------------------------*
 | Function which reads in the given fiber value             Thon 01/15 |
 *----------------------------------------------------------------------*/
void Mat::GrowthVolumetric::ReadFiber(
    Input::LineDefinition* linedef, std::string specifier, Core::LinAlg::Matrix<3, 1>& fiber_vector)
{
  std::vector<double> fiber1;
  linedef->extract_double_vector(std::move(specifier), fiber1);
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

FOUR_C_NAMESPACE_CLOSE
