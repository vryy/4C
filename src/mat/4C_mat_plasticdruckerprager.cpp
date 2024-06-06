/*----------------------------------------------------------------------*/
/*! \file
\brief Contains the functions to establish local material law stress-strain law for
isotropic material for a 3D hex element following Drucker Prager plasticity model
Reference:
EA de Souza Neto, D Peric, DRJ Owen. Computational Methods of Plasticity: Theory and Applications,
John Wiley & Sons, Ltd, 2008
\level 3
*/
/*----------------------------------------------------------------------*/
#include "4C_mat_plasticdruckerprager.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_FADmatrix_utils.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_stvenantkirchhoff.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_utils_local_newton.hpp"

FOUR_C_NAMESPACE_OPEN


Mat::PAR::PlasticDruckerPrager::PlasticDruckerPrager(Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : Parameter(matdata),
      youngs_(matdata->Get<double>("YOUNG")),
      poissonratio_(matdata->Get<double>("NUE")),
      density_(matdata->Get<double>("DENS")),
      isohard_(matdata->Get<double>("ISOHARD")),
      abstol_(matdata->Get<double>("TOL")),
      cohesion_(matdata->Get<double>("C")),
      eta_(matdata->Get<double>("ETA")),
      xi_(matdata->Get<double>("XI")),
      etabar_(matdata->Get<double>("ETABAR")),
      itermax_(matdata->Get<int>("MAXITER"))
{
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::PlasticDruckerPrager::create_material()
{
  return Teuchos::rcp(new Mat::PlasticDruckerPrager(this));
}
Mat::PlasticDruckerPragerType Mat::PlasticDruckerPragerType::instance_;

Core::Communication::ParObject* Mat::PlasticDruckerPragerType::Create(const std::vector<char>& data)
{
  Mat::PlasticDruckerPrager* plastic = new Mat::PlasticDruckerPrager();
  plastic->Unpack(data);
  return plastic;
}

Mat::PlasticDruckerPrager::PlasticDruckerPrager() : params_(nullptr) {}

Mat::PlasticDruckerPrager::PlasticDruckerPrager(Mat::PAR::PlasticDruckerPrager* params)
    : params_(params)
{
}

void Mat::PlasticDruckerPrager::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();
  AddtoPack(data, matid);
  int histsize = Initialized() ? strainpllast_.size() : 0;
  AddtoPack(data, histsize);
  for (int var = 0; var < histsize; ++var)
  {
    AddtoPack(data, strainpllast_.at(var));
    AddtoPack(data, strainbarpllast_.at(var));
  }
}

void Mat::PlasticDruckerPrager::Unpack(const std::vector<char>& data)
{
  isinit_ = true;
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::PlasticDruckerPrager*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

    int histsize;
    ExtractfromPack(position, data, histsize);

    if (histsize == 0) isinit_ = false;

    strainpllast_ = std::vector<Core::LinAlg::Matrix<NUM_STRESS_3D, 1>>();
    strainplcurr_ = std::vector<Core::LinAlg::Matrix<NUM_STRESS_3D, 1>>();
    strainbarpllast_ = std::vector<double>();
    strainbarplcurr_ = std::vector<double>();
    for (int var = 0; var < histsize; ++var)
    {
      Core::LinAlg::Matrix<NUM_STRESS_3D, 1> tmp_vect(true);
      double tmp_scalar = 0.0;

      ExtractfromPack(position, data, tmp_vect);
      strainpllast_.push_back(tmp_vect);

      ExtractfromPack(position, data, tmp_scalar);
      strainbarpllast_.push_back(tmp_scalar);

      strainplcurr_.push_back(tmp_vect);
      strainbarplcurr_.push_back(tmp_scalar);
    }
    if (position != data.size())
      FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
  }
}
void Mat::PlasticDruckerPrager::Setup(int numgp, Input::LineDefinition* linedef)
{
  strainpllast_.resize(numgp);
  strainplcurr_.resize(numgp);

  strainbarpllast_.resize(numgp);
  strainbarplcurr_.resize(numgp);

  isinit_ = true;
}

void Mat::PlasticDruckerPrager::Update()
{
  strainpllast_ = strainplcurr_;
  strainbarpllast_ = strainbarplcurr_;

  std::for_each(strainplcurr_.begin(), strainplcurr_.end(), [](auto& item) { item.Clear(); });
  std::fill(strainbarplcurr_.begin(), strainbarplcurr_.end(), 0.0);
}

void Mat::PlasticDruckerPrager::setup_cmat(Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat)
{
  double young = params_->youngs_;

  double nu = params_->poissonratio_;

  Mat::StVenantKirchhoff::FillCmat(cmat, young, nu);
}

void Mat::PlasticDruckerPrager::setup_cmat_elasto_plastic_cone(
    Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat, double Dgamma, double G, double Kappa,
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& devstrain, double xi, double Hiso, double eta,
    double etabar)
{
  cmat.Clear();

  Core::LinAlg::Matrix<NUM_STRESS_3D, 1> id2(true);
  Core::LinAlg::Voigt::IdentityMatrix(id2);
  Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> id4sharp(true);
  for (int i = 0; i < 3; i++) id4sharp(i, i) = 1.0;
  for (int i = 3; i < 6; i++) id4sharp(i, i) = 0.5;

  const double normdevstrain =
      sqrt(devstrain(0) * devstrain(0) + devstrain(1) * devstrain(1) + devstrain(2) * devstrain(2) +
           2 * (devstrain(3) * devstrain(3) + devstrain(4) * devstrain(4) +
                   devstrain(5) * devstrain(5)));
  const double epfac = 2 * G * (1 - (Dgamma / sqrt(2) / normdevstrain));

  cmat.Update(epfac, id4sharp, 1.0);

  double epfac1 = 0.0;
  double epfac2 = 0.0;
  double epfac3 = 0.0;
  double epfac4 = 0.0;
  epfac1 = epfac / (-3.0);
  cmat.MultiplyNT(epfac1, id2, id2, 1.0);

  double A = 0.0;
  A = 1 / (G + Kappa * etabar * eta + xi * xi * Hiso);
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1> D(true);
  D.Update(1 / normdevstrain, devstrain);
  epfac2 = 2 * G * (Dgamma / (sqrt(2) * normdevstrain) - G * A);

  for (int k = 0; k < 6; ++k)
  {
    for (int i = 0; i < 6; ++i)
    {
      cmat(i, k) += epfac2 * D(i) * D(k);
    }
  }

  epfac3 = -sqrt(2) * G * A * Kappa;

  for (int k = 0; k < 6; ++k)
  {
    for (int i = 0; i < 6; ++i)
    {
      cmat(k, i) += epfac3 * (eta * D(k) * id2(i) + etabar * id2(k) * D(i));
    }
  }

  epfac4 = Kappa * (1 - Kappa * eta * etabar * A);

  for (int k = 0; k < 6; ++k)
  {
    for (int i = 0; i < 6; ++i)
    {
      cmat(i, k) += epfac4 * id2(k) * id2(i);
    }
  }
}

void Mat::PlasticDruckerPrager::setup_cmat_elasto_plastic_apex(
    Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat, double Kappa,
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& devstrain, double xi, double Hiso, double eta,
    double etabar)
{
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1> id2(true);
  for (int i = 0; i < 3; i++) id2(i) = 1.0;
  double epfac = 0.0;
  epfac = Kappa * (1 - Kappa / (Kappa + xi / eta * xi / etabar * Hiso));
  cmat.Clear();
  cmat.MultiplyNT(epfac, id2, id2, 0.0);
}

template <typename ScalarT>
void Mat::PlasticDruckerPrager::EvaluateFAD(const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<6, 1, ScalarT>* linstrain, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 1, ScalarT>* stress, Core::LinAlg::Matrix<6, 6>* cmat, const int gp,
    const int eleGID)
{
  Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1> plstrain(true);

  ScalarT young = params_->youngs_;
  ScalarT nu = params_->poissonratio_;
  ScalarT Hiso = params_->isohard_;
  ScalarT cohesion = params_->cohesion_;
  ScalarT eta = params_->eta_;
  ScalarT xi = params_->xi_;
  ScalarT etabar = params_->etabar_;
  const int itermax = params_->itermax_;
  ScalarT G = 0.0;
  G = young / (2.0 * (1.0 + nu));
  ScalarT kappa = 0.0;
  kappa = young / (3.0 * (1.0 - 2.0 * nu));
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1> id2(true);
  for (int i = 0; i < 3; i++) id2(i) = 1.0;
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT> strain(*linstrain);

  Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT> strain_p(false);
  for (int i = 0; i < 6; i++) strain_p(i, 0) = strainpllast_.at(gp)(i, 0);
  ScalarT strainbar_p = 0.0;
  strainbar_p = (strainbarpllast_.at(gp));
  if (strainbarpllast_.at(gp) < 0.0)
    FOUR_C_THROW("accumulated plastic strain has to be equal to or greater than zero!");

  for (int i = 3; i < 6; ++i) strain(i) /= 2.0;
  for (int i = 3; i < 6; ++i) strain_p(i) /= 2.0;

  Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT> strain_e(true);
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT> trialstrain_e(false);
  trialstrain_e.Update(1.0, strain, (-1.0), strain_p);
  ScalarT tracestrain = trialstrain_e(0) + trialstrain_e(1) + trialstrain_e(2);
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT> volumetricstrain(false);
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT> id2Scalar(true);
  for (int i = 0; i < NUM_STRESS_3D; ++i) id2Scalar(i) = static_cast<ScalarT>(id2(i));
  volumetricstrain.Update((tracestrain / 3.0), id2Scalar);
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT> devstrain(false);
  devstrain.Update(1.0, trialstrain_e, (-1.0), volumetricstrain);

  ScalarT p = kappa * tracestrain;
  ScalarT p_trial = p;
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT> devstress(false);
  devstress.Update((2.0 * G), devstrain);

  ScalarT J2 = 1.0 / 2.0 *
                   (devstress(0) * devstress(0) + devstress(1) * devstress(1) +
                       devstress(2) * devstress(2)) +
               devstress(3) * devstress(3) + devstress(4) * devstress(4) +
               devstress(5) * devstress(5);
  ScalarT Phi_trial = 0.0;
  Phi_trial = std::sqrt(J2) + eta * p - xi * cohesion - xi * Hiso * strainbar_p;
  ScalarT Dgamma = 0.0;
  ScalarT dstrainv = 0.0;
  if (Phi_trial / abs(cohesion) > params_->abstol_)
  {
    auto returnToConeFunctAndDeriv = [this, &G, &kappa, &Phi_trial](ScalarT Dgamma_init)
    { return this->return_to_cone_funct_and_deriv(Dgamma_init, G, kappa, Phi_trial); };

    const double tol = params_->abstol_;
    Dgamma =
        Core::UTILS::SolveLocalNewton(returnToConeFunctAndDeriv, Dgamma, tol * cohesion, itermax);
    strainbar_p = (strainbarpllast_.at(gp)) + xi * Dgamma;
    devstress.Scale(1.0 - (G * Dgamma / std::sqrt(J2)));
    p = p_trial - kappa * etabar * Dgamma;
    if ((std::sqrt(J2) - G * Dgamma) / abs(cohesion) < params_->abstol_)
    {
      strainbar_p = (strainbarpllast_.at(gp));
      auto returnToApexFunctAndDeriv = [this, &p_trial, &kappa, &strainbar_p](ScalarT dstrainv_init)
      { return this->return_to_apex_funct_and_deriv(dstrainv_init, p_trial, kappa, strainbar_p); };

      const double tol = params_->abstol_;
      dstrainv = Core::UTILS::SolveLocalNewton(
          returnToApexFunctAndDeriv, dstrainv, tol * cohesion, itermax);
      strainbar_p = (strainbarpllast_.at(gp)) + xi / eta * dstrainv;
      p = p_trial - kappa * dstrainv;
      for (int i = 0; i < 6; i++) devstress(i) = 0.0;
    }
    Stress(p, devstress, *stress);
    strain_e.Update(1 / G / 2, devstress, p / kappa / 3, id2Scalar);
    for (int i = 3; i < 6; ++i) strain_e(i) *= 2.0;
    for (int i = 3; i < 6; ++i) strain(i) *= 2.0;
    strain_p.Update(1.0, strain, -1.0, strain_e);

    strainplcurr_.at(gp) = Core::FADUtils::CastToDouble(strain_p);
    strainbarplcurr_.at(gp) = Core::FADUtils::CastToDouble(strainbar_p);
  }
  else
  {
    Stress(p, devstress, *stress);
    strain_e.Update(trialstrain_e);
    for (int i = 3; i < 6; ++i) strain_e(i) *= 2.0;
    for (int i = 3; i < 6; ++i) strain(i) *= 2.0;

    strainplcurr_.at(gp) = strainpllast_.at(gp);
    strainbarplcurr_.at(gp) = strainbarpllast_.at(gp);
  }
  if (Phi_trial > 0)
  {
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> devstraindouble =
        Core::FADUtils::CastToDouble(devstrain);
    if (dstrainv != 0.0)
    {
      setup_cmat_elasto_plastic_apex(*cmat, Core::FADUtils::CastToDouble(kappa), devstraindouble,
          Core::FADUtils::CastToDouble(xi), Core::FADUtils::CastToDouble(Hiso),
          Core::FADUtils::CastToDouble(eta), Core::FADUtils::CastToDouble(etabar));
    }
    else
    {
      setup_cmat_elasto_plastic_cone(*cmat, Core::FADUtils::CastToDouble(Dgamma),
          Core::FADUtils::CastToDouble(G), Core::FADUtils::CastToDouble(kappa), devstraindouble,
          Core::FADUtils::CastToDouble(xi), Core::FADUtils::CastToDouble(Hiso),
          Core::FADUtils::CastToDouble(eta), Core::FADUtils::CastToDouble(etabar));
    }
  }
  else
  {
    setup_cmat(*cmat);
  }
}

void Mat::PlasticDruckerPrager::register_output_data_names(
    std::unordered_map<std::string, int>& names_and_size) const
{
  names_and_size["accumulated_plastic_strain"] = 1;
  names_and_size["plastic_strain"] = 6;
}

bool Mat::PlasticDruckerPrager::EvaluateOutputData(
    const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const
{
  if (name == "accumulated_plastic_strain")
  {
    for (std::size_t gp = 0; gp < strainbarplcurr_.size(); ++gp)
    {
      data(gp, 0) = strainbarplcurr_.at(int(gp));
    }
    return true;
  }
  if (name == "plastic_strain")
  {
    for (std::size_t gp = 0; gp < strainplcurr_.size(); ++gp)
    {
      const double* values = strainplcurr_.at(gp).A();
      for (std::size_t i = 0; i < 6; ++i)
      {
        data(gp, i) = values[i];
      }
    }
    return true;
  }
  return false;
}

template <typename T>
void Mat::PlasticDruckerPrager::Stress(const T p,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1, T>& devstress,
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1, T>& stress)
{
  stress.Update(devstress);
  for (int i = 0; i < 3; ++i) stress(i) += p;
}

template <typename T>
std::pair<T, T> Mat::PlasticDruckerPrager::return_to_cone_funct_and_deriv(
    T Dgamma, T G, T kappa, T Phi_trial)
{
  double Hiso = params_->isohard_;
  double eta = params_->eta_;
  double xi = params_->xi_;
  double etabar = params_->etabar_;
  T Res = Phi_trial - Dgamma * (G + eta * kappa * etabar) - (xi * xi * Dgamma * Hiso);
  T d = -G - (kappa * etabar * eta) - (xi * xi * Hiso);
  return {Res, d};
}

template <typename T>
std::pair<T, T> Mat::PlasticDruckerPrager::return_to_apex_funct_and_deriv(
    T dstrainv, T p, T kappa, T strainbar_p)
{
  double Hiso = params_->isohard_;
  double eta = params_->eta_;
  double xi = params_->xi_;
  double cohesion = params_->cohesion_;
  double etabar = params_->etabar_;
  double alpha = xi / eta;
  double beta = xi / etabar;
  T Res =
      beta * cohesion + beta * strainbar_p * Hiso - p + dstrainv * (alpha * beta * Hiso + kappa);
  T d = xi * xi / eta / etabar * Hiso + kappa;
  return {Res, d};
}

template void Mat::PlasticDruckerPrager::EvaluateFAD<double>(const Core::LinAlg::Matrix<3, 3>*,
    const Core::LinAlg::Matrix<6, 1, double>*, Teuchos::ParameterList&,
    Core::LinAlg::Matrix<6, 1, double>*, Core::LinAlg::Matrix<6, 6>*, int gp, int eleGID);
template void Mat::PlasticDruckerPrager::EvaluateFAD<FAD>(const Core::LinAlg::Matrix<3, 3>*,
    const Core::LinAlg::Matrix<6, 1, FAD>*, Teuchos::ParameterList&,
    Core::LinAlg::Matrix<6, 1, FAD>*, Core::LinAlg::Matrix<6, 6>*, int gp, int eleGID);
template void Mat::PlasticDruckerPrager::Stress<double>(const double p,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1, double>& devstress,
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1, double>& stress);
template void Mat::PlasticDruckerPrager::Stress<FAD>(const FAD p,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1, FAD>& devstress,
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1, FAD>& stress);
template std::pair<double, double>
Mat::PlasticDruckerPrager::return_to_cone_funct_and_deriv<double>(
    double Dgamma, double G, double kappa, double Phi_trial);
template std::pair<FAD, FAD> Mat::PlasticDruckerPrager::return_to_cone_funct_and_deriv<FAD>(
    FAD Dgamma, FAD G, FAD kappa, FAD Phi_trial);
template std::pair<double, double>
Mat::PlasticDruckerPrager::return_to_apex_funct_and_deriv<double>(
    double dstrainv, double p, double kappa, double strainbar_p);
template std::pair<FAD, FAD> Mat::PlasticDruckerPrager::return_to_apex_funct_and_deriv<FAD>(
    FAD dstrainv, FAD p, FAD kappa, FAD strainbar_p);

FOUR_C_NAMESPACE_CLOSE
