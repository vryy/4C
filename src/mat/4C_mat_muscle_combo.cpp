/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of the Combo active skeletal muscle material with variable time-dependent
activations

\level 3

*/
/*----------------------------------------------------------------------*/

#include "4C_mat_muscle_combo.hpp"

#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_mat_muscle_utils.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"
#include "4C_matelast_aniso_structuraltensor_strategy.hpp"
#include "4C_utils_exceptions.hpp"

#include <variant>

FOUR_C_NAMESPACE_OPEN

namespace
{
  using ActivationMapType = std::unordered_map<int, std::vector<std::pair<double, double>>>;

  Mat::PAR::MuscleCombo::ActivationParameterVariant GetActivationParams(
      const Core::Mat::PAR::Parameter::Data& matdata,
      const Inpar::Mat::ActivationType& activation_type)
  {
    if (activation_type == Inpar::Mat::ActivationType::function_of_space_time)
    {
      auto actFunctId = matdata.parameters.get<int>("FUNCTID");
      if (actFunctId <= 0) FOUR_C_THROW("Function id must be positive");
      return actFunctId;
    }
    else if (activation_type == Inpar::Mat::ActivationType::map)
    {
      return matdata.parameters.get<const ActivationMapType>("MAPFILE");
    }
    else
      return std::monostate{};
  }

  struct ActivationParamsVisitor
  {
    Mat::MuscleCombo::ActivationEvaluatorVariant operator()(const int function_id) const
    {
      return &Global::Problem::Instance()->FunctionById<Core::UTILS::FunctionOfSpaceTime>(
          function_id - 1);
    }

    Mat::MuscleCombo::ActivationEvaluatorVariant operator()(const ActivationMapType& map) const
    {
      return &map;
    }

    Mat::MuscleCombo::ActivationEvaluatorVariant operator()(const std::monostate& /*unused*/) const
    {
      FOUR_C_THROW(
          "Error in ActivationParamsVisitor. You're calling it with the default-constructed "
          "state.");
    }
  };

  struct ActivationEvalVisitor
  {
    double operator()(const ActivationMapType*& map) const
    {
      // use one-based element ids in the pattern file (corresponding to the ones in the dat-file)
      return Mat::UTILS::Muscle::EvaluateTimeSpaceDependentActiveStressByMap(
          Popt_, *map, t_tot_, eleGID_ + 1);
    }

    double operator()(const Core::UTILS::FunctionOfSpaceTime*& function) const
    {
      return Mat::UTILS::Muscle::EvaluateTimeSpaceDependentActiveStressByFunct(
          Popt_, *function, t_tot_, element_center_reference_coordinates_);
    }

    double operator()(const std::monostate& /*unused*/) const
    {
      FOUR_C_THROW(
          "Error in ActivationEvalVisitor. You're calling it with the default-constructed state.");
    }

    const double& Popt_;
    const double& t_tot_;
    const Core::LinAlg::Matrix<3, 1>& element_center_reference_coordinates_;
    const int& eleGID_;
  };
}  // namespace


Mat::PAR::MuscleCombo::MuscleCombo(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      alpha_(matdata.parameters.get<double>("ALPHA")),
      beta_(matdata.parameters.get<double>("BETA")),
      gamma_(matdata.parameters.get<double>("GAMMA")),
      kappa_(matdata.parameters.get<double>("KAPPA")),
      omega0_(matdata.parameters.get<double>("OMEGA0")),
      Popt_(matdata.parameters.get<double>("POPT")),
      lambdaMin_(matdata.parameters.get<double>("LAMBDAMIN")),
      lambdaOpt_(matdata.parameters.get<double>("LAMBDAOPT")),
      activationType_(
          static_cast<Inpar::Mat::ActivationType>(matdata.parameters.get<int>("ACTEVALTYPE"))),
      activationParams_(GetActivationParams(matdata, activationType_)),
      density_(matdata.parameters.get<double>("DENS"))
{
  // error handling for parameter ranges
  // passive material parameters
  if (alpha_ <= 0.0) FOUR_C_THROW("Material parameter ALPHA must be greater zero");
  if (beta_ <= 0.0) FOUR_C_THROW("Material parameter BETA must be greater zero");
  if (gamma_ <= 0.0) FOUR_C_THROW("Material parameter GAMMA must be greater zero");
  if (omega0_ < 0.0 || omega0_ > 1.0) FOUR_C_THROW("Material parameter OMEGA0 must be in [0;1]");

  // active material parameters
  if (Popt_ < 0.0)
  {
    FOUR_C_THROW("Material parameter POPT must be postive or zero");
  }

  // stretch dependent parameters
  if (lambdaMin_ <= 0.0) FOUR_C_THROW("Material parameter LAMBDAMIN must be postive");
  if (lambdaOpt_ <= 0.0) FOUR_C_THROW("Material parameter LAMBDAOPT must be postive");

  // density
  if (density_ < 0.0) FOUR_C_THROW("DENS should be positive");
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::MuscleCombo::create_material()
{
  return Teuchos::rcp(new Mat::MuscleCombo(this));
}

Mat::MuscleComboType Mat::MuscleComboType::instance_;

Core::Communication::ParObject* Mat::MuscleComboType::Create(const std::vector<char>& data)
{
  auto* muscle_combo = new Mat::MuscleCombo();
  muscle_combo->unpack(data);
  return muscle_combo;
}

Mat::MuscleCombo::MuscleCombo()
    : params_(nullptr),
      anisotropy_(),
      anisotropy_extension_(true, 0.0, 0,
          Teuchos::rcp<Mat::Elastic::StructuralTensorStrategyBase>(
              new Mat::Elastic::StructuralTensorStrategyStandard(nullptr)),
          {0}),
      activation_evaluator_(std::monostate{})
{
}

Mat::MuscleCombo::MuscleCombo(Mat::PAR::MuscleCombo* params)
    : params_(params),
      anisotropy_(),
      anisotropy_extension_(true, 0.0, 0,
          Teuchos::rcp<Mat::Elastic::StructuralTensorStrategyBase>(
              new Mat::Elastic::StructuralTensorStrategyStandard(nullptr)),
          {0}),
      activation_evaluator_(std::monostate{})
{
  // register anisotropy extension to global anisotropy
  anisotropy_.register_anisotropy_extension(anisotropy_extension_);

  // initialize fiber directions and structural tensor
  anisotropy_extension_.register_needed_tensors(
      Mat::FiberAnisotropyExtension<1>::FIBER_VECTORS |
      Mat::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR |
      Mat::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR_STRESS);

  // cannot set activation_function here, because function manager did not yet read functions
}

void Mat::MuscleCombo::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  anisotropy_extension_.pack_anisotropy(data);
}

void Mat::MuscleCombo::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // make sure we have a pristine material
  params_ = nullptr;

  // matid and recover params_
  int matid;
  extract_from_pack(position, data, matid);

  if (Global::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
      {
        params_ = static_cast<Mat::PAR::MuscleCombo*>(mat);
        activation_evaluator_ = std::visit(ActivationParamsVisitor(), params_->activationParams_);
      }
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  anisotropy_extension_.unpack_anisotropy(data, position);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

void Mat::MuscleCombo::setup(int numgp, Input::LineDefinition* linedef)
{
  // Read anisotropy
  anisotropy_.set_number_of_gauss_points(numgp);
  anisotropy_.read_anisotropy_from_element(linedef);

  activation_evaluator_ = std::visit(ActivationParamsVisitor(), params_->activationParams_);
}

void Mat::MuscleCombo::Update(Core::LinAlg::Matrix<3, 3> const& defgrd, int const gp,
    Teuchos::ParameterList& params, int const eleGID)
{
}

void Mat::MuscleCombo::evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 1>* stress, Core::LinAlg::Matrix<6, 6>* cmat, const int gp,
    const int eleGID)
{
  Core::LinAlg::Matrix<6, 1> Sc_stress(true);
  Core::LinAlg::Matrix<6, 6> ccmat(true);

  // get passive material parameters
  const double alpha = params_->alpha_;
  const double beta = params_->beta_;
  const double gamma = params_->gamma_;
  const double kappa = params_->kappa_;
  const double omega0 = params_->omega0_;

  // compute matrices
  // right Cauchy Green tensor C
  Core::LinAlg::Matrix<3, 3> C(false);                     // matrix notation
  C.MultiplyTN(*defgrd, *defgrd);                          // C = F^T F
  Core::LinAlg::Matrix<6, 1> Cv(false);                    // Voigt notation
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(C, Cv);  // Cv

  // inverse right Cauchy Green tensor C^-1
  Core::LinAlg::Matrix<3, 3> invC(false);                        // matrix notation
  invC.Invert(C);                                                // invC = C^-1
  Core::LinAlg::Matrix<6, 1> invCv(false);                       // Voigt notation
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(invC, invCv);  // invCv

  // structural tensor M, i.e. dyadic product of fibre directions
  Core::LinAlg::Matrix<3, 3> M = anisotropy_extension_.get_structural_tensor(gp, 0);
  Core::LinAlg::Matrix<6, 1> Mv(false);                    // Voigt notation
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(M, Mv);  // Mv
  // structural tensor L = omega0/3*Identity + omegap*M
  Core::LinAlg::Matrix<3, 3> L(M);
  L.Scale(1.0 - omega0);  // omegap*M
  for (unsigned i = 0; i < 3; ++i) L(i, i) += omega0 / 3.0;

  // product invC*L
  Core::LinAlg::Matrix<3, 3> invCL(false);
  invCL.MultiplyNN(invC, L);

  // product invC*L*invC
  Core::LinAlg::Matrix<3, 3> invCLinvC(false);  // matrix notation
  invCLinvC.MultiplyNN(invCL, invC);
  Core::LinAlg::Matrix<6, 1> invCLinvCv(false);  // Voigt notation
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(invCLinvC, invCLinvCv);

  // stretch in fibre direction lambdaM
  // lambdaM = sqrt(C:M) = sqrt(tr(C^T M)), see Holzapfel2000, p.14
  double lambdaM = Mat::UTILS::Muscle::FiberStretch(C, M);

  // computation of active nominal stress Pa, and derivative derivPa
  double intPa = 0.0;
  double Pa = 0.0;
  double derivPa = 0.0;
  if (params_->Popt_ != 0)
  {  // if active material
    evaluate_active_nominal_stress(params, eleGID, lambdaM, intPa, Pa, derivPa);
  }  // else: intPa, Pa and derivPa remain 0.0

  // computation of activation level omegaa and derivative \frac{\partial omegaa}{\partial C}
  double omegaa = 0.0;
  double derivOmegaa = 0.0;
  double derivDerivOmegaa = 0.0;

  // prefactor eta and its first derivative w.r.t. lambdaM
  double eta = 0;
  double dEta = 0;

  // compute activation level and derivative only if active nominal stress is not zero
  if (Pa >= 1E-12)
  {
    evaluate_activation_level(lambdaM, intPa, Pa, derivPa, omegaa, derivOmegaa, derivDerivOmegaa);

    // passive part of invariant I and its first and second derivatives w.r.t. lambdaM
    double Ip =
        (omega0 / 3) * (std::pow(lambdaM, 2) + 2 / lambdaM) + (1 - omega0) * std::pow(lambdaM, 2);
    double dIp =
        -2 * lambdaM * (omega0 - 1) + (omega0 / 3) * (2 * lambdaM - 2 / std::pow(lambdaM, 2));

    // invariant I and its first derivative w.r.t. lambdaM
    double I = Ip + omegaa * std::pow(lambdaM, 2);
    double dI = dIp + 2 * omegaa * lambdaM + 2 * std::pow(lambdaM, 2) * derivOmegaa;

    // prefactor eta and its first derivative w.r.t. lambdaM
    eta = exp(alpha * (I - 1)) * lambdaM * derivOmegaa;
    dEta = exp(alpha * (I - 1)) * lambdaM *
           (derivOmegaa * (1 / lambdaM + alpha * dI) + derivDerivOmegaa);
  }

  // compute derivative \frac{\partial omegaa}{\partial C} in Voigt notation
  Core::LinAlg::Matrix<6, 1> domegaadCv(Mv);
  domegaadCv.Scale(derivOmegaa * 0.5 / lambdaM);

  // compute helper matrices for further calculation
  Core::LinAlg::Matrix<3, 3> LomegaaM(L);
  LomegaaM.Update(omegaa, M, 1.0);  // LomegaaM = L + omegaa*M
  Core::LinAlg::Matrix<6, 1> LomegaaMv(false);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(LomegaaM, LomegaaMv);

  Core::LinAlg::Matrix<3, 3> LfacomegaaM(L);  // LfacomegaaM = L + fac*M
  LfacomegaaM.Update((1.0 + omegaa * alpha * std::pow(lambdaM, 2)) / (alpha * std::pow(lambdaM, 2)),
      M, 1.0);  // + fac*M
  Core::LinAlg::Matrix<6, 1> LfacomegaaMv(false);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(LfacomegaaM, LfacomegaaMv);

  Core::LinAlg::Matrix<3, 3> transpCLomegaaM(false);
  transpCLomegaaM.MultiplyTN(1.0, C, LomegaaM);  // C^T*(L+omegaa*M)
  Core::LinAlg::Matrix<6, 1> transpCLomegaaMv(false);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(transpCLomegaaM, transpCLomegaaMv);

  // generalized invariants including active material properties
  double detC = C.Determinant();  // detC = det(C)
  // I = C:(L+omegaa*M) = tr(C^T (L+omegaa*M)) since A:B = tr(A^T B) for real matrices
  double I = transpCLomegaaM(0, 0) + transpCLomegaaM(1, 1) + transpCLomegaaM(2, 2);
  // J = cof(C):L = tr(cof(C)^T L) = tr(adj(C) L) = tr(det(C) C^-1 L) = det(C)*tr(C^-1 L)
  double J = detC * (invCL(0, 0) + invCL(1, 1) + invCL(2, 2));
  // exponential prefactors
  double expalpha = std::exp(alpha * (I - 1.0));
  double expbeta = std::exp(beta * (J - 1.0));

  // compute second Piola-Kirchhoff stress
  Core::LinAlg::Matrix<3, 3> stressM(false);
  stressM.Update(expalpha, LomegaaM, 0.0);  // add contributions
  stressM.Update(-expbeta * detC, invCLinvC, 1.0);
  stressM.Update(J * expbeta - std::pow(detC, -kappa), invC, 1.0);
  stressM.Update(0.5 * eta, M, 1.0);
  stressM.Scale(0.5 * gamma);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(
      stressM, Sc_stress);  // convert to Voigt notation and update stress

  // compute cmat
  ccmat.MultiplyNT(alpha * expalpha, LomegaaMv, LomegaaMv, 1.0);  // add contributions
  ccmat.MultiplyNT(alpha * std::pow(lambdaM, 2) * expalpha, LfacomegaaMv, domegaadCv, 1.0);
  ccmat.MultiplyNT(beta * expbeta * std::pow(detC, 2.), invCLinvCv, invCLinvCv, 1.0);
  ccmat.MultiplyNT(-(beta * J + 1.) * expbeta * detC, invCv, invCLinvCv, 1.0);
  ccmat.MultiplyNT(-(beta * J + 1.) * expbeta * detC, invCLinvCv, invCv, 1.0);
  ccmat.MultiplyNT(
      (beta * J + 1.) * J * expbeta + kappa * std::pow(detC, -kappa), invCv, invCv, 1.0);
  // adds scalar * (invC boeppel invC) to cmat, see Holzapfel2000, p. 254
  Mat::add_holzapfel_product(ccmat, invCv, -(J * expbeta - std::pow(detC, -kappa)));
  // adds -expbeta*detC * dinvCLinvCdCv to cmats
  Mat::add_derivative_of_inva_b_inva_product(-expbeta * detC, invCv, invCLinvCv, ccmat);
  // additional term for corrected derivative of strain energy function
  ccmat.MultiplyNT(dEta / (8 * lambdaM), Mv, Mv, 1.0);
  ccmat.Scale(gamma);

  // update stress and material tangent with the computed stress and cmat values
  stress->Update(1.0, Sc_stress, 1.0);
  cmat->Update(1.0, ccmat, 1.0);
}

void Mat::MuscleCombo::evaluate_active_nominal_stress(Teuchos::ParameterList& params,
    const int eleGID, const double lambdaM, double& intPa, double& Pa, double& derivPa)
{
  // save current simulation time
  double t_tot = params.get<double>("total time", -1);
  if (abs(t_tot + 1.0) < 1e-14) FOUR_C_THROW("No total time given for muscle Combo material!");
  // save (time) step size
  double timestep = params.get<double>("delta time", -1);
  if (abs(timestep + 1.0) < 1e-14)
    FOUR_C_THROW("No time step size given for muscle Combo material!");

  // get active parameters from params_
  const double lambdaMin = params_->lambdaMin_;
  const double lambdaOpt = params_->lambdaOpt_;

  const double Popt = params_->Popt_;

  // get element center coordinates in reference configuration
  const auto& element_center_reference_coordinates =
      params.get<Core::LinAlg::Matrix<3, 1>>("elecenter_coords_ref");

  // compute force-time-space dependency Poptft
  const double Poptft =
      std::visit(ActivationEvalVisitor{Popt, t_tot, element_center_reference_coordinates, eleGID},
          activation_evaluator_);

  // compute the force-stretch dependency fxi, its integral in the boundaries lambdaMin to lambdaM,
  // and its derivative w.r.t. lambdaM
  double intFxi = Mat::UTILS::Muscle::EvaluateIntegralForceStretchDependencyEhret(
      lambdaM, lambdaMin, lambdaOpt);
  double fxi =
      Mat::UTILS::Muscle::EvaluateForceStretchDependencyEhret(lambdaM, lambdaMin, lambdaOpt);
  double dFxidLamdaM = Mat::UTILS::Muscle::EvaluateDerivativeForceStretchDependencyEhret(
      lambdaM, lambdaMin, lambdaOpt);

  // compute active nominal stress Pa, its integral in the boundaries lambdaMin to lambdaM,
  // and its derivative w.r.t. lambdaM
  intPa = Poptft * intFxi;
  Pa = Poptft * fxi;
  derivPa = Poptft * dFxidLamdaM;
}

void Mat::MuscleCombo::evaluate_activation_level(const double lambdaM, const double intPa,
    const double Pa, const double derivPa, double& omegaa, double& derivOmegaa,
    double& derivDerivOmegaa)
{
  // get passive material parameters
  const double alpha = params_->alpha_;
  const double gamma = params_->gamma_;
  const double omega0 = params_->omega0_;

  // passive part of invariant I and its first and second derivatives w.r.t. lambdaM
  double Ip = (omega0 / 3.0) * (std::pow(lambdaM, 2) + 2.0 / lambdaM) +
              (1.0 - omega0) * std::pow(lambdaM, 2);
  double dIp = (omega0 / 3.0) * (2.0 * lambdaM - 2.0 / std::pow(lambdaM, 2)) +
               2.0 * (1.0 - omega0) * lambdaM;
  double ddIp = (omega0 / 3.0) * (2.0 + 4.0 / std::pow(lambdaM, 3)) + 2.0 * (1.0 - omega0);

  // helper tau and its first and second derivatives w.r.t. lambdaM
  double tau = alpha * (1 - Ip);
  double dTau = -alpha * dIp;
  double ddTau = -alpha * ddIp;

  // helper phi and its first and second derivatives w.r.t. lambdaM
  double phi = 1 + ((4 * alpha) / gamma) * intPa * std::exp(tau);
  double dPhi = (4 * alpha / gamma) * (intPa * std::exp(tau) * dTau + std::exp(tau) * Pa);
  double ddPhi = (4 * alpha / gamma) * std::exp(tau) *
                 (2 * Pa * dTau + intPa * std::pow(dTau, 2) + intPa * ddTau + derivPa);

  // computation of activation level omegaa
  omegaa = std::log(phi) / (alpha * std::pow(lambdaM, 2));

  // computation of partial derivative of omegaa w.r.t. lambdaM
  derivOmegaa = -2 * std::log(phi) / (alpha * std::pow(lambdaM, 3)) +
                dPhi / (phi * alpha * std::pow(lambdaM, 2));
  derivDerivOmegaa = 6 * std::log(phi) / (alpha * std::pow(lambdaM, 4)) -
                     4 * dPhi / (phi * alpha * std::pow(lambdaM, 3)) -
                     std::pow(dPhi, 2) / (phi * phi * alpha * std::pow(lambdaM, 2)) +
                     ddPhi / (phi * alpha * std::pow(lambdaM, 2));
}
FOUR_C_NAMESPACE_CLOSE
