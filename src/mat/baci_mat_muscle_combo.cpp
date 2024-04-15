/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of the Combo active skeletal muscle material with variable time-dependent
activations

\level 3

*/
/*----------------------------------------------------------------------*/

#include "baci_mat_muscle_combo.hpp"

#include "baci_global_data.hpp"
#include "baci_io_linedefinition.hpp"
#include "baci_linalg_fixedsizematrix_voigt_notation.hpp"
#include "baci_mat_muscle_utils.hpp"
#include "baci_mat_par_bundle.hpp"
#include "baci_mat_service.hpp"
#include "baci_matelast_aniso_structuraltensor_strategy.hpp"
#include "baci_utils_exceptions.hpp"

#include <variant>

FOUR_C_NAMESPACE_OPEN

namespace
{
  using ActivationMapType = std::unordered_map<int, std::vector<std::pair<double, double>>>;

  MAT::PAR::MuscleCombo::ActivationParameterVariant GetActivationParams(
      const Teuchos::RCP<MAT::PAR::Material>& matdata,
      const INPAR::MAT::ActivationType& activation_type)
  {
    if (activation_type == INPAR::MAT::ActivationType::function_of_space_time)
    {
      auto actFunctId = *matdata->Get<int>("ACTFUNCT");
      if (actFunctId <= 0) dserror("Function id must be positive");
      return actFunctId;
    }
    else if (activation_type == INPAR::MAT::ActivationType::map_from_csv)
    {
      return *matdata->Get<const ActivationMapType>("ACTCSV");
    }
    else
      return std::monostate{};
  }

  struct ActivationParamsVisitor
  {
    MAT::MuscleCombo::ActivationEvaluatorVariant operator()(const int function_id) const
    {
      return &GLOBAL::Problem::Instance()->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(
          function_id - 1);
    }

    MAT::MuscleCombo::ActivationEvaluatorVariant operator()(
        const ActivationMapType& activation_csv_data_map) const
    {
      return &activation_csv_data_map;
    }

    MAT::MuscleCombo::ActivationEvaluatorVariant operator()(const std::monostate& /*unused*/) const
    {
      dserror(
          "Error in ActivationParamsVisitor. You're calling it with the default-constructed "
          "state.");
    }
  };

  struct ActivationEvalVisitor
  {
    double operator()(const ActivationMapType*& csvMap) const
    {
      // use one-based element ids in the csv files (corresponding to the ones in the dat-file)
      return MAT::UTILS::MUSCLE::EvaluateTimeSpaceDependentActiveStressByFunct(
          Popt_, *csvMap, t_tot_, eleGID_ + 1);
    }

    double operator()(const CORE::UTILS::FunctionOfSpaceTime*& spaceTimeFunction) const
    {
      return MAT::UTILS::MUSCLE::EvaluateTimeSpaceDependentActiveStressByFunct(
          Popt_, *spaceTimeFunction, t_tot_, element_center_reference_coordinates_);
    }

    double operator()(const std::monostate& /*unused*/) const
    {
      dserror(
          "Error in ActivationEvalVisitor. You're calling it with the default-constructed state.");
    }

    const double& Popt_;
    const double& t_tot_;
    const CORE::LINALG::Matrix<3, 1>& element_center_reference_coordinates_;
    const int& eleGID_;
  };
}  // namespace


MAT::PAR::MuscleCombo::MuscleCombo(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      alpha_(*matdata->Get<double>("ALPHA")),
      beta_(*matdata->Get<double>("BETA")),
      gamma_(*matdata->Get<double>("GAMMA")),
      kappa_(*matdata->Get<double>("KAPPA")),
      omega0_(*matdata->Get<double>("OMEGA0")),
      Popt_(*matdata->Get<double>("POPT")),
      lambdaMin_(*matdata->Get<double>("LAMBDAMIN")),
      lambdaOpt_(*matdata->Get<double>("LAMBDAOPT")),
      activationType_(static_cast<INPAR::MAT::ActivationType>(*matdata->Get<int>("ACTEVALTYPE"))),
      activationParams_(GetActivationParams(matdata, activationType_)),
      density_(*matdata->Get<double>("DENS"))
{
  printf("ACTEVALTYOPE IS %d", *matdata->Get<int>("ACTEVALTYPE"));
  // error handling for parameter ranges
  // passive material parameters
  if (alpha_ <= 0.0) dserror("Material parameter ALPHA must be greater zero");
  if (beta_ <= 0.0) dserror("Material parameter BETA must be greater zero");
  if (gamma_ <= 0.0) dserror("Material parameter GAMMA must be greater zero");
  if (omega0_ < 0.0 || omega0_ > 1.0) dserror("Material parameter OMEGA0 must be in [0;1]");

  // active material parameters
  if (Popt_ < 0.0)
  {
    dserror("Material parameter POPT must be postive or zero");
  }

  // stretch dependent parameters
  if (lambdaMin_ <= 0.0) dserror("Material parameter LAMBDAMIN must be postive");
  if (lambdaOpt_ <= 0.0) dserror("Material parameter LAMBDAOPT must be postive");

  // density
  if (density_ < 0.0) dserror("DENS should be positive");
}

Teuchos::RCP<MAT::Material> MAT::PAR::MuscleCombo::CreateMaterial()
{
  return Teuchos::rcp(new MAT::MuscleCombo(this));
}

MAT::MuscleComboType MAT::MuscleComboType::instance_;

CORE::COMM::ParObject* MAT::MuscleComboType::Create(const std::vector<char>& data)
{
  auto* muscle_combo = new MAT::MuscleCombo();
  muscle_combo->Unpack(data);
  return muscle_combo;
}

MAT::MuscleCombo::MuscleCombo()
    : params_(nullptr),
      anisotropy_(),
      anisotropyExtension_(true, 0.0, 0,
          Teuchos::rcp<MAT::ELASTIC::StructuralTensorStrategyBase>(
              new MAT::ELASTIC::StructuralTensorStrategyStandard(nullptr)),
          {0}),
      activationEvaluator_(std::monostate{})
{
}

MAT::MuscleCombo::MuscleCombo(MAT::PAR::MuscleCombo* params)
    : params_(params),
      anisotropy_(),
      anisotropyExtension_(true, 0.0, 0,
          Teuchos::rcp<MAT::ELASTIC::StructuralTensorStrategyBase>(
              new MAT::ELASTIC::StructuralTensorStrategyStandard(nullptr)),
          {0}),
      activationEvaluator_(std::monostate{})
{
  // register anisotropy extension to global anisotropy
  anisotropy_.RegisterAnisotropyExtension(anisotropyExtension_);

  // initialize fiber directions and structural tensor
  anisotropyExtension_.RegisterNeededTensors(
      MAT::FiberAnisotropyExtension<1>::FIBER_VECTORS |
      MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR |
      MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR_STRESS);

  // cannot set activation_function here, because function manager did not yet read functions
}

void MAT::MuscleCombo::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  anisotropyExtension_.PackAnisotropy(data);
}

void MAT::MuscleCombo::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // make sure we have a pristine material
  params_ = nullptr;

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);

  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
      {
        params_ = static_cast<MAT::PAR::MuscleCombo*>(mat);
        activationEvaluator_ = std::visit(ActivationParamsVisitor(), params_->activationParams_);
      }
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  anisotropyExtension_.UnpackAnisotropy(data, position);

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

void MAT::MuscleCombo::Setup(int numgp, INPUT::LineDefinition* linedef)
{
  // Read anisotropy
  anisotropy_.SetNumberOfGaussPoints(numgp);
  anisotropy_.ReadAnisotropyFromElement(linedef);

  activationEvaluator_ = std::visit(ActivationParamsVisitor(), params_->activationParams_);
}

void MAT::MuscleCombo::Update(CORE::LINALG::Matrix<3, 3> const& defgrd, int const gp,
    Teuchos::ParameterList& params, int const eleGID)
{
}

void MAT::MuscleCombo::Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
    const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<6, 1>* stress, CORE::LINALG::Matrix<6, 6>* cmat, const int gp,
    const int eleGID)
{
  CORE::LINALG::Matrix<6, 1> Sc_stress(true);
  CORE::LINALG::Matrix<6, 6> ccmat(true);

  // get passive material parameters
  const double alpha = params_->alpha_;
  const double beta = params_->beta_;
  const double gamma = params_->gamma_;
  const double kappa = params_->kappa_;
  const double omega0 = params_->omega0_;

  // compute matrices
  // right Cauchy Green tensor C
  CORE::LINALG::Matrix<3, 3> C(false);                   // matrix notation
  C.MultiplyTN(*defgrd, *defgrd);                        // C = F^T F
  CORE::LINALG::Matrix<6, 1> Cv(false);                  // Voigt notation
  CORE::LINALG::VOIGT::Stresses::MatrixToVector(C, Cv);  // Cv

  // inverse right Cauchy Green tensor C^-1
  CORE::LINALG::Matrix<3, 3> invC(false);                      // matrix notation
  invC.Invert(C);                                              // invC = C^-1
  CORE::LINALG::Matrix<6, 1> invCv(false);                     // Voigt notation
  CORE::LINALG::VOIGT::Stresses::MatrixToVector(invC, invCv);  // invCv

  // structural tensor M, i.e. dyadic product of fibre directions
  CORE::LINALG::Matrix<3, 3> M = anisotropyExtension_.GetStructuralTensor(gp, 0);
  CORE::LINALG::Matrix<6, 1> Mv(false);                  // Voigt notation
  CORE::LINALG::VOIGT::Stresses::MatrixToVector(M, Mv);  // Mv
  // structural tensor L = omega0/3*Identity + omegap*M
  CORE::LINALG::Matrix<3, 3> L(M);
  L.Scale(1.0 - omega0);  // omegap*M
  for (unsigned i = 0; i < 3; ++i) L(i, i) += omega0 / 3.0;

  // product invC*L
  CORE::LINALG::Matrix<3, 3> invCL(false);
  invCL.MultiplyNN(invC, L);

  // product invC*L*invC
  CORE::LINALG::Matrix<3, 3> invCLinvC(false);  // matrix notation
  invCLinvC.MultiplyNN(invCL, invC);
  CORE::LINALG::Matrix<6, 1> invCLinvCv(false);  // Voigt notation
  CORE::LINALG::VOIGT::Stresses::MatrixToVector(invCLinvC, invCLinvCv);

  // stretch in fibre direction lambdaM
  // lambdaM = sqrt(C:M) = sqrt(tr(C^T M)), see Holzapfel2000, p.14
  double lambdaM = MAT::UTILS::MUSCLE::FiberStretch(C, M);

  // computation of active nominal stress Pa, and derivative derivPa
  double intPa = 0.0;
  double Pa = 0.0;
  double derivPa = 0.0;
  if (params_->Popt_ != 0)
  {  // if active material
    EvaluateActiveNominalStress(params, eleGID, lambdaM, intPa, Pa, derivPa);
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
    EvaluateActivationLevel(lambdaM, intPa, Pa, derivPa, omegaa, derivOmegaa, derivDerivOmegaa);

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
  CORE::LINALG::Matrix<6, 1> domegaadCv(Mv);
  domegaadCv.Scale(derivOmegaa * 0.5 / lambdaM);

  // compute helper matrices for further calculation
  CORE::LINALG::Matrix<3, 3> LomegaaM(L);
  LomegaaM.Update(omegaa, M, 1.0);  // LomegaaM = L + omegaa*M
  CORE::LINALG::Matrix<6, 1> LomegaaMv(false);
  CORE::LINALG::VOIGT::Stresses::MatrixToVector(LomegaaM, LomegaaMv);

  CORE::LINALG::Matrix<3, 3> LfacomegaaM(L);  // LfacomegaaM = L + fac*M
  LfacomegaaM.Update((1.0 + omegaa * alpha * std::pow(lambdaM, 2)) / (alpha * std::pow(lambdaM, 2)),
      M, 1.0);  // + fac*M
  CORE::LINALG::Matrix<6, 1> LfacomegaaMv(false);
  CORE::LINALG::VOIGT::Stresses::MatrixToVector(LfacomegaaM, LfacomegaaMv);

  CORE::LINALG::Matrix<3, 3> transpCLomegaaM(false);
  transpCLomegaaM.MultiplyTN(1.0, C, LomegaaM);  // C^T*(L+omegaa*M)
  CORE::LINALG::Matrix<6, 1> transpCLomegaaMv(false);
  CORE::LINALG::VOIGT::Stresses::MatrixToVector(transpCLomegaaM, transpCLomegaaMv);

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
  CORE::LINALG::Matrix<3, 3> stressM(false);
  stressM.Update(expalpha, LomegaaM, 0.0);  // add contributions
  stressM.Update(-expbeta * detC, invCLinvC, 1.0);
  stressM.Update(J * expbeta - std::pow(detC, -kappa), invC, 1.0);
  stressM.Update(0.5 * eta, M, 1.0);
  stressM.Scale(0.5 * gamma);
  CORE::LINALG::VOIGT::Stresses::MatrixToVector(
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
  MAT::AddtoCmatHolzapfelProduct(ccmat, invCv, -(J * expbeta - std::pow(detC, -kappa)));
  // adds -expbeta*detC * dinvCLinvCdCv to cmats
  MAT::AddDerivInvABInvBProduct(-expbeta * detC, invCv, invCLinvCv, ccmat);
  // additional term for corrected derivative of strain energy function
  ccmat.MultiplyNT(dEta / (8 * lambdaM), Mv, Mv, 1.0);
  ccmat.Scale(gamma);

  // update stress and material tangent with the computed stress and cmat values
  stress->Update(1.0, Sc_stress, 1.0);
  cmat->Update(1.0, ccmat, 1.0);
}

void MAT::MuscleCombo::EvaluateActiveNominalStress(Teuchos::ParameterList& params, const int eleGID,
    const double lambdaM, double& intPa, double& Pa, double& derivPa)
{
  // save current simulation time
  double t_tot = params.get<double>("total time", -1);
  if (abs(t_tot + 1.0) < 1e-14) dserror("No total time given for muscle Combo material!");
  // save (time) step size
  double timestep = params.get<double>("delta time", -1);
  if (abs(timestep + 1.0) < 1e-14) dserror("No time step size given for muscle Combo material!");

  // get active parameters from params_
  const double lambdaMin = params_->lambdaMin_;
  const double lambdaOpt = params_->lambdaOpt_;

  const double Popt = params_->Popt_;

  // get element center coordinates in reference configuration
  const auto& element_center_reference_coordinates =
      params.get<CORE::LINALG::Matrix<3, 1>>("elecenter_coords_ref");

  // compute force-time-space dependency Poptft
  const double Poptft =
      std::visit(ActivationEvalVisitor{Popt, t_tot, element_center_reference_coordinates, eleGID},
          activationEvaluator_);

  // compute the force-stretch dependency fxi, its integral in the boundaries lambdaMin to lambdaM,
  // and its derivative w.r.t. lambdaM
  double intFxi = MAT::UTILS::MUSCLE::EvaluateIntegralForceStretchDependencyEhret(
      lambdaM, lambdaMin, lambdaOpt);
  double fxi =
      MAT::UTILS::MUSCLE::EvaluateForceStretchDependencyEhret(lambdaM, lambdaMin, lambdaOpt);
  double dFxidLamdaM = MAT::UTILS::MUSCLE::EvaluateDerivativeForceStretchDependencyEhret(
      lambdaM, lambdaMin, lambdaOpt);

  // compute active nominal stress Pa, its integral in the boundaries lambdaMin to lambdaM,
  // and its derivative w.r.t. lambdaM
  intPa = Poptft * intFxi;
  Pa = Poptft * fxi;
  derivPa = Poptft * dFxidLamdaM;
}

void MAT::MuscleCombo::EvaluateActivationLevel(const double lambdaM, const double intPa,
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
