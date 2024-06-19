/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of the Weickenmeier active skeletal muscle material (generalized active strain
approach)

\level 3

*/
/*----------------------------------------------------------------------*/

#include "4C_mat_muscle_weickenmeier.hpp"

#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_mat_muscle_utils.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"
#include "4C_matelast_aniso_structuraltensor_strategy.hpp"

FOUR_C_NAMESPACE_OPEN


Mat::PAR::MuscleWeickenmeier::MuscleWeickenmeier(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      alpha_(matdata.parameters.get<double>("ALPHA")),
      beta_(matdata.parameters.get<double>("BETA")),
      gamma_(matdata.parameters.get<double>("GAMMA")),
      kappa_(matdata.parameters.get<double>("KAPPA")),
      omega0_(matdata.parameters.get<double>("OMEGA0")),
      Na_(matdata.parameters.get<double>("ACTMUNUM")),
      muTypesNum_(matdata.parameters.get<int>("MUTYPESNUM")),
      I_((matdata.parameters.get<std::vector<double>>("INTERSTIM"))),
      rho_((matdata.parameters.get<std::vector<double>>("FRACACTMU"))),
      F_((matdata.parameters.get<std::vector<double>>("FTWITCH"))),
      T_((matdata.parameters.get<std::vector<double>>("TTWITCH"))),
      lambdaMin_(matdata.parameters.get<double>("LAMBDAMIN")),
      lambdaOpt_(matdata.parameters.get<double>("LAMBDAOPT")),
      dotLambdaMMin_(matdata.parameters.get<double>("DOTLAMBDAMIN")),
      ke_(matdata.parameters.get<double>("KE")),
      kc_(matdata.parameters.get<double>("KC")),
      de_(matdata.parameters.get<double>("DE")),
      dc_(matdata.parameters.get<double>("DC")),
      actTimesNum_(matdata.parameters.get<int>("ACTTIMESNUM")),
      actTimes_((matdata.parameters.get<std::vector<double>>("ACTTIMES"))),
      actIntervalsNum_(matdata.parameters.get<int>("ACTINTERVALSNUM")),
      actValues_((matdata.parameters.get<std::vector<double>>("ACTVALUES"))),
      density_(matdata.parameters.get<double>("DENS"))
{
  // error handling for parameter ranges
  // passive material parameters
  if (alpha_ <= 0.0) FOUR_C_THROW("Material parameter ALPHA must be greater zero");
  if (beta_ <= 0.0) FOUR_C_THROW("Material parameter BETA must be greater zero");
  if (gamma_ <= 0.0) FOUR_C_THROW("Material parameter GAMMA must be greater zero");
  if (omega0_ < 0.0 || omega0_ > 1.0) FOUR_C_THROW("Material parameter OMEGA0 must be in [0;1]");

  // active material parameters
  // stimulation frequency dependent parameters
  if (Na_ < 0.0)
  {
    FOUR_C_THROW("Material parameter ACTMUNUM must be postive or zero");
  }

  double sumrho = 0.0;
  for (int iMU = 0; iMU < muTypesNum_; ++iMU)
  {
    if (I_[iMU] < 0.0) FOUR_C_THROW("Material parameter INTERSTIM must be postive or zero");
    if (rho_[iMU] < 0.0) FOUR_C_THROW("Material parameter FRACACTMU must be postive or zero");

    sumrho += rho_[iMU];
    if (F_[iMU] < 0.0) FOUR_C_THROW("Material parameter FTWITCH must be postive or zero");
    if (T_[iMU] < 0.0) FOUR_C_THROW("Material parameter TTWITCH must be postive or zero");
  }

  if (muTypesNum_ > 1 && sumrho != 1.0) FOUR_C_THROW("Sum of fractions of MU types must equal one");

  // stretch dependent parameters
  if (lambdaMin_ <= 0.0) FOUR_C_THROW("Material parameter LAMBDAMIN must be postive");
  if (lambdaOpt_ <= 0.0) FOUR_C_THROW("Material parameter LAMBDAOPT must be postive");

  // velocity dependent parameters
  if (ke_ < 0.0) FOUR_C_THROW("Material parameter KE should be postive or zero");
  if (kc_ < 0.0) FOUR_C_THROW("Material parameter KC should be postive or zero");
  if (de_ < 0.0) FOUR_C_THROW("Material parameter DE should be postive or zero");
  if (dc_ < 0.0) FOUR_C_THROW("Material parameter DC should be postive or zero");

  // prescribed activation in time intervals
  if (actTimesNum_ != int(actTimes_.size()))
    FOUR_C_THROW("Number of activation times ACTTIMES must equal ACTTIMESNUM");
  if (actIntervalsNum_ != int(actValues_.size()))
    FOUR_C_THROW("Number of activation values ACTVALUES must equal ACTINTERVALSNUM");
  if (actTimesNum_ != actIntervalsNum_ + 1)
    FOUR_C_THROW("ACTTIMESNUM must be one smaller than ACTINTERVALSNUM");

  // density
  if (density_ < 0.0) FOUR_C_THROW("DENS should be positive");
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::MuscleWeickenmeier::create_material()
{
  return Teuchos::rcp(new Mat::MuscleWeickenmeier(this));
}

Mat::MuscleWeickenmeierType Mat::MuscleWeickenmeierType::instance_;

Core::Communication::ParObject* Mat::MuscleWeickenmeierType::Create(const std::vector<char>& data)
{
  auto* muscle_weickenmeier = new Mat::MuscleWeickenmeier();
  muscle_weickenmeier->unpack(data);
  return muscle_weickenmeier;
}

Mat::MuscleWeickenmeier::MuscleWeickenmeier()
    : params_(nullptr),
      lambda_m_old_(1.0),
      anisotropy_(),
      anisotropy_extension_(true, 0.0, 0,
          Teuchos::rcp<Mat::Elastic::StructuralTensorStrategyBase>(
              new Mat::Elastic::StructuralTensorStrategyStandard(nullptr)),
          {0})
{
}

Mat::MuscleWeickenmeier::MuscleWeickenmeier(Mat::PAR::MuscleWeickenmeier* params)
    : params_(params),
      lambda_m_old_(1.0),
      anisotropy_(),
      anisotropy_extension_(true, 0.0, 0,
          Teuchos::rcp<Mat::Elastic::StructuralTensorStrategyBase>(
              new Mat::Elastic::StructuralTensorStrategyStandard(nullptr)),
          {0})
{
  // initialize lambdaMOld_
  lambda_m_old_ = 1.0;

  // register anisotropy extension to global anisotropy
  anisotropy_.register_anisotropy_extension(anisotropy_extension_);

  // initialize fiber directions and structural tensor
  anisotropy_extension_.register_needed_tensors(
      Mat::FiberAnisotropyExtension<1>::FIBER_VECTORS |
      Mat::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR);
}

void Mat::MuscleWeickenmeier::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  add_to_pack(data, lambda_m_old_);

  anisotropy_extension_.pack_anisotropy(data);
}

void Mat::MuscleWeickenmeier::unpack(const std::vector<char>& data)
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
        params_ = static_cast<Mat::PAR::MuscleWeickenmeier*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  extract_from_pack(position, data, lambda_m_old_);

  anisotropy_extension_.unpack_anisotropy(data, position);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

void Mat::MuscleWeickenmeier::setup(int numgp, Input::LineDefinition* linedef)
{
  // Read anisotropy
  anisotropy_.set_number_of_gauss_points(numgp);
  anisotropy_.read_anisotropy_from_element(linedef);
}

void Mat::MuscleWeickenmeier::Update(Core::LinAlg::Matrix<3, 3> const& defgrd, int const gp,
    Teuchos::ParameterList& params, int const eleGID)
{
  // compute the current fibre stretch using the deformation gradient and the structural tensor
  // right Cauchy Green tensor C= F^T F
  Core::LinAlg::Matrix<3, 3> C(false);
  C.MultiplyTN(defgrd, defgrd);

  // structural tensor M, i.e. dyadic product of fibre directions
  const Core::LinAlg::Matrix<3, 3>& M = anisotropy_extension_.get_structural_tensor(gp, 0);

  // save the current fibre stretch in lambdaMOld_
  lambda_m_old_ = Mat::UTILS::Muscle::FiberStretch(C, M);
}

void Mat::MuscleWeickenmeier::evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,
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
  double Pa = 0.0;
  double derivPa = 0.0;
  if (params_->muTypesNum_ != 0)
  {  // if active material
    evaluate_active_nominal_stress(params, lambdaM, Pa, derivPa);
  }  // else: Pa and derivPa remain 0.0

  // computation of activation level omegaa and derivative w.r.t. fiber stretch
  double omegaa = 0.0;
  double derivOmegaa = 0.0;
  // compute activation level and derivative only if active nominal stress is not zero
  // if active nominal stress is zero, material is purely passive, thus activation level and
  // derivative are zero
  if (Pa != 0.0)
  {
    evaluate_activation_level(lambdaM, Pa, derivPa, omegaa, derivOmegaa);
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
  LfacomegaaM.Update(
      (1.0 + omegaa * alpha * std::pow(lambdaM, 2.)) / (alpha * std::pow(lambdaM, 2.)), M,
      1.0);  // + fac*M
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
  stressM.Scale(0.5 * gamma);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(
      stressM, Sc_stress);  // convert to Voigt notation and update stress

  // compute cmat
  ccmat.MultiplyNT(alpha * expalpha, LomegaaMv, LomegaaMv, 1.0);  // add contributions
  ccmat.MultiplyNT(alpha * std::pow(lambdaM, 2.) * expalpha, LfacomegaaMv, domegaadCv, 1.0);
  ccmat.MultiplyNT(beta * expbeta * std::pow(detC, 2.), invCLinvCv, invCLinvCv, 1.0);
  ccmat.MultiplyNT(-(beta * J + 1.) * expbeta * detC, invCv, invCLinvCv, 1.0);
  ccmat.MultiplyNT(-(beta * J + 1.) * expbeta * detC, invCLinvCv, invCv, 1.0);
  ccmat.MultiplyNT(
      (beta * J + 1.) * J * expbeta + kappa * std::pow(detC, -kappa), invCv, invCv, 1.0);
  // adds scalar * (invC boeppel invC) to cmat, see Holzapfel2000, p. 254
  Mat::add_holzapfel_product(ccmat, invCv, -(J * expbeta - std::pow(detC, -kappa)));
  // adds -expbeta*detC * dinvCLinvCdCv to cmats
  Mat::add_derivative_of_inva_b_inva_product(-expbeta * detC, invCv, invCLinvCv, ccmat);
  ccmat.Scale(gamma);

  // update stress and material tangent with the computed stress and cmat values
  stress->Update(1.0, Sc_stress, 1.0);
  cmat->Update(1.0, ccmat, 1.0);
}

void Mat::MuscleWeickenmeier::evaluate_active_nominal_stress(
    Teuchos::ParameterList& params, const double lambdaM, double& Pa, double& derivPa)
{
  // save current simulation time
  double t_tot = params.get<double>("total time", -1);
  if (std::abs(t_tot + 1.0) < 1e-14)
    FOUR_C_THROW("No total time given for muscle Weickenmeier material!");
  // save (time) step size
  double timestep = params.get<double>("delta time", -1);
  if (std::abs(timestep + 1.0) < 1e-14)
    FOUR_C_THROW("No time step size given for muscle Weickenmeier material!");

  // approximate first time derivative of lambdaM through BW Euler
  // dotLambdaM = (lambdaM_n - lambdaM_{n-1})/dt
  double dotLambdaM = (lambdaM - lambda_m_old_) / timestep;

  // approximate second time derivative of lambdaM through BW Euler
  // dDotLambdaMdLambdaM = 1/dt approximated through BW Euler
  double dDotLambdaMdLambdaM = 1 / timestep;

  // get active microstructural parameters from params_
  const double Na = params_->Na_;
  const int muTypesNum = params_->muTypesNum_;
  const auto& I = params_->I_;
  const auto& rho = params_->rho_;
  const auto& F = params_->F_;
  const auto& T = params_->T_;

  const double lambdaMin = params_->lambdaMin_;
  const double lambdaOpt = params_->lambdaOpt_;

  const double dotLambdaMMin = params_->dotLambdaMMin_;
  const double ke = params_->ke_;
  const double kc = params_->kc_;
  const double de = params_->de_;
  const double dc = params_->dc_;

  const int actIntervalsNum = params_->actIntervalsNum_;
  const auto& actTimes = params_->actTimes_;
  const auto& actValues = params_->actValues_;

  // compute force-time/stimulation frequency dependency Poptft
  double Poptft = Mat::UTILS::Muscle::EvaluateTimeDependentActiveStressEhret(
      Na, muTypesNum, rho, I, F, T, actIntervalsNum, actTimes, actValues, t_tot);

  // compute force-stretch dependency fxi
  double fxi =
      Mat::UTILS::Muscle::EvaluateForceStretchDependencyEhret(lambdaM, lambdaMin, lambdaOpt);

  // compute force-velocity dependency fv
  double fv = Mat::UTILS::Muscle::EvaluateForceVelocityDependencyBoel(
      dotLambdaM, dotLambdaMMin, de, dc, ke, kc);

  // compute active nominal stress Pa
  Pa = Poptft * fxi * fv;

  // compute derivative of force-stretch dependency fxi and of force-velocity dependency fv
  // w.r.t. lambdaM
  double dFxidLamdaM = 0.0;
  double dFvdLambdaM = 0.0;
  if (Pa != 0)
  {
    dFxidLamdaM = Mat::UTILS::Muscle::EvaluateDerivativeForceStretchDependencyEhret(
        lambdaM, lambdaMin, lambdaOpt);
    dFvdLambdaM = Mat::UTILS::Muscle::EvaluateDerivativeForceVelocityDependencyBoel(
        dotLambdaM, dDotLambdaMdLambdaM, dotLambdaMMin, de, dc, ke, kc);
  }

  // compute derivative of active nominal stress Pa w.r.t. lambdaM
  derivPa = Poptft * (fv * dFxidLamdaM + fxi * dFvdLambdaM);
}

void Mat::MuscleWeickenmeier::evaluate_activation_level(const double lambdaM, const double Pa,
    const double derivPa, double& omegaa, double& derivOmegaa)
{
  // get passive material parameters
  const double alpha = params_->alpha_;
  const double gamma = params_->gamma_;
  const double omega0 = params_->omega0_;

  // passive part of invariant I and its first and second derivatives w.r.t. lambdaM
  double Ip = (omega0 / 3.0) * (std::pow(lambdaM, 2.) + 2.0 / lambdaM) +
              (1.0 - omega0) * std::pow(lambdaM, 2.);
  double derivIp = (omega0 / 3.0) * (2.0 * lambdaM - 2.0 / std::pow(lambdaM, 2.)) +
                   2.0 * (1.0 - omega0) * lambdaM;
  double derivderivIp = (omega0 / 3.0) * (2.0 + 4.0 / std::pow(lambdaM, 3.)) + 2.0 * (1.0 - omega0);

  // argument for Lambert W function
  const double xi =
      Pa * ((2.0 * alpha * lambdaM) / gamma) *
          std::exp(0.5 * alpha * (2.0 - 2.0 * Ip + lambdaM * derivIp)) +
      0.5 * alpha * lambdaM * derivIp * std::exp(0.5 * alpha * lambdaM * derivIp);  // argument xi

  // solution W0 of principal branch of Lambert W function approximated with Halley's method
  double W0 = 1.0;           // starting guess for solution
  const double tol = 1e-15;  // tolerance for numeric approximation 10^-15
  const int maxiter = 100;   // maximal number of iterations
  Mat::UTILS::Muscle::EvaluateLambert(xi, W0, tol, maxiter);

  // derivatives of xi and W0 w.r.t. lambdaM used for activation level computation
  double derivXi =
      (2.0 * alpha / gamma * std::exp(0.5 * alpha * (2.0 - 2.0 * Ip + lambdaM * derivIp))) *
          (Pa + lambdaM * derivPa +
              0.5 * alpha * Pa * lambdaM * (lambdaM * derivderivIp - derivIp)) +
      0.5 * alpha * (1.0 + 0.5 * alpha) * std::exp(0.5 * alpha * lambdaM * derivIp) *
          (derivIp + lambdaM * derivderivIp);
  double derivLambert = derivXi / ((1.0 + W0) * std::exp(W0));

  // computation of activation level omegaa
  omegaa = W0 / (alpha * std::pow(lambdaM, 2.)) - derivIp / (2.0 * lambdaM);

  // computation of partial derivative of omegaa w.r.t. lambdaM
  derivOmegaa = derivLambert / (alpha * lambdaM * lambdaM) -
                2.0 * W0 / (alpha * lambdaM * lambdaM * lambdaM) - derivderivIp / (2.0 * lambdaM) +
                derivIp / (2.0 * lambdaM * lambdaM);
}
FOUR_C_NAMESPACE_CLOSE
