/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of the Weickenmeier active skeletal muscle material

\level 3

*/
/*----------------------------------------------------------------------*/

#include "baci_mat_muscle_weickenmeier.H"

#include "baci_lib_globalproblem.H"
#include "baci_lib_linedefinition.H"
#include "baci_lib_voigt_notation.H"
#include "baci_mat_muscle_utils.H"
#include "baci_mat_par_bundle.H"
#include "baci_mat_service.H"
#include "baci_matelast_aniso_structuraltensor_strategy.H"


MAT::PAR::Muscle_Weickenmeier::Muscle_Weickenmeier(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      alpha_(matdata->GetDouble("ALPHA")),
      beta_(matdata->GetDouble("BETA")),
      gamma_(matdata->GetDouble("GAMMA")),
      kappa_(matdata->GetDouble("KAPPA")),
      omega0_(matdata->GetDouble("OMEGA0")),
      Na_(matdata->GetDouble("ACTMUNUM")),
      muTypesNum_(matdata->GetInt("MUTYPESNUM")),
      I_(*(matdata->Get<std::vector<double>>("INTERSTIM"))),
      rho_(*(matdata->Get<std::vector<double>>("FRACACTMU"))),
      F_(*(matdata->Get<std::vector<double>>("FTWITCH"))),
      T_(*(matdata->Get<std::vector<double>>("TTWITCH"))),
      lambdaMin_(matdata->GetDouble("LAMBDAMIN")),
      lambdaOpt_(matdata->GetDouble("LAMBDAOPT")),
      dotLambdaMMin_(matdata->GetDouble("DOTLAMBDAMIN")),
      ke_(matdata->GetDouble("KE")),
      kc_(matdata->GetDouble("KC")),
      de_(matdata->GetDouble("DE")),
      dc_(matdata->GetDouble("DC")),
      actTimesNum_(matdata->GetInt("ACTTIMESNUM")),
      actTimes_(*(matdata->Get<std::vector<double>>("ACTTIMES"))),
      actIntervalsNum_(matdata->GetInt("ACTINTERVALSNUM")),
      actValues_(*(matdata->Get<std::vector<double>>("ACTVALUES"))),
      density_(matdata->GetDouble("DENS"))
{
  // error handling for parameter ranges
  // passive material parameters
  if (alpha_ <= 0.0) dserror("Material parameter ALPHA must be greater zero");
  if (beta_ <= 0.0) dserror("Material parameter BETA must be greater zero");
  if (gamma_ <= 0.0) dserror("Material parameter GAMMA must be greater zero");
  if (omega0_ < 0.0 || omega0_ > 1.0) dserror("Material parameter OMEGA0 must be in [0;1]");

  // active material parameters
  // stimulation frequency dependent parameters
  if (Na_ < 0.0)
  {
    dserror("Material parameter ACTMUNUM must be postive or zero");
  }

  double sumrho = 0.0;
  for (int iMU = 0; iMU < muTypesNum_; ++iMU)
  {
    if (I_[iMU] < 0.0) dserror("Material parameter INTERSTIM must be postive or zero");
    if (rho_[iMU] < 0.0) dserror("Material parameter FRACACTMU must be postive or zero");

    sumrho += rho_[iMU];
    if (F_[iMU] < 0.0) dserror("Material parameter FTWITCH must be postive or zero");
    if (T_[iMU] < 0.0) dserror("Material parameter TTWITCH must be postive or zero");
  }

  if (muTypesNum_ > 1 && sumrho != 1.0) dserror("Sum of fractions of MU types must equal one");

  // stretch dependent parameters
  if (lambdaMin_ <= 0.0) dserror("Material parameter LAMBDAMIN must be postive");
  if (lambdaOpt_ <= 0.0) dserror("Material parameter LAMBDAOPT must be postive");

  // velocity dependent parameters
  if (ke_ < 0.0) dserror("Material parameter KE should be postive or zero");
  if (kc_ < 0.0) dserror("Material parameter KC should be postive or zero");
  if (de_ < 0.0) dserror("Material parameter DE should be postive or zero");
  if (dc_ < 0.0) dserror("Material parameter DC should be postive or zero");

  // prescribed activation in time intervals
  if (actTimesNum_ != int(actTimes_.size()))
    dserror("Number of activation times ACTTIMES must equal ACTTIMESNUM");
  if (actIntervalsNum_ != int(actValues_.size()))
    dserror("Number of activation values ACTVALUES must equal ACTINTERVALSNUM");
  if (actTimesNum_ != actIntervalsNum_ + 1)
    dserror("ACTTIMESNUM must be one smaller than ACTINTERVALSNUM");

  // density
  if (density_ < 0.0) dserror("DENS should be positive");
}

Teuchos::RCP<MAT::Material> MAT::PAR::Muscle_Weickenmeier::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Muscle_Weickenmeier(this));
}

MAT::Muscle_WeickenmeierType MAT::Muscle_WeickenmeierType::instance_;

DRT::ParObject* MAT::Muscle_WeickenmeierType::Create(const std::vector<char>& data)
{
  auto* muscle_weickenmeier = new MAT::Muscle_Weickenmeier();
  muscle_weickenmeier->Unpack(data);
  return muscle_weickenmeier;
}

MAT::Muscle_Weickenmeier::Muscle_Weickenmeier()
    : params_(nullptr),
      lambdaMOld_(1.0),
      anisotropy_(),
      anisotropyExtension_(true, 0.0, 0,
          Teuchos::rcp<MAT::ELASTIC::StructuralTensorStrategyBase>(
              new MAT::ELASTIC::StructuralTensorStrategyStandard(nullptr)),
          {0})
{
}

MAT::Muscle_Weickenmeier::Muscle_Weickenmeier(MAT::PAR::Muscle_Weickenmeier* params)
    : params_(params),
      lambdaMOld_(1.0),
      anisotropy_(),
      anisotropyExtension_(true, 0.0, 0,
          Teuchos::rcp<MAT::ELASTIC::StructuralTensorStrategyBase>(
              new MAT::ELASTIC::StructuralTensorStrategyStandard(nullptr)),
          {0})
{
  // initialize lambdaMOld_
  lambdaMOld_ = 1.0;

  // register anisotropy extension to global anisotropy
  anisotropy_.RegisterAnisotropyExtension(anisotropyExtension_);

  // initialize fiber directions and structural tensor
  anisotropyExtension_.RegisterNeededTensors(MAT::FiberAnisotropyExtension<1>::FIBER_VECTORS |
                                             MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR);
}

void MAT::Muscle_Weickenmeier::Pack(DRT::PackBuffer& data) const
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

  AddtoPack(data, lambdaMOld_);

  anisotropyExtension_.PackAnisotropy(data);
}

void MAT::Muscle_Weickenmeier::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId())
    dserror(
        "Wrong instance type data. The extracted type id is %d, while the UniqueParObjectId is %d",
        type, UniqueParObjectId());


  // make sure we have a pristine material
  params_ = nullptr;

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);

  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Muscle_Weickenmeier*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  ExtractfromPack(position, data, lambdaMOld_);

  anisotropyExtension_.UnpackAnisotropy(data, position);

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

void MAT::Muscle_Weickenmeier::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // Read anisotropy
  anisotropy_.SetNumberOfGaussPoints(numgp);
  anisotropy_.ReadAnisotropyFromElement(linedef);
}

void MAT::Muscle_Weickenmeier::Update(CORE::LINALG::Matrix<3, 3> const& defgrd, int const gp,
    Teuchos::ParameterList& params, int const eleGID)
{
  // compute the current fibre stretch using the deformation gradient and the structural tensor
  // right Cauchy Green tensor C= F^T F
  CORE::LINALG::Matrix<3, 3> C(false);
  C.MultiplyTN(defgrd, defgrd);

  // structural tensor M, i.e. dyadic product of fibre directions
  const CORE::LINALG::Matrix<3, 3>& M = anisotropyExtension_.GetStructuralTensor(gp, 0);

  // save the current fibre stretch in lambdaMOld_
  lambdaMOld_ = MAT::UTILS::MUSCLE::FiberStretch(C, M);
}

void MAT::Muscle_Weickenmeier::Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
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
  CORE::LINALG::Matrix<3, 3> C(false);              // matrix notation
  C.MultiplyTN(*defgrd, *defgrd);                   // C = F^T F
  CORE::LINALG::Matrix<6, 1> Cv(false);             // Voigt notation
  ::UTILS::VOIGT::Stresses::MatrixToVector(C, Cv);  // Cv

  // inverse right Cauchy Green tensor C^-1
  CORE::LINALG::Matrix<3, 3> invC(false);                 // matrix notation
  invC.Invert(C);                                         // invC = C^-1
  CORE::LINALG::Matrix<6, 1> invCv(false);                // Voigt notation
  ::UTILS::VOIGT::Stresses::MatrixToVector(invC, invCv);  // invCv

  // structural tensor M, i.e. dyadic product of fibre directions
  CORE::LINALG::Matrix<3, 3> M = anisotropyExtension_.GetStructuralTensor(gp, 0);
  CORE::LINALG::Matrix<6, 1> Mv(false);             // Voigt notation
  ::UTILS::VOIGT::Stresses::MatrixToVector(M, Mv);  // Mv

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
  ::UTILS::VOIGT::Stresses::MatrixToVector(invCLinvC, invCLinvCv);

  // stretch in fibre direction lambdaM
  // lambdaM = sqrt(C:M) = sqrt(tr(C^T M)), see Holzapfel2000, p.14
  double lambdaM = MAT::UTILS::MUSCLE::FiberStretch(C, M);

  // computation of active nominal stress Pa, and derivative derivPa
  double Pa = 0.0;
  double derivPa = 0.0;
  if (params_->muTypesNum_ != 0)
  {  // if active material
    EvaluateActiveNominalStress(params, lambdaM, Pa, derivPa);
  }  // else: Pa and derivPa remain 0.0

  // computation of activation level omegaa and derivative w.r.t. fiber stretch
  double omegaa = 0.0;
  double derivOmegaa = 0.0;
  // compute activation level and derivative only if active nominal stress is not zero
  // if active nominal stress is zero, material is purely passive, thus activation level and
  // derivative are zero
  if (Pa != 0.0)
  {
    EvaluateActivationLevel(lambdaM, Pa, derivPa, omegaa, derivOmegaa);
  }
  // compute derivative \frac{\partial omegaa}{\partial C} in Voigt notation
  CORE::LINALG::Matrix<6, 1> domegaadCv(Mv);
  domegaadCv.Scale(derivOmegaa);

  // compute helper matrices for further calculation
  CORE::LINALG::Matrix<3, 3> LomegaaM(L);
  LomegaaM.Update(omegaa, M, 1.0);  // LomegaaM = L + omegaa*M
  CORE::LINALG::Matrix<6, 1> LomegaaMv(false);
  ::UTILS::VOIGT::Stresses::MatrixToVector(LomegaaM, LomegaaMv);

  CORE::LINALG::Matrix<3, 3> LfacomegaaM(L);  // LfacomegaaM = L + fac*M
  LfacomegaaM.Update(
      (1.0 + omegaa * alpha * std::pow(lambdaM, 2.)) / (alpha * std::pow(lambdaM, 2.)), M,
      1.0);  // + fac*M
  CORE::LINALG::Matrix<6, 1> LfacomegaaMv(false);
  ::UTILS::VOIGT::Stresses::MatrixToVector(LfacomegaaM, LfacomegaaMv);

  CORE::LINALG::Matrix<3, 3> transpCLomegaaM(false);
  transpCLomegaaM.MultiplyTN(1.0, C, LomegaaM);  // C^T*(L+omegaa*M)
  CORE::LINALG::Matrix<6, 1> transpCLomegaaMv(false);
  ::UTILS::VOIGT::Stresses::MatrixToVector(transpCLomegaaM, transpCLomegaaMv);

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
  stressM.Scale(0.5 * gamma);
  ::UTILS::VOIGT::Stresses::MatrixToVector(
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
  MAT::AddtoCmatHolzapfelProduct(ccmat, invCv, -(J * expbeta - std::pow(detC, -kappa)));
  // adds -expbeta*detC * dinvCLinvCdCv to cmats
  MAT::AddDerivInvABInvBProduct(-expbeta * detC, invCv, invCLinvCv, ccmat);
  ccmat.Scale(gamma);

  // update stress and material tangent with the computed stress and cmat values
  stress->Update(1.0, Sc_stress, 1.0);
  cmat->Update(1.0, ccmat, 1.0);
}

void MAT::Muscle_Weickenmeier::EvaluateActiveNominalStress(
    Teuchos::ParameterList& params, const double lambdaM, double& Pa, double& derivPa)
{
  // save current simulation time
  double t_tot = params.get<double>("total time", -1);
  if (abs(t_tot + 1.0) < 1e-14) dserror("No total time given for muscle Weickenmeier material!");
  // save (time) step size
  double timestep = params.get<double>("delta time", -1);
  if (abs(timestep + 1.0) < 1e-14)
    dserror("No time step size given for muscle Weickenmeier material!");

  // approximate first time derivative of lambdaM through BW Euler
  // dotLambdaM = (lambdaM_n - lambdaM_{n-1})/dt
  double dotLambdaM = (lambdaM - lambdaMOld_) / timestep;

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
  double Poptft = MAT::UTILS::MUSCLE::EvaluateTimeDependentActiveStressEhret(
      Na, muTypesNum, rho, I, F, T, actIntervalsNum, actTimes, actValues, t_tot);

  // compute force-stretch dependency fxi
  double fxi =
      MAT::UTILS::MUSCLE::EvaluateForceStretchDependencyEhret(lambdaM, lambdaMin, lambdaOpt);

  // compute force-velocity dependency fv
  double fv = MAT::UTILS::MUSCLE::EvaluateForceVelocityDependencyBoel(
      dotLambdaM, dotLambdaMMin, de, dc, ke, kc);

  // compute active nominal stress Pa
  Pa = Poptft * fxi * fv;

  // compute derivative of force-stretch dependency fxi and of force-velocity dependency fv
  // w.r.t. lambdaM
  double dFxidLamdaM = 0.0;
  double dFvdLambdaM = 0.0;
  if (Pa != 0)
  {
    dFxidLamdaM = MAT::UTILS::MUSCLE::EvaluateDerivativeForceStretchDependencyEhret(
        lambdaM, lambdaMin, lambdaOpt);
    dFvdLambdaM = MAT::UTILS::MUSCLE::EvaluateDerivativeForceVelocityDependencyBoel(
        dotLambdaM, dDotLambdaMdLambdaM, dotLambdaMMin, de, dc, ke, kc);
  }

  // compute derivative of active nominal stress Pa w.r.t. lambdaM
  derivPa = Poptft * (fv * dFxidLamdaM + fxi * dFvdLambdaM);
}

void MAT::Muscle_Weickenmeier::EvaluateActivationLevel(const double lambdaM, const double Pa,
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
  MAT::UTILS::MUSCLE::EvaluateLambert(xi, W0, tol, maxiter);

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
  derivOmegaa = derivLambert / (2.0 * alpha * std::pow(lambdaM, 3.)) -
                W0 / (alpha * std::pow(lambdaM, 4)) - derivderivIp / (4.0 * std::pow(lambdaM, 2.)) +
                derivIp / (4.0 * std::pow(lambdaM, 3.));
}