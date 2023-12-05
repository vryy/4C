/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of the Giantesio active strain skeletal muscle material (active strain
approach)

\level 3

*/
/*----------------------------------------------------------------------*/

#include "baci_mat_muscle_giantesio.H"

#include "baci_io_linedefinition.H"
#include "baci_lib_globalproblem.H"
#include "baci_linalg_fixedsizematrix_voigt_notation.H"
#include "baci_linalg_four_tensor.H"
#include "baci_mat_muscle_utils.H"
#include "baci_mat_par_bundle.H"
#include "baci_mat_service.H"
#include "baci_matelast_aniso_structuraltensor_strategy.H"
#include "baci_utils_local_newton.H"


namespace
{
  /*!
   * @brief Returns the fourth order tensor
   *        result_ijkl = scalar * (tensor_A_ijkm * matrix_B_mn * matrix_C_nl
   *               + matrix_A_im * tensor_B_mjkn * matrix_C_nl
   *               + matrix_A_im * matrix_B_mn * tensor_C_njkl)
   */
  CORE::LINALG::FourTensor<3> SumMultiplyTmmMtmMmt(const CORE::LINALG::Matrix<3, 3>& matrix_A,
      const CORE::LINALG::Matrix<3, 3>& matrix_B, const CORE::LINALG::Matrix<3, 3>& matrix_C,
      const CORE::LINALG::FourTensor<3>& tensor_A, const CORE::LINALG::FourTensor<3>& tensor_B,
      const CORE::LINALG::FourTensor<3>& tensor_C, const double& scalar)
  {
    CORE::LINALG::Matrix<3, 3> mBmC(false);
    mBmC.MultiplyNN(matrix_B, matrix_C);

    CORE::LINALG::Matrix<3, 3> mAmB(false);
    mAmB.MultiplyNN(matrix_A, matrix_B);

    // part1 = tensor_A * matrix_B * matrix_C
    CORE::LINALG::FourTensor<3> part1(true);
    MAT::MultiplyMatrixFourTensorBySecondIndex<3>(part1, mBmC, tensor_A, false);

    // part2 = matrix_A * tensor_B * matrix_C
    CORE::LINALG::FourTensor<3> mAtB(true);  // matrix_A * tensor_B
    MAT::MultiplyMatrixFourTensor<3>(mAtB, matrix_A, tensor_B, false);
    CORE::LINALG::FourTensor<3> part2(true);  // (matrix_A * tensor_B) * matrix_C
    MAT::MultiplyMatrixFourTensorBySecondIndex<3>(part2, matrix_C, mAtB, false);

    // part3 = matrix_A * matrix_B * tensor_C
    CORE::LINALG::FourTensor<3> part3(true);
    MAT::MultiplyMatrixFourTensor<3>(part3, mAmB, tensor_C, false);

    // result = part1 + part2 + part3
    CORE::LINALG::FourTensor<3> tAmBmCmAtBmCmAmBtC(true);
    for (unsigned i = 0; i < 3; ++i)
      for (unsigned j = 0; j < 3; ++j)
        for (unsigned k = 0; k < 3; ++k)
          for (unsigned l = 0; l < 3; ++l)
            tAmBmCmAtBmCmAmBtC(i, j, k, l) +=
                scalar * (part1(i, j, k, l) + part2(i, j, k, l) + part3(i, j, k, l));

    return tAmBmCmAtBmCmAmBtC;
  }

  /*!
   * @brief Adds the triple multiplication to the result matrix
   *        NNN_ij += scalar * left_ij middle_jk right_kl
   */
  void MultiplyNNN(CORE::LINALG::Matrix<3, 3>& NNN, const double& scalar,
      const CORE::LINALG::Matrix<3, 3>& left, const CORE::LINALG::Matrix<3, 3>& middle,
      const CORE::LINALG::Matrix<3, 3>& right)
  {
    CORE::LINALG::Matrix<3, 3> NN(true);
    NN.MultiplyNN(left, middle);
    NNN.MultiplyNN(scalar, NN, right, 1.0);
  }

  /*!
   * @brief Compute the derivative of the Cauchy Green tensor C w.r.t. C
   */
  CORE::LINALG::FourTensor<3> ComputeDCDC()
  {
    CORE::LINALG::Matrix<6, 6> dCdCv(true);
    for (int i = 0; i < 3; i++) dCdCv(i, i) = 1.0;
    for (int i = 3; i < 6; i++) dCdCv(i, i) = 0.5;

    CORE::LINALG::FourTensor<3> dCdC(true);
    MAT::SetupFourTensor(dCdC, dCdCv);

    return dCdC;
  }

  /*!
   * @brief Compute the structural tensor L
   */
  CORE::LINALG::Matrix<3, 3> ComputeStructuralTensorL(
      const CORE::LINALG::Matrix<3, 3>& M, const double& omega0)
  {
    CORE::LINALG::Matrix<3, 3> L(M);
    L.Scale(1.0 - omega0);  // omegap*M
    for (unsigned i = 0; i < 3; ++i) L(i, i) += omega0 / 3.0;

    return L;
  }

  struct StressAndDeriv
  {
    CORE::LINALG::Matrix<3, 3> Se;
    CORE::LINALG::FourTensor<3> dSedC;
  };

  /*!
   * @brief Compute the elastic second Piola Kirchhoff stress tensor Se and its derivative dSedC
   * w.r.t. the right Cauchy Green tensor C
   */
  StressAndDeriv ComputeSeAndDSedC(const double& alpha, const double& beta, const double& gamma,
      const double& omega0, const CORE::LINALG::Matrix<3, 3>& M,
      const CORE::LINALG::Matrix<3, 3>& invFa, const CORE::LINALG::FourTensor<3>& dinvFadC,
      const CORE::LINALG::Matrix<3, 3>& C, const CORE::LINALG::FourTensor<3>& dCdC)
  {
    // structural tensor L = omega0/3*Identity + omegap*M
    CORE::LINALG::Matrix<3, 3> L = ComputeStructuralTensorL(M, omega0);
    CORE::LINALG::Matrix<6, 1> Lv(false);  // Voigt notation
    CORE::LINALG::VOIGT::Stresses::MatrixToVector(L, Lv);

    // elastic right Cauchy Green tensor Ce = Fe^T Fe
    // = Fa^-T C Fa^-1 = Fa^-1 C Fa^-1 (with Fa^-1 sym.)
    CORE::LINALG::Matrix<3, 3> Ce(true);
    MultiplyNNN(Ce, 1.0, invFa, C, invFa);

    // derivative of Ce w.r.t C
    // dCedC_ijkl = dinvFadC_iakl (C invFa)_aj + invFa_ia dCdC_abkl invFa_bj + (invFa C)_ab
    // dinvFadC_bjkl
    CORE::LINALG::FourTensor<3> dCedC(true);
    dCedC = SumMultiplyTmmMtmMmt(invFa, C, invFa, dinvFadC, dCdC, dinvFadC, 1.0);
    CORE::LINALG::Matrix<6, 6> dCedCv(true);
    MAT::Setup6x6VoigtMatrix(dCedCv, dCedC);

    // inverse of the elastic right Cauchy Green tensor Ce
    CORE::LINALG::Matrix<3, 3> invCe(true);
    invCe.Invert(Ce);
    CORE::LINALG::Matrix<6, 1> invCev(true);
    CORE::LINALG::VOIGT::Stresses::MatrixToVector(invCe, invCev);

    // product invCe*L
    CORE::LINALG::Matrix<3, 3> invCeL(true);
    invCeL.MultiplyNN(invCe, L);

    // product invCe*L*invCe
    CORE::LINALG::Matrix<3, 3> invCeLinvCe(true);
    invCeLinvCe.MultiplyNN(invCeL, invCe);
    CORE::LINALG::Matrix<6, 1> invCeLinvCev(true);
    CORE::LINALG::VOIGT::Stresses::MatrixToVector(invCeLinvCe, invCeLinvCev);

    // product Ce^T*L
    CORE::LINALG::Matrix<3, 3> transpCeL(true);
    transpCeL.MultiplyTN(Ce, L);

    // generalized elastic invariants including only passive properties
    double detCe = Ce.Determinant();  // determinant of Ce
    // Ie = C:L = tr(C^T L) since A:B = tr(A^T B) for real matrices
    double Ie = transpCeL(0, 0) + transpCeL(1, 1) + transpCeL(2, 2);
    // Je = cof(Ce):L = tr(cof(Ce)^T L) = tr(adj(Ce) L) = tr(det(Ce) Ce^-1 L) = det(Ce)*tr(Ce^-1 L)
    double Je = detCe * (invCeL(0, 0) + invCeL(1, 1) + invCeL(2, 2));
    // exponential prefactors
    double expalpha = std::exp(alpha * (Ie - 1.0));
    double expbeta = std::exp(beta * (Je - 1.0));

    // elastic second Piola Kirchhoff stress tensor Se
    // Se = (0.5 * gamma * (expalpha * L + expbeta * (Je * invCe - invCeLinvCe)))
    CORE::LINALG::Matrix<3, 3> Se(true);
    Se.Update(expalpha, L);
    Se.Update(-expbeta * detCe, invCeLinvCe, 1.0);
    Se.Update(Je * expbeta, invCe, 1.0);
    Se.Scale(0.5 * gamma);

    // derivative of Se w.r.t Ce
    // taken from Weickenmeier et al.
    CORE::LINALG::Matrix<6, 6> dSedCev(true);
    dSedCev.MultiplyNT(alpha * expalpha, Lv, Lv, 1.0);  // add contributions
    dSedCev.MultiplyNT(beta * expbeta * std::pow(detCe, 2), invCeLinvCev, invCeLinvCev, 1.0);
    dSedCev.MultiplyNT(-(beta * Je + 1.) * expbeta * detCe, invCev, invCeLinvCev, 1.0);
    dSedCev.MultiplyNT(-(beta * Je + 1.) * expbeta * detCe, invCeLinvCev, invCev, 1.0);
    dSedCev.MultiplyNT((beta * Je + 1.) * Je * expbeta, invCev, invCev, 1.0);
    // adds scalar * (invC boeppel invC) to cmat, see Holzapfel2000, p. 254
    MAT::AddtoCmatHolzapfelProduct(dSedCev, invCev, -Je * expbeta);
    // adds -expbeta * detCe * dinvCLinvCdCv to cmat
    MAT::AddDerivInvABInvBProduct(-expbeta * detCe, invCev, invCeLinvCev, dSedCev);
    dSedCev.Scale(gamma / 2);

    CORE::LINALG::FourTensor<3> dSedCe(true);
    MAT::SetupFourTensor(dSedCe, dSedCev);

    // derivative of Se w.r.t C
    // dSedC_ijkl = dSedCe_ijab dCedC_abkl
    CORE::LINALG::FourTensor<3> dSedC(true);
    MAT::MultiplyFourTensorFourTensor<3>(dSedC, dSedCe, dCedC, true);

    StressAndDeriv Se_dSedC = {.Se = Se, .dSedC = dSedC};

    return Se_dSedC;
  }
}  // namespace

MAT::PAR::Muscle_Giantesio::Muscle_Giantesio(Teuchos::RCP<MAT::PAR::Material> matdata)
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

Teuchos::RCP<MAT::Material> MAT::PAR::Muscle_Giantesio::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Muscle_Giantesio(this));
}

MAT::Muscle_GiantesioType MAT::Muscle_GiantesioType::instance_;

CORE::COMM::ParObject* MAT::Muscle_GiantesioType::Create(const std::vector<char>& data)
{
  auto* muscle_giantesio = new MAT::Muscle_Giantesio();
  muscle_giantesio->Unpack(data);
  return muscle_giantesio;
}

MAT::Muscle_Giantesio::Muscle_Giantesio()
    : params_(nullptr),
      lambdaMOld_(1.0),
      omegaaOld_(-1.0),
      anisotropy_(),
      anisotropyExtension_(true, 0.0, 0,
          Teuchos::rcp<MAT::ELASTIC::StructuralTensorStrategyBase>(
              new MAT::ELASTIC::StructuralTensorStrategyStandard(nullptr)),
          {0})
{
}

MAT::Muscle_Giantesio::Muscle_Giantesio(MAT::PAR::Muscle_Giantesio* params)
    : params_(params),
      lambdaMOld_(1.0),
      omegaaOld_(-1.0),
      anisotropy_(),
      anisotropyExtension_(true, 0.0, 0,
          Teuchos::rcp<MAT::ELASTIC::StructuralTensorStrategyBase>(
              new MAT::ELASTIC::StructuralTensorStrategyStandard(nullptr)),
          {0})
{
  // initialize lambdaMOld_ and omegaOld_
  lambdaMOld_ = 1.0;
  omegaaOld_ = -1.0;

  // register anisotropy extension to global anisotropy
  anisotropy_.RegisterAnisotropyExtension(anisotropyExtension_);

  // initialize fiber directions and structural tensor
  anisotropyExtension_.RegisterNeededTensors(MAT::FiberAnisotropyExtension<1>::FIBER_VECTORS |
                                             MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR);
}

void MAT::Muscle_Giantesio::Pack(CORE::COMM::PackBuffer& data) const
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

  AddtoPack(data, lambdaMOld_);
  AddtoPack(data, omegaaOld_);

  anisotropyExtension_.PackAnisotropy(data);
}

void MAT::Muscle_Giantesio::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

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
        params_ = static_cast<MAT::PAR::Muscle_Giantesio*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  ExtractfromPack(position, data, lambdaMOld_);
  ExtractfromPack(position, data, omegaaOld_);

  anisotropyExtension_.UnpackAnisotropy(data, position);

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

void MAT::Muscle_Giantesio::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // Read anisotropy
  anisotropy_.SetNumberOfGaussPoints(numgp);
  anisotropy_.ReadAnisotropyFromElement(linedef);
}

void MAT::Muscle_Giantesio::Update(CORE::LINALG::Matrix<3, 3> const& defgrd, int const gp,
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

void MAT::Muscle_Giantesio::Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,
    const CORE::LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<6, 1>* stress, CORE::LINALG::Matrix<6, 6>* cmat, const int gp,
    const int eleGID)
{
  stress->Clear();
  cmat->Clear();

  // get passive material parameters
  const double alpha = params_->alpha_;
  const double beta = params_->beta_;
  const double gamma = params_->gamma_;
  const double kappa = params_->kappa_;
  const double omega0 = params_->omega0_;

  // save current simulation time
  double currentTime = params.get<double>("total time", -1);
  if (std::abs(currentTime + 1.0) < 1e-14)
    dserror("No total time given for muscle Giantesio material!");

  // save (time) step size
  double timeStepSize = params.get<double>("delta time", -1);
  if (std::abs(timeStepSize + 1.0) < 1e-14)
    dserror("No time step size given for muscle Giantesio material!");

  // compute matrices
  // right Cauchy Green tensor C
  CORE::LINALG::Matrix<3, 3> C(false);  // matrix notation
  C.MultiplyTN(*defgrd, *defgrd);       // C = F^T F

  // determinant of C
  double detC = C.Determinant();

  // derivative of C w.r.t C
  CORE::LINALG::FourTensor<3> dCdC = ComputeDCDC();

  // inverse right Cauchy Green tensor C^-1
  CORE::LINALG::Matrix<3, 3> invC(true);                       // matrix notation
  invC.Invert(C);                                              // invC = C^-1
  CORE::LINALG::Matrix<6, 1> invCv(true);                      // Voigt notation
  CORE::LINALG::VOIGT::Stresses::MatrixToVector(invC, invCv);  // invCv

  // structural tensor M, i.e. dyadic product of fibre directions
  CORE::LINALG::Matrix<3, 3> M = anisotropyExtension_.GetStructuralTensor(gp, 0);

  // stretch in fibre direction lambdaM
  double lambdaM = MAT::UTILS::MUSCLE::FiberStretch(C, M);

  // derivative of lambdaM w.r.t. C
  CORE::LINALG::Matrix<3, 3> dlambdaMdC = MAT::UTILS::MUSCLE::DFiberStretch_DC(lambdaM, C, M);
  CORE::LINALG::Matrix<6, 1> dlambdaMdCv(true);
  CORE::LINALG::VOIGT::Stresses::MatrixToVector(dlambdaMdC, dlambdaMdCv);

  // contraction velocity dotLambdaM
  double dotLambdaM =
      MAT::UTILS::MUSCLE::ContractionVelocityBWEuler(lambdaM, lambdaMOld_, timeStepSize);

  // compute activation level omegaa and derivative w.r.t. fiber stretch if material is active
  CORE::UTILS::ValuesFunctAndFunctDerivs omegaaAndDerivs = {
      .val_funct = 0.0, .val_deriv_funct = 0.0, .val_deriv_deriv_funct = 0.0};

  if (IsActive(currentTime))
  {
    omegaaAndDerivs = EvaluateActivationLevelAndDerivatives(lambdaM, dotLambdaM, currentTime);
    omegaaOld_ = omegaaAndDerivs.val_funct;  // update omegaaOld
  }
  else
  {
    omegaaOld_ = -1;  // reset omegaaOld to -1
  }

  // --------------------------------------------------------------
  // first derivative of omegaa w.r.t. C
  // domegaadC = domegaadlambdaM dlambdaMdC
  CORE::LINALG::Matrix<3, 3> domegaadC(dlambdaMdC);
  domegaadC.Scale(omegaaAndDerivs.val_deriv_funct);
  CORE::LINALG::Matrix<6, 1> domegaadCv(false);
  CORE::LINALG::VOIGT::Stresses::MatrixToVector(domegaadC, domegaadCv);

  // second derivative of omegaa w.r.t. C
  // ddomegaaddC = (ddomegaaddlambdaM - 1/lambdaM domegaadlambdaM) dlambdaM_ij dlambdaM_kl
  CORE::LINALG::Matrix<6, 6> ddomegaaddCv(false);
  ddomegaaddCv.MultiplyNT(
      omegaaAndDerivs.val_deriv_deriv_funct - omegaaAndDerivs.val_deriv_funct / lambdaM,
      dlambdaMdCv, dlambdaMdCv);

  // --------------------------------------------------------------
  // active deformation gradient Fa
  CORE::LINALG::Matrix<3, 3> Fa = ActDefGrad(omegaaAndDerivs.val_funct, M);

  // first derivative of Fa w.r.t. omegaa
  CORE::LINALG::Matrix<3, 3> dFadomegaa = DActDefGrad_DActLevel(omegaaAndDerivs.val_funct, M);

  // second derivative of Fa w.r.t. omegaa
  CORE::LINALG::Matrix<3, 3> ddFaddomegaa = DDActDefGrad_DDActLevel(omegaaAndDerivs.val_funct, M);

  // determinant of Fa
  double detFa = Fa.Determinant();

  // --------------------------------------------------------------
  // inverse of active deformation gradient Fa^{-1}
  CORE::LINALG::Matrix<3, 3> invFa = InvActDefGrad(Fa);

  // first derivative of Fa^{-1} w.r.t. omegaa
  CORE::LINALG::Matrix<3, 3> dinvFadomegaa = DInvActDefGrad_DActLevel(Fa, dFadomegaa);

  // first derivative of Fa^{-1} w.r.t. C
  // dinvFadC_ijkl = dinvFadomegaa_ij domegaadC_kl
  CORE::LINALG::FourTensor<3> dinvFadC(true);
  MAT::AddDyadicProductMatrixMatrix(dinvFadC, dinvFadomegaa, domegaadC);

  // --------------------------------------------------------------
  // elastic second Piola Kirchhoff stress tensor Se and its derivative w.r.t C
  StressAndDeriv Se_dSedC =
      ComputeSeAndDSedC(alpha, beta, gamma, omega0, M, invFa, dinvFadC, C, dCdC);

  // --------------------------------------------------------------
  // compute second Piola Kirchhoff stress components S1, S2, Svol
  // note: S3 = invF * ddetFadF * Psie = 0 since ddetFawa = 0 and ddetFadF = ddetFadwa dwadF

  // S1 = detFa * invFa * Se * invFa^T
  // in index notation: S1_sl = detFa * invFa_sr * Se_rq * invFa_rl (with invFa sym.)
  CORE::LINALG::Matrix<3, 3> S1(true);
  MultiplyNNN(S1, detFa, invFa, Se_dSedC.Se, invFa);
  CORE::LINALG::Matrix<6, 1> S1v(true);
  CORE::LINALG::VOIGT::Stresses::MatrixToVector(S1, S1v);

  // S2 = detFa * invF * dPsiedFe * (F * dinvFadF)  = 2 detFa (C Fa^-1 Se) dFa^-1dC
  // or: S2 = -2 S1 : (C Fa^-1 dFadomegaa) domegaadC
  // in index notation: S2_sl = -2 S1_ij (C Fa^-1 dFadomegaa)_ij domegaadC_sl
  CORE::LINALG::Matrix<3, 3> CinvFadFadomegaa(true);
  MultiplyNNN(CinvFadFadomegaa, 1.0, C, invFa, dFadomegaa);

  CORE::LINALG::Matrix<6, 1> S2v(domegaadCv);
  double scalar_contraction = -2.0 * MAT::ContractMatrixMatrix(CinvFadFadomegaa, S1);
  S2v.Scale(scalar_contraction);

  // Svol = - gamma/2 * detC^-kappa
  CORE::LINALG::Matrix<6, 1> Svolv(invCv);
  Svolv.Scale(-0.5 * gamma * std ::pow(detC, -kappa));

  // --------------------------------------------------------------
  // compute material tangent components cmat1, cmat2, cmatvol

  // compute cmat1 = 2 dS1dC
  CORE::LINALG::Matrix<6, 6> cmat1v(true);

  // derivative of S1 = detFa * invFa * Se * invFa^T  w.r.t. C
  // = detFa * (dinvFadC * (Se * dinvFa^T) + invFa * dSedC * invFa^T + (invFa * Se) * dinvFadC)
  CORE::LINALG::FourTensor<3> dS1dC(true);
  dS1dC =
      SumMultiplyTmmMtmMmt(invFa, Se_dSedC.Se, invFa, dinvFadC, Se_dSedC.dSedC, dinvFadC, detFa);
  MAT::Setup6x6VoigtMatrix(cmat1v, dS1dC);
  cmat1v.Scale(2.0);

  // compute cmat2 = 2 dS2dC
  CORE::LINALG::Matrix<6, 6> cmat2v(true);

  // derivative of S2 w.r.t. C
  // in index notation: S2_sl = -2 S1_ij (C Fa^-1 dFadomegaa)_ij domegaadC_sl
  // dS2dC_slpq = -4 [ S1_ij d(C Fa^-1 dFadomegaa)dC_ijpq domegaadC_sl +
  //                 [ S1_ij (C Fa^-1 dFadomegaa)_ij ddomegaaddC_slpq ]
  // dS2dC_slpq =  2 [ dscalardC_pq domegaadC_sl + [ scalar_contraction ddomegaaddC_slpq ]

  // derivative of the scalar (-2 S1_ij (C Fa^-1 dFadomegaa)_ij) w.r.t. C
  // dscalardC_pq = -2 (dS1dC_ijpq (C invFa dFadomegaa)_ij + S1_ij d(C invFa dFadomegaa)dC_ijpq )
  CORE::LINALG::Matrix<3, 3> H(true);
  H.MultiplyNN(dinvFadomegaa, dFadomegaa);
  H.MultiplyNN(1.0, invFa, ddFaddomegaa, 1.0);
  CORE::LINALG::Matrix<3, 3> CH(true);
  CH.MultiplyNN(C, H);
  double helper = MAT::ContractMatrixMatrix(S1, CH);

  CORE::LINALG::Matrix<3, 3> dscalardC(domegaadC);
  dscalardC.Scale(helper);
  MultiplyNNN(dscalardC, 1.0, S1, dFadomegaa, invFa);
  MAT::AddContractionMatrixFourTensor(dscalardC, CinvFadFadomegaa, dS1dC);
  CORE::LINALG::Matrix<6, 1> dscalardCv(true);
  CORE::LINALG::VOIGT::Stresses::MatrixToVector(dscalardC, dscalardCv);  // not yet scaled with -2.0

  // compute cmat2 = 2 dS2dC = 2 [domegaadC_sl dscalardC_pq + scalar * ddomegaaddC_slpq]
  cmat2v.MultiplyNT(-4.0, domegaadCv, dscalardCv);             // 2 * dscalardC_pq domegaadC_sl
  cmat2v.Update(2.0 * scalar_contraction, ddomegaaddCv, 1.0);  // 2 * scalar * ddomegaaddC_slpq

  // compute cmatvol = 2 dS2dC
  CORE::LINALG::Matrix<6, 6> cmatvolv(true);
  cmatvolv.MultiplyNT(gamma * kappa * std::pow(detC, -kappa), invCv, invCv, 1.0);
  MAT::AddtoCmatHolzapfelProduct(cmatvolv, invCv, gamma * std::pow(detC, -kappa));

  // --------------------------------------------------------------
  // update constituent stress and material tangent
  stress->Update(1.0, S1v, 1.0);
  stress->Update(1.0, S2v, 1.0);
  stress->Update(1.0, Svolv, 1.0);

  cmat->Update(1.0, cmat1v, 1.0);
  cmat->Update(1.0, cmat2v, 1.0);
  cmat->Update(1.0, cmatvolv, 1.0);
}

CORE::UTILS::ValuesFunctAndFunctDerivs MAT::Muscle_Giantesio::EvaluateActivationLevelAndDerivatives(
    const double& lambdaM, const double& dotLambdaM, const double& currentTime)
{
  // setup function
  std::function<double(double)> FunctionSolveActivationLevelEquation =
      [this, &dotLambdaM, &currentTime](double lambda_i)
  { return this->SolveActivationLevelEquation(lambda_i, dotLambdaM, currentTime); };

  // step size for finite differences (error will be of O(h^2))
  // rule of thumb: h = x0*sqrt(ulp); but omegaa is only precise to 1e-12
  const double h = lambdaM * std::sqrt(1e-12);

  // compute activation level omegaa and its first and second derivative w.r.t. lambdaM using
  // central differences
  CORE::UTILS::ValuesFunctAndFunctDerivs omegaaAndDerivs =
      CORE::UTILS::EvaluateFunctionAndDerivativesCentralDifferences(
          FunctionSolveActivationLevelEquation, lambdaM, h);

  return omegaaAndDerivs;
}

double MAT::Muscle_Giantesio::SolveActivationLevelEquation(
    const double& lambdaM, const double& dotLambdaM, const double& currentTime)
{
  // compute right hand side of equation
  double rhs = EvaluateRhsActivationLevelEquation(lambdaM, dotLambdaM, currentTime);

  // setup the activation level equation as f(omegaa) = lhs-rhs
  // and its derivative as df/domegaa(omegaa) = dlhs/domegaa
  auto actLevelEquationAndDeriv = [this, &lambdaM, &rhs](double omegaa_init)
  { return this->EvaluateActivationLevelEquationAndDeriv(omegaa_init, lambdaM, rhs); };

  // determine starting guess for newton solver
  double omegaa_init;
  if (omegaaOld_ < 0.0)
  {
    // setup parameters for bisection method for the approximation of the starting guess for the
    // newton solver in case omegaa is not available from a previous timestep
    const double tol_bisec = 1e-2;  // abs(f(x)) <= 0.01
    const int maxiter_bisec = 100;
    const double omegaa_a_init = 0.0;  // omegaa >= 0
    const double omegaa_b_init = 1.0;  // omegaa <= 1

    // approximate the starting guess for the newton solver via bisection method
    omegaa_init = CORE::UTILS::Bisection([&](double omegaa_init)
        { return std::get<0>(actLevelEquationAndDeriv(omegaa_init)); },
        omegaa_a_init, omegaa_b_init, tol_bisec, maxiter_bisec);
  }
  else
    // use omegaa from the previous timestep as a starting guess
    omegaa_init = omegaaOld_;

  // setup parameters for newton method
  const double tol_newton = 1e-12;  // reasonably small
  const int maxiter_newton = 200;

  // compute activation level as solution of activation level equation
  double omegaa = CORE::UTILS::SolveLocalNewton(
      actLevelEquationAndDeriv, omegaa_init, tol_newton, maxiter_newton);

  return omegaa;
}

std::tuple<double, double> MAT::Muscle_Giantesio::EvaluateActivationLevelEquationAndDeriv(
    double omegaa, const double& lambdaM, const double& rhs)
{
  // get active microstructural parameters from params_
  const double alpha = params_->alpha_;
  const double beta = params_->beta_;
  const double omega0 = params_->omega0_;

  // elastic invariants Ie and Je
  double Ie = (-omega0 * (std::pow(omegaa - 1., 3) / std::pow(lambdaM, 3) + 1.) + 1.5) *
              std::pow(lambdaM, 2) / (1.5 * std::pow(omegaa - 1., 2));
  double Je = (-omega0 * (std::pow(lambdaM, 3) / std::pow(omegaa - 1., 3) + 1.) + 1.5) *
              std::pow(omegaa - 1., 2) / (1.5 * std::pow(lambdaM, 2));

  // derivatives of elastic invariants Ie and Je w.r.t. lambdaM
  double derivIe = (std::pow(lambdaM, 3) * (2 * omega0 - 3) - omega0 * std::pow(omegaa - 1., 3)) /
                   (1.5 * lambdaM * std::pow(omegaa - 1., 3));
  double derivJe = (std::pow(lambdaM, 3) * omega0 - (2 * omega0 - 3) * std::pow(omegaa - 1., 3)) /
                   (1.5 * std::pow(lambdaM, 2) * std::pow(omegaa - 1., 2));

  // compute left hand side of the activation level equation
  double lhs = std::exp(alpha * (Ie - 1)) / alpha + std::exp(beta * (Je - 1)) / beta;

  // compute derivative of left hand side of the activation level equation
  double derivlhs = derivIe * std::exp(alpha * (Ie - 1)) + derivJe * std::exp(beta * (Je - 1));

  // compute activation level equation and its derivative
  return {lhs - rhs, derivlhs};
}

double MAT::Muscle_Giantesio::EvaluateRhsActivationLevelEquation(
    const double& lambdaM, const double& dotLambdaM, const double& currentTime)
{
  // get active microstructural parameters from params_
  const double alpha = params_->alpha_;
  const double beta = params_->beta_;
  const double gamma = params_->gamma_;
  const double omega0 = params_->omega0_;

  // compute active nominaal stress and its integral w.r.t. lambdaM
  double intPa = EvaluateActiveNominalStressIntegral(lambdaM, dotLambdaM, currentTime);

  // passive part of invariants Ip and Jp w.r.t. lambdaM
  const double Ip = (omega0 / std::pow(lambdaM, 3) - omega0 + 1.5) * std::pow(lambdaM, 2) / 1.5;
  const double Jp = (omega0 * std::pow(lambdaM, 3) - omega0 + 1.5) / std::pow(lambdaM, 2) / 1.5;

  // compute right hand side of equation
  double rhs = 4.0 * intPa / gamma + std::exp(alpha * (Ip - 1.0)) / alpha +
               std::exp(beta * (Jp - 1.0)) / beta;

  return rhs;
}

double MAT::Muscle_Giantesio::EvaluateActiveNominalStressIntegral(
    const double& lambdaM, const double& dotLambdaM, const double& currentTime)
{
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
      Na, muTypesNum, rho, I, F, T, actIntervalsNum, actTimes, actValues, currentTime);

  // compute integral of the force-stretch dependency fxi in the boundaries lambdaMin to lambdaM
  double intFxi = MAT::UTILS::MUSCLE::EvaluateIntegralForceStretchDependencyEhret(
      lambdaM, lambdaMin, lambdaOpt);

  // compute force-velocity dependency fv
  double fv = MAT::UTILS::MUSCLE::EvaluateForceVelocityDependencyBoel(
      dotLambdaM, dotLambdaMMin, de, dc, ke, kc);

  // compute integral of the active nominal stress Pa
  double intPa = Poptft * intFxi * fv;

  return intPa;
}

CORE::LINALG::Matrix<3, 3> MAT::Muscle_Giantesio::ActDefGrad(
    const double omegaa, const CORE::LINALG::Matrix<3, 3>& M)
{
  // active deformation gradient Fa
  CORE::LINALG::Matrix<3, 3> Fa(M);
  Fa.Scale(1 - omegaa - std::pow(1. - omegaa, -0.5));
  for (unsigned j = 0; j < 3; ++j) Fa(j, j) += std::pow(1. - omegaa, -0.5);
  return Fa;
}

CORE::LINALG::Matrix<3, 3> MAT::Muscle_Giantesio::DActDefGrad_DActLevel(
    const double omegaa, const CORE::LINALG::Matrix<3, 3>& M)
{
  // first derivative of Fa w.r.t. omegaa
  CORE::LINALG::Matrix<3, 3> dFadomegaa(M);
  if (omegaa != 0)
  {
    dFadomegaa.Scale(-1.0 - 0.5 * std::pow(1. - omegaa, -1.5));
    for (unsigned j = 0; j < 3; ++j) dFadomegaa(j, j) += 0.5 * std::pow(1. - omegaa, -1.5);
  }
  else
  {
    dFadomegaa.Clear();
  }
  return dFadomegaa;
}

CORE::LINALG::Matrix<3, 3> MAT::Muscle_Giantesio::DDActDefGrad_DDActLevel(
    const double omegaa, const CORE::LINALG::Matrix<3, 3>& M)
{
  // second derivative of Fa w.r.t. omegaa
  CORE::LINALG::Matrix<3, 3> ddFaddomegaa(M);
  if (omegaa != 0)
  {
    ddFaddomegaa.Scale(-0.75 * std::pow(1. - omegaa, -2.5));
    for (unsigned j = 0; j < 3; ++j) ddFaddomegaa(j, j) += 0.75 * std::pow(1. - omegaa, -2.5);
  }
  else
  {
    ddFaddomegaa.Clear();
  }
  return ddFaddomegaa;
}

CORE::LINALG::Matrix<3, 3> MAT::Muscle_Giantesio::InvActDefGrad(
    const CORE::LINALG::Matrix<3, 3>& Fa)
{
  // inverse of active deformation gradient Fa^{-1}
  CORE::LINALG::Matrix<3, 3> invFa(true);
  invFa.Invert(Fa);
  return invFa;
}

CORE::LINALG::Matrix<3, 3> MAT::Muscle_Giantesio::DInvActDefGrad_DActLevel(
    const CORE::LINALG::Matrix<3, 3>& Fa, const CORE::LINALG::Matrix<3, 3>& dFadomegaa)
{
  CORE::LINALG::Matrix<3, 3> invFa = InvActDefGrad(Fa);

  // first derivative of Fa^{-1} w.r.t. omegaa
  // dinvFadomegaa_ij = - invFa_ik dFadomegaa_kl invFa_lj
  CORE::LINALG::Matrix<3, 3> dinvFadomegaa(true);
  MultiplyNNN(dinvFadomegaa, -1.0, invFa, dFadomegaa, invFa);

  return dinvFadomegaa;
}

bool MAT::Muscle_Giantesio::IsActive(const double& currentTime)
{
  bool isActive = false;

  int times_index = 0;
  while (currentTime >= params_->actTimes_[times_index]) times_index++;
  int values_index = times_index - 1;

  // values_index = -1 -> currentTime = 0 -> choose first value
  if (values_index == -1) values_index = 0;

  if (values_index >= 0)
  {
    if (params_->actValues_[values_index] > 1e-14) isActive = true;
  }

  return isActive;
};