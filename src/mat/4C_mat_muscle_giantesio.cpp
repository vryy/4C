/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of the Giantesio active strain skeletal muscle material (active strain
approach)

\level 3

*/
/*----------------------------------------------------------------------*/

#include "4C_mat_muscle_giantesio.hpp"

#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_linalg_four_tensor.hpp"
#include "4C_mat_muscle_utils.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"
#include "4C_matelast_aniso_structuraltensor_strategy.hpp"
#include "4C_utils_local_newton.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  /*!
   * @brief Returns the fourth order tensor
   *        result_ijkl = scalar * (tensor_A_ijkm * matrix_B_mn * matrix_C_nl
   *               + matrix_A_im * tensor_B_mjkn * matrix_C_nl
   *               + matrix_A_im * matrix_B_mn * tensor_C_njkl)
   */
  Core::LinAlg::FourTensor<3> SumMultiplyTmmMtmMmt(const Core::LinAlg::Matrix<3, 3>& matrix_A,
      const Core::LinAlg::Matrix<3, 3>& matrix_B, const Core::LinAlg::Matrix<3, 3>& matrix_C,
      const Core::LinAlg::FourTensor<3>& tensor_A, const Core::LinAlg::FourTensor<3>& tensor_B,
      const Core::LinAlg::FourTensor<3>& tensor_C, const double& scalar)
  {
    Core::LinAlg::Matrix<3, 3> mBmC(false);
    mBmC.MultiplyNN(matrix_B, matrix_C);

    Core::LinAlg::Matrix<3, 3> mAmB(false);
    mAmB.MultiplyNN(matrix_A, matrix_B);

    // part1 = tensor_A * matrix_B * matrix_C
    Core::LinAlg::FourTensor<3> part1(true);
    Mat::multiply_matrix_four_tensor_by_second_index<3>(part1, mBmC, tensor_A, false);

    // part2 = matrix_A * tensor_B * matrix_C
    Core::LinAlg::FourTensor<3> mAtB(true);  // matrix_A * tensor_B
    Mat::multiply_matrix_four_tensor<3>(mAtB, matrix_A, tensor_B, false);
    Core::LinAlg::FourTensor<3> part2(true);  // (matrix_A * tensor_B) * matrix_C
    Mat::multiply_matrix_four_tensor_by_second_index<3>(part2, matrix_C, mAtB, false);

    // part3 = matrix_A * matrix_B * tensor_C
    Core::LinAlg::FourTensor<3> part3(true);
    Mat::multiply_matrix_four_tensor<3>(part3, mAmB, tensor_C, false);

    // result = part1 + part2 + part3
    Core::LinAlg::FourTensor<3> tAmBmCmAtBmCmAmBtC(true);
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
  void MultiplyNNN(Core::LinAlg::Matrix<3, 3>& NNN, const double& scalar,
      const Core::LinAlg::Matrix<3, 3>& left, const Core::LinAlg::Matrix<3, 3>& middle,
      const Core::LinAlg::Matrix<3, 3>& right)
  {
    Core::LinAlg::Matrix<3, 3> NN(true);
    NN.MultiplyNN(left, middle);
    NNN.MultiplyNN(scalar, NN, right, 1.0);
  }

  /*!
   * @brief Compute the derivative of the Cauchy Green tensor C w.r.t. C
   */
  Core::LinAlg::FourTensor<3> ComputeDCDC()
  {
    Core::LinAlg::Matrix<6, 6> dCdCv(true);
    for (int i = 0; i < 3; i++) dCdCv(i, i) = 1.0;
    for (int i = 3; i < 6; i++) dCdCv(i, i) = 0.5;

    Core::LinAlg::FourTensor<3> dCdC(true);
    Mat::setup_four_tensor_from_6x6_voigt_matrix(dCdC, dCdCv);

    return dCdC;
  }

  /*!
   * @brief Compute the structural tensor L
   */
  Core::LinAlg::Matrix<3, 3> ComputeStructuralTensorL(
      const Core::LinAlg::Matrix<3, 3>& M, const double& omega0)
  {
    Core::LinAlg::Matrix<3, 3> L(M);
    L.Scale(1.0 - omega0);  // omegap*M
    for (unsigned i = 0; i < 3; ++i) L(i, i) += omega0 / 3.0;

    return L;
  }

  struct StressAndDeriv
  {
    Core::LinAlg::Matrix<3, 3> Se;
    Core::LinAlg::FourTensor<3> dSedC;
  };

  /*!
   * @brief Compute the elastic second Piola Kirchhoff stress tensor Se and its derivative dSedC
   * w.r.t. the right Cauchy Green tensor C
   */
  StressAndDeriv ComputeSeAndDSedC(const double& alpha, const double& beta, const double& gamma,
      const double& omega0, const Core::LinAlg::Matrix<3, 3>& M,
      const Core::LinAlg::Matrix<3, 3>& invFa, const Core::LinAlg::FourTensor<3>& dinvFadC,
      const Core::LinAlg::Matrix<3, 3>& C, const Core::LinAlg::FourTensor<3>& dCdC)
  {
    // structural tensor L = omega0/3*Identity + omegap*M
    Core::LinAlg::Matrix<3, 3> L = ComputeStructuralTensorL(M, omega0);
    Core::LinAlg::Matrix<6, 1> Lv(false);  // Voigt notation
    Core::LinAlg::Voigt::Stresses::matrix_to_vector(L, Lv);

    // elastic right Cauchy Green tensor Ce = Fe^T Fe
    // = Fa^-T C Fa^-1 = Fa^-1 C Fa^-1 (with Fa^-1 sym.)
    Core::LinAlg::Matrix<3, 3> Ce(true);
    MultiplyNNN(Ce, 1.0, invFa, C, invFa);

    // derivative of Ce w.r.t C
    // dCedC_ijkl = dinvFadC_iakl (C invFa)_aj + invFa_ia dCdC_abkl invFa_bj + (invFa C)_ab
    // dinvFadC_bjkl
    Core::LinAlg::FourTensor<3> dCedC(true);
    dCedC = SumMultiplyTmmMtmMmt(invFa, C, invFa, dinvFadC, dCdC, dinvFadC, 1.0);
    Core::LinAlg::Matrix<6, 6> dCedCv(true);
    Mat::setup_6x6_voigt_matrix_from_four_tensor(dCedCv, dCedC);

    // inverse of the elastic right Cauchy Green tensor Ce
    Core::LinAlg::Matrix<3, 3> invCe(true);
    invCe.Invert(Ce);
    Core::LinAlg::Matrix<6, 1> invCev(true);
    Core::LinAlg::Voigt::Stresses::matrix_to_vector(invCe, invCev);

    // product invCe*L
    Core::LinAlg::Matrix<3, 3> invCeL(true);
    invCeL.MultiplyNN(invCe, L);

    // product invCe*L*invCe
    Core::LinAlg::Matrix<3, 3> invCeLinvCe(true);
    invCeLinvCe.MultiplyNN(invCeL, invCe);
    Core::LinAlg::Matrix<6, 1> invCeLinvCev(true);
    Core::LinAlg::Voigt::Stresses::matrix_to_vector(invCeLinvCe, invCeLinvCev);

    // product Ce^T*L
    Core::LinAlg::Matrix<3, 3> transpCeL(true);
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
    Core::LinAlg::Matrix<3, 3> Se(true);
    Se.Update(expalpha, L);
    Se.Update(-expbeta * detCe, invCeLinvCe, 1.0);
    Se.Update(Je * expbeta, invCe, 1.0);
    Se.Scale(0.5 * gamma);

    // derivative of Se w.r.t Ce
    // taken from Weickenmeier et al.
    Core::LinAlg::Matrix<6, 6> dSedCev(true);
    dSedCev.MultiplyNT(alpha * expalpha, Lv, Lv, 1.0);  // add contributions
    dSedCev.MultiplyNT(beta * expbeta * std::pow(detCe, 2), invCeLinvCev, invCeLinvCev, 1.0);
    dSedCev.MultiplyNT(-(beta * Je + 1.) * expbeta * detCe, invCev, invCeLinvCev, 1.0);
    dSedCev.MultiplyNT(-(beta * Je + 1.) * expbeta * detCe, invCeLinvCev, invCev, 1.0);
    dSedCev.MultiplyNT((beta * Je + 1.) * Je * expbeta, invCev, invCev, 1.0);
    // adds scalar * (invC boeppel invC) to cmat, see Holzapfel2000, p. 254
    Mat::add_holzapfel_product(dSedCev, invCev, -Je * expbeta);
    // adds -expbeta * detCe * dinvCLinvCdCv to cmat
    Mat::add_derivative_of_inva_b_inva_product(-expbeta * detCe, invCev, invCeLinvCev, dSedCev);
    dSedCev.Scale(gamma / 2);

    Core::LinAlg::FourTensor<3> dSedCe(true);
    Mat::setup_four_tensor_from_6x6_voigt_matrix(dSedCe, dSedCev);

    // derivative of Se w.r.t C
    // dSedC_ijkl = dSedCe_ijab dCedC_abkl
    Core::LinAlg::FourTensor<3> dSedC(true);
    Mat::multiply_four_tensor_four_tensor<3>(dSedC, dSedCe, dCedC, true);

    StressAndDeriv Se_dSedC = {.Se = Se, .dSedC = dSedC};

    return Se_dSedC;
  }
}  // namespace

Mat::PAR::MuscleGiantesio::MuscleGiantesio(const Core::Mat::PAR::Parameter::Data& matdata)
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

Teuchos::RCP<Core::Mat::Material> Mat::PAR::MuscleGiantesio::create_material()
{
  return Teuchos::rcp(new Mat::MuscleGiantesio(this));
}

Mat::MuscleGiantesioType Mat::MuscleGiantesioType::instance_;

Core::Communication::ParObject* Mat::MuscleGiantesioType::Create(const std::vector<char>& data)
{
  auto* muscle_giantesio = new Mat::MuscleGiantesio();
  muscle_giantesio->unpack(data);
  return muscle_giantesio;
}

Mat::MuscleGiantesio::MuscleGiantesio()
    : params_(nullptr),
      lambda_m_old_(1.0),
      omegaa_old_(-1.0),
      anisotropy_(),
      anisotropy_extension_(true, 0.0, 0,
          Teuchos::rcp<Mat::Elastic::StructuralTensorStrategyBase>(
              new Mat::Elastic::StructuralTensorStrategyStandard(nullptr)),
          {0})
{
}

Mat::MuscleGiantesio::MuscleGiantesio(Mat::PAR::MuscleGiantesio* params)
    : params_(params),
      lambda_m_old_(1.0),
      omegaa_old_(-1.0),
      anisotropy_(),
      anisotropy_extension_(true, 0.0, 0,
          Teuchos::rcp<Mat::Elastic::StructuralTensorStrategyBase>(
              new Mat::Elastic::StructuralTensorStrategyStandard(nullptr)),
          {0})
{
  // initialize lambdaMOld_ and omegaOld_
  lambda_m_old_ = 1.0;
  omegaa_old_ = -1.0;

  // register anisotropy extension to global anisotropy
  anisotropy_.register_anisotropy_extension(anisotropy_extension_);

  // initialize fiber directions and structural tensor
  anisotropy_extension_.register_needed_tensors(
      Mat::FiberAnisotropyExtension<1>::FIBER_VECTORS |
      Mat::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR);
}

void Mat::MuscleGiantesio::pack(Core::Communication::PackBuffer& data) const
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
  add_to_pack(data, omegaa_old_);

  anisotropy_extension_.pack_anisotropy(data);
}

void Mat::MuscleGiantesio::unpack(const std::vector<char>& data)
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
        params_ = static_cast<Mat::PAR::MuscleGiantesio*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  extract_from_pack(position, data, lambda_m_old_);
  extract_from_pack(position, data, omegaa_old_);

  anisotropy_extension_.unpack_anisotropy(data, position);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

void Mat::MuscleGiantesio::setup(int numgp, Input::LineDefinition* linedef)
{
  // Read anisotropy
  anisotropy_.set_number_of_gauss_points(numgp);
  anisotropy_.read_anisotropy_from_element(linedef);
}

void Mat::MuscleGiantesio::Update(Core::LinAlg::Matrix<3, 3> const& defgrd, int const gp,
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

void Mat::MuscleGiantesio::evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 1>* stress, Core::LinAlg::Matrix<6, 6>* cmat, const int gp,
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
    FOUR_C_THROW("No total time given for muscle Giantesio material!");

  // save (time) step size
  double timeStepSize = params.get<double>("delta time", -1);
  if (std::abs(timeStepSize + 1.0) < 1e-14)
    FOUR_C_THROW("No time step size given for muscle Giantesio material!");

  // compute matrices
  // right Cauchy Green tensor C
  Core::LinAlg::Matrix<3, 3> C(false);  // matrix notation
  C.MultiplyTN(*defgrd, *defgrd);       // C = F^T F

  // determinant of C
  double detC = C.Determinant();

  // derivative of C w.r.t C
  Core::LinAlg::FourTensor<3> dCdC = ComputeDCDC();

  // inverse right Cauchy Green tensor C^-1
  Core::LinAlg::Matrix<3, 3> invC(true);                         // matrix notation
  invC.Invert(C);                                                // invC = C^-1
  Core::LinAlg::Matrix<6, 1> invCv(true);                        // Voigt notation
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(invC, invCv);  // invCv

  // structural tensor M, i.e. dyadic product of fibre directions
  Core::LinAlg::Matrix<3, 3> M = anisotropy_extension_.get_structural_tensor(gp, 0);
  double lambdaM = Mat::UTILS::Muscle::FiberStretch(C, M);
  // derivative of lambdaM w.r.t. C
  Core::LinAlg::Matrix<3, 3> dlambdaMdC = Mat::UTILS::Muscle::DFiberStretch_DC(lambdaM, C, M);
  Core::LinAlg::Matrix<6, 1> dlambdaMdCv(true);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(dlambdaMdC, dlambdaMdCv);

  // contraction velocity dotLambdaM
  double dotLambdaM =
      Mat::UTILS::Muscle::ContractionVelocityBWEuler(lambdaM, lambda_m_old_, timeStepSize);

  // compute activation level omegaa and derivative w.r.t. fiber stretch if material is active
  Core::UTILS::ValuesFunctAndFunctDerivs omegaaAndDerivs = {
      .val_funct = 0.0, .val_deriv_funct = 0.0, .val_deriv_deriv_funct = 0.0};

  if (is_active(currentTime))
  {
    omegaaAndDerivs = evaluate_activation_level_and_derivatives(lambdaM, dotLambdaM, currentTime);
    omegaa_old_ = omegaaAndDerivs.val_funct;  // update omegaaOld
  }
  else
  {
    omegaa_old_ = -1;  // reset omegaaOld to -1
  }

  // --------------------------------------------------------------
  // first derivative of omegaa w.r.t. C
  // domegaadC = domegaadlambdaM dlambdaMdC
  Core::LinAlg::Matrix<3, 3> domegaadC(dlambdaMdC);
  domegaadC.Scale(omegaaAndDerivs.val_deriv_funct);
  Core::LinAlg::Matrix<6, 1> domegaadCv(false);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(domegaadC, domegaadCv);

  // second derivative of omegaa w.r.t. C
  // ddomegaaddC = (ddomegaaddlambdaM - 1/lambdaM domegaadlambdaM) dlambdaM_ij dlambdaM_kl
  Core::LinAlg::Matrix<6, 6> ddomegaaddCv(false);
  ddomegaaddCv.MultiplyNT(
      omegaaAndDerivs.val_deriv_deriv_funct - omegaaAndDerivs.val_deriv_funct / lambdaM,
      dlambdaMdCv, dlambdaMdCv);

  // --------------------------------------------------------------
  // active deformation gradient Fa
  Core::LinAlg::Matrix<3, 3> Fa = act_def_grad(omegaaAndDerivs.val_funct, M);

  // first derivative of Fa w.r.t. omegaa
  Core::LinAlg::Matrix<3, 3> dFadomegaa = d_act_def_grad_d_act_level(omegaaAndDerivs.val_funct, M);

  // second derivative of Fa w.r.t. omegaa
  Core::LinAlg::Matrix<3, 3> ddFaddomegaa =
      dd_act_def_grad_dd_act_level(omegaaAndDerivs.val_funct, M);

  // determinant of Fa
  double detFa = Fa.Determinant();

  // --------------------------------------------------------------
  // inverse of active deformation gradient Fa^{-1}
  Core::LinAlg::Matrix<3, 3> invFa = inv_act_def_grad(Fa);

  // first derivative of Fa^{-1} w.r.t. omegaa
  Core::LinAlg::Matrix<3, 3> dinvFadomegaa = d_inv_act_def_grad_d_act_level(Fa, dFadomegaa);

  // first derivative of Fa^{-1} w.r.t. C
  // dinvFadC_ijkl = dinvFadomegaa_ij domegaadC_kl
  Core::LinAlg::FourTensor<3> dinvFadC(true);
  Mat::add_dyadic_product_matrix_matrix(dinvFadC, dinvFadomegaa, domegaadC);

  // --------------------------------------------------------------
  // elastic second Piola Kirchhoff stress tensor Se and its derivative w.r.t C
  StressAndDeriv Se_dSedC =
      ComputeSeAndDSedC(alpha, beta, gamma, omega0, M, invFa, dinvFadC, C, dCdC);

  // --------------------------------------------------------------
  // compute second Piola Kirchhoff stress components S1, S2, Svol
  // note: S3 = invF * ddetFadF * Psie = 0 since ddetFawa = 0 and ddetFadF = ddetFadwa dwadF

  // S1 = detFa * invFa * Se * invFa^T
  // in index notation: S1_sl = detFa * invFa_sr * Se_rq * invFa_rl (with invFa sym.)
  Core::LinAlg::Matrix<3, 3> S1(true);
  MultiplyNNN(S1, detFa, invFa, Se_dSedC.Se, invFa);
  Core::LinAlg::Matrix<6, 1> S1v(true);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(S1, S1v);

  // S2 = detFa * invF * dPsiedFe * (F * dinvFadF)  = 2 detFa (C Fa^-1 Se) dFa^-1dC
  // or: S2 = -2 S1 : (C Fa^-1 dFadomegaa) domegaadC
  // in index notation: S2_sl = -2 S1_ij (C Fa^-1 dFadomegaa)_ij domegaadC_sl
  Core::LinAlg::Matrix<3, 3> CinvFadFadomegaa(true);
  MultiplyNNN(CinvFadFadomegaa, 1.0, C, invFa, dFadomegaa);

  Core::LinAlg::Matrix<6, 1> S2v(domegaadCv);
  double scalar_contraction = -2.0 * Mat::contract_matrix_matrix(CinvFadFadomegaa, S1);
  S2v.Scale(scalar_contraction);

  // Svol = - gamma/2 * detC^-kappa
  Core::LinAlg::Matrix<6, 1> Svolv(invCv);
  Svolv.Scale(-0.5 * gamma * std ::pow(detC, -kappa));

  // --------------------------------------------------------------
  // compute material tangent components cmat1, cmat2, cmatvol

  // compute cmat1 = 2 dS1dC
  Core::LinAlg::Matrix<6, 6> cmat1v(true);

  // derivative of S1 = detFa * invFa * Se * invFa^T  w.r.t. C
  // = detFa * (dinvFadC * (Se * dinvFa^T) + invFa * dSedC * invFa^T + (invFa * Se) * dinvFadC)
  Core::LinAlg::FourTensor<3> dS1dC(true);
  dS1dC =
      SumMultiplyTmmMtmMmt(invFa, Se_dSedC.Se, invFa, dinvFadC, Se_dSedC.dSedC, dinvFadC, detFa);
  Mat::setup_6x6_voigt_matrix_from_four_tensor(cmat1v, dS1dC);
  cmat1v.Scale(2.0);

  // compute cmat2 = 2 dS2dC
  Core::LinAlg::Matrix<6, 6> cmat2v(true);

  // derivative of S2 w.r.t. C
  // in index notation: S2_sl = -2 S1_ij (C Fa^-1 dFadomegaa)_ij domegaadC_sl
  // dS2dC_slpq = -4 [ S1_ij d(C Fa^-1 dFadomegaa)dC_ijpq domegaadC_sl +
  //                 [ S1_ij (C Fa^-1 dFadomegaa)_ij ddomegaaddC_slpq ]
  // dS2dC_slpq =  2 [ dscalardC_pq domegaadC_sl + [ scalar_contraction ddomegaaddC_slpq ]

  // derivative of the scalar (-2 S1_ij (C Fa^-1 dFadomegaa)_ij) w.r.t. C
  // dscalardC_pq = -2 (dS1dC_ijpq (C invFa dFadomegaa)_ij + S1_ij d(C invFa dFadomegaa)dC_ijpq )
  Core::LinAlg::Matrix<3, 3> H(true);
  H.MultiplyNN(dinvFadomegaa, dFadomegaa);
  H.MultiplyNN(1.0, invFa, ddFaddomegaa, 1.0);
  Core::LinAlg::Matrix<3, 3> CH(true);
  CH.MultiplyNN(C, H);
  double helper = Mat::contract_matrix_matrix(S1, CH);

  Core::LinAlg::Matrix<3, 3> dscalardC(domegaadC);
  dscalardC.Scale(helper);
  MultiplyNNN(dscalardC, 1.0, S1, dFadomegaa, invFa);
  Mat::add_contraction_matrix_four_tensor(dscalardC, CinvFadFadomegaa, dS1dC);
  Core::LinAlg::Matrix<6, 1> dscalardCv(true);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(
      dscalardC, dscalardCv);  // not yet scaled with -2.0

  // compute cmat2 = 2 dS2dC = 2 [domegaadC_sl dscalardC_pq + scalar * ddomegaaddC_slpq]
  cmat2v.MultiplyNT(-4.0, domegaadCv, dscalardCv);             // 2 * dscalardC_pq domegaadC_sl
  cmat2v.Update(2.0 * scalar_contraction, ddomegaaddCv, 1.0);  // 2 * scalar * ddomegaaddC_slpq

  // compute cmatvol = 2 dS2dC
  Core::LinAlg::Matrix<6, 6> cmatvolv(true);
  cmatvolv.MultiplyNT(gamma * kappa * std::pow(detC, -kappa), invCv, invCv, 1.0);
  Mat::add_holzapfel_product(cmatvolv, invCv, gamma * std::pow(detC, -kappa));

  // --------------------------------------------------------------
  // update constituent stress and material tangent
  stress->Update(1.0, S1v, 1.0);
  stress->Update(1.0, S2v, 1.0);
  stress->Update(1.0, Svolv, 1.0);

  cmat->Update(1.0, cmat1v, 1.0);
  cmat->Update(1.0, cmat2v, 1.0);
  cmat->Update(1.0, cmatvolv, 1.0);
}

Core::UTILS::ValuesFunctAndFunctDerivs
Mat::MuscleGiantesio::evaluate_activation_level_and_derivatives(
    const double& lambdaM, const double& dotLambdaM, const double& currentTime)
{
  // setup function
  std::function<double(double)> FunctionSolveActivationLevelEquation =
      [this, &dotLambdaM, &currentTime](double lambda_i)
  { return this->solve_activation_level_equation(lambda_i, dotLambdaM, currentTime); };

  // step size for finite differences (error will be of O(h^2))
  // rule of thumb: h = x0*sqrt(ulp); but omegaa is only precise to 1e-12
  const double h = lambdaM * std::sqrt(1e-12);

  // compute activation level omegaa and its first and second derivative w.r.t. lambdaM using
  // central differences
  Core::UTILS::ValuesFunctAndFunctDerivs omegaaAndDerivs =
      Core::UTILS::evaluate_function_and_derivatives_central_differences(
          FunctionSolveActivationLevelEquation, lambdaM, h);

  return omegaaAndDerivs;
}

double Mat::MuscleGiantesio::solve_activation_level_equation(
    const double& lambdaM, const double& dotLambdaM, const double& currentTime)
{
  // compute right hand side of equation
  double rhs = evaluate_rhs_activation_level_equation(lambdaM, dotLambdaM, currentTime);

  // setup the activation level equation as f(omegaa) = lhs-rhs
  // and its derivative as df/domegaa(omegaa) = dlhs/domegaa
  auto actLevelEquationAndDeriv = [this, &lambdaM, &rhs](double omegaa_init)
  { return this->evaluate_activation_level_equation_and_deriv(omegaa_init, lambdaM, rhs); };

  // determine starting guess for newton solver
  double omegaa_init;
  if (omegaa_old_ < 0.0)
  {
    // setup parameters for bisection method for the approximation of the starting guess for the
    // newton solver in case omegaa is not available from a previous timestep
    const double tol_bisec = 1e-2;  // abs(f(x)) <= 0.01
    const int maxiter_bisec = 100;
    const double omegaa_a_init = 0.0;  // omegaa >= 0
    const double omegaa_b_init = 1.0;  // omegaa <= 1

    // approximate the starting guess for the newton solver via bisection method
    omegaa_init = Core::UTILS::bisection([&](double omegaa_init)
        { return std::get<0>(actLevelEquationAndDeriv(omegaa_init)); },
        omegaa_a_init, omegaa_b_init, tol_bisec, maxiter_bisec);
  }
  else
    // use omegaa from the previous timestep as a starting guess
    omegaa_init = omegaa_old_;

  // setup parameters for newton method
  const double tol_newton = 1e-12;  // reasonably small
  const int maxiter_newton = 200;

  // compute activation level as solution of activation level equation
  double omegaa = Core::UTILS::solve_local_newton(
      actLevelEquationAndDeriv, omegaa_init, tol_newton, maxiter_newton);

  return omegaa;
}

std::tuple<double, double> Mat::MuscleGiantesio::evaluate_activation_level_equation_and_deriv(
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

double Mat::MuscleGiantesio::evaluate_rhs_activation_level_equation(
    const double& lambdaM, const double& dotLambdaM, const double& currentTime)
{
  // get active microstructural parameters from params_
  const double alpha = params_->alpha_;
  const double beta = params_->beta_;
  const double gamma = params_->gamma_;
  const double omega0 = params_->omega0_;

  // compute active nominaal stress and its integral w.r.t. lambdaM
  double intPa = evaluate_active_nominal_stress_integral(lambdaM, dotLambdaM, currentTime);

  // passive part of invariants Ip and Jp w.r.t. lambdaM
  const double Ip = (omega0 / std::pow(lambdaM, 3) - omega0 + 1.5) * std::pow(lambdaM, 2) / 1.5;
  const double Jp = (omega0 * std::pow(lambdaM, 3) - omega0 + 1.5) / std::pow(lambdaM, 2) / 1.5;

  // compute right hand side of equation
  double rhs = 4.0 * intPa / gamma + std::exp(alpha * (Ip - 1.0)) / alpha +
               std::exp(beta * (Jp - 1.0)) / beta;

  return rhs;
}

double Mat::MuscleGiantesio::evaluate_active_nominal_stress_integral(
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
  double Poptft = Mat::UTILS::Muscle::EvaluateTimeDependentActiveStressEhret(
      Na, muTypesNum, rho, I, F, T, actIntervalsNum, actTimes, actValues, currentTime);

  // compute integral of the force-stretch dependency fxi in the boundaries lambdaMin to lambdaM
  double intFxi = Mat::UTILS::Muscle::EvaluateIntegralForceStretchDependencyEhret(
      lambdaM, lambdaMin, lambdaOpt);

  // compute force-velocity dependency fv
  double fv = Mat::UTILS::Muscle::EvaluateForceVelocityDependencyBoel(
      dotLambdaM, dotLambdaMMin, de, dc, ke, kc);

  // compute integral of the active nominal stress Pa
  double intPa = Poptft * intFxi * fv;

  return intPa;
}

Core::LinAlg::Matrix<3, 3> Mat::MuscleGiantesio::act_def_grad(
    const double omegaa, const Core::LinAlg::Matrix<3, 3>& M)
{
  // active deformation gradient Fa
  Core::LinAlg::Matrix<3, 3> Fa(M);
  Fa.Scale(1 - omegaa - std::pow(1. - omegaa, -0.5));
  for (unsigned j = 0; j < 3; ++j) Fa(j, j) += std::pow(1. - omegaa, -0.5);
  return Fa;
}

Core::LinAlg::Matrix<3, 3> Mat::MuscleGiantesio::d_act_def_grad_d_act_level(
    const double omegaa, const Core::LinAlg::Matrix<3, 3>& M)
{
  // first derivative of Fa w.r.t. omegaa
  Core::LinAlg::Matrix<3, 3> dFadomegaa(M);
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

Core::LinAlg::Matrix<3, 3> Mat::MuscleGiantesio::dd_act_def_grad_dd_act_level(
    const double omegaa, const Core::LinAlg::Matrix<3, 3>& M)
{
  // second derivative of Fa w.r.t. omegaa
  Core::LinAlg::Matrix<3, 3> ddFaddomegaa(M);
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

Core::LinAlg::Matrix<3, 3> Mat::MuscleGiantesio::inv_act_def_grad(
    const Core::LinAlg::Matrix<3, 3>& Fa)
{
  // inverse of active deformation gradient Fa^{-1}
  Core::LinAlg::Matrix<3, 3> invFa(true);
  invFa.Invert(Fa);
  return invFa;
}

Core::LinAlg::Matrix<3, 3> Mat::MuscleGiantesio::d_inv_act_def_grad_d_act_level(
    const Core::LinAlg::Matrix<3, 3>& Fa, const Core::LinAlg::Matrix<3, 3>& dFadomegaa)
{
  Core::LinAlg::Matrix<3, 3> invFa = inv_act_def_grad(Fa);

  // first derivative of Fa^{-1} w.r.t. omegaa
  // dinvFadomegaa_ij = - invFa_ik dFadomegaa_kl invFa_lj
  Core::LinAlg::Matrix<3, 3> dinvFadomegaa(true);
  MultiplyNNN(dinvFadomegaa, -1.0, invFa, dFadomegaa, invFa);

  return dinvFadomegaa;
}

bool Mat::MuscleGiantesio::is_active(const double& currentTime)
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
FOUR_C_NAMESPACE_CLOSE
