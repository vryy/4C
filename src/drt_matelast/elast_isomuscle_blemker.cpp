/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the isochoric part of the Blemker active skeletal muscle material

\level 3
*/
/*----------------------------------------------------------------------*/

#include "elast_isomuscle_blemker.H"
#include "elast_aniso_structuraltensor_strategy.H"
#include "drt_globalproblem.H"
#include "standardtypes_cpp.H"
#include "drt_linedefinition.H"
#include "voigt_notation.H"
#include "elasthyper_service.H"
#include "material_service.H"
#include "matpar_material.H"


MAT::ELASTIC::PAR::IsoMuscleBlemker::IsoMuscleBlemker(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata),
      G1_(matdata->GetDouble("G1")),
      G2_(matdata->GetDouble("G2")),
      P1_(matdata->GetDouble("P1")),
      P2_(matdata->GetDouble("P2")),
      sigma_max_(matdata->GetDouble("SIGMAMAX")),
      lambda_ofl_(matdata->GetDouble("LAMBDAOFL")),
      lambda_star_(matdata->GetDouble("LAMBDASTAR")),
      alpha_(matdata->GetDouble("ALPHA")),
      beta_(matdata->GetDouble("BETA")),
      t_act_start_(matdata->GetDouble("ACTSTARTTIME"))
{
  // error handling for parameter ranges
  if (G1_ < 0.0) dserror("Material parameter G1 must be positive or zero");
  if (G2_ < 0.0) dserror("Material parameter G2 must be positive or zero");
  if (P1_ <= 0.0) dserror("Material parameter P1 must be greater zero");
  if (P2_ <= 0.0) dserror("Material parameter P2 must be greater zero");
  if (sigma_max_ < 0.0) dserror("Material parameter SIGMA_MAX must be positive or zero");
  if (lambda_ofl_ <= 0.0) dserror("Material parameter LAMBDA_OFL must be greater zero");
  if (lambda_star_ <= 0.0) dserror("Material parameter LAMBDA_STAR must be greater zero");
  if (alpha_ < 0.0) dserror("Material parameter ALPHA must be positive or zero");
  if (beta_ < 0.0) dserror("Material parameter BETA must be positive or zero");
}


MAT::ELASTIC::IsoMuscleBlemker::IsoMuscleBlemker(MAT::ELASTIC::PAR::IsoMuscleBlemker* params)
    : params_(params),
      anisotropyExtension_(true, 0.0, false,
          Teuchos::rcp<MAT::ELASTIC::StructuralTensorStrategyBase>(
              new MAT::ELASTIC::StructuralTensorStrategyStandard(nullptr)),
          {0})

{
  // initialize fiber directions and structural tensor
  anisotropyExtension_.RegisterNeededTensors(
      MAT::FiberAnisotropyExtension<1>::FIBER_VECTORS |
      MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR |
      MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR_STRESS);
}

void MAT::ELASTIC::IsoMuscleBlemker::PackSummand(DRT::PackBuffer& data) const
{
  anisotropyExtension_.PackAnisotropy(data);
}

void MAT::ELASTIC::IsoMuscleBlemker::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  anisotropyExtension_.UnpackAnisotropy(data, position);
}

void MAT::ELASTIC::IsoMuscleBlemker::RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy)
{
  anisotropy.RegisterAnisotropyExtension(anisotropyExtension_);
}

void MAT::ELASTIC::IsoMuscleBlemker::AddStressAnisoModified(const LINALG::Matrix<6, 1>& rcg,
    const LINALG::Matrix<6, 1>& icg, LINALG::Matrix<6, 6>& cmat, LINALG::Matrix<6, 1>& stress,
    double I3, const int gp, const int eleGID, Teuchos::ParameterList& params)
{
  // right Cauchy Green tensor C in matrix notation
  LINALG::Matrix<3, 3> C(true);
  UTILS::VOIGT::Strains::VectorToMatrix(rcg, C);

  // compute volume ratio J
  double J = std::sqrt(I3);
  double incJ = std::pow(J, -2.0 / 3.0);

  // compute modified right Cauchy-Green Tensor
  LINALG::Matrix<3, 3> modC(false);
  modC.Update(incJ, C, 0.0);

  // structural tensor M, i.e. dyadic product of fibre directions
  LINALG::Matrix<3, 3> M = anisotropyExtension_.GetStructuralTensor(gp, 0);
  LINALG::Matrix<6, 1> Mv = anisotropyExtension_.GetStructuralTensor_stress(gp, 0);

  // compute modified invariants modI1, modI4 and modI5
  double modI1 = modC(0, 0) + modC(1, 1) + modC(2, 2);

  LINALG::Matrix<3, 3> modCM(true);
  modCM.MultiplyNN(1.0, modC, M);                            // modC*M
  double modI4 = (modCM(0, 0) + modCM(1, 1) + modCM(2, 2));  // modI4 = tr(modC*M)

  LINALG::Matrix<3, 3> modC2(true);
  modC2.MultiplyNN(1.0, modC, modC);  // modC*modC
  LINALG::Matrix<3, 3> modC2M(true);
  modC2M.MultiplyNN(1.0, modC2, M);                             // (modC*modC)*M
  double modI5 = (modC2M(0, 0) + modC2M(1, 1) + modC2M(2, 2));  // modI5 = tr((modC*modC)*M)

  // compute time-dependent activation level
  const double& alpha = params_->alpha_;
  const double& beta = params_->beta_;
  const double& t_tot = params.get<double>("total time");  // current simulation time
  double alpha_act = EvaluateTimeDependentActivationLevel(alpha, beta, t_tot);

  // compute total fiber cauchy stress and derivative w.r.t. fibre stretch
  double sigma_fiber_total = 0.0;
  double deriv_sigma_fiber_total = 0.0;
  EvaluateTotalFiberCauchyStressAndDerivative(
      modI4, alpha_act, sigma_fiber_total, deriv_sigma_fiber_total);

  // helper variables for computation of 2nd Piola Kirchhoff stress and elasticity tensor
  double H1 = (modI1 * modI4 - modI5) / (2.0 * std::sqrt(modI4));
  // prevents singularitys in cross-fiber-shear-free states
  if ((H1 - 1.0) < EPS15) H1 = 1.0 + EPS15;
  double H2 = std::sqrt(H1 * H1 - 1.0);
  double H3 = modI1 / (2.0 * std::sqrt(modI4)) - H1 / (2.0 * modI4);
  double B2 = std::acosh(H1);

  // get material parameters
  const double& G1 = params_->G1_;
  const double& G2 = params_->G2_;

  // scalar prefactors, gamma_i = derivative of strain-energy function w.r.t. modified
  // invariant i; Psi_iso = W1(modI4, modI5) + W2(mod_I1, modI4, modI5) + W3(modI4)
  double gamma1 = 2.0 * G2 * (B2 / H2) * std::sqrt(modI4);
  double gamma4_1 = -4.0 * G1 * modI5 / std::pow(modI4, 3.0);  // dW1/dÍ4
  double gamma4_2 = 4.0 * G2 * (B2 / H2) * H3;                 // dW2/dÍ4
  double gamma4_3 = sigma_fiber_total / modI4;                 // dW3/dÍ4
  double gamma4 = gamma4_1 + gamma4_2 + gamma4_3;
  double gamma5 = 2.0 * G1 / std::pow(modI4, 2.0) - 2.0 * G2 * (B2 / H2) / std::sqrt(modI4);

  // matrix terms for stress evaluation
  // unitary 3x3 matrix
  LINALG::Matrix<3, 3> Id3(false);
  MAT::IdentityMatrix(Id3);
  LINALG::Matrix<6, 1> Id3v(false);
  UTILS::VOIGT::IdentityMatrix(Id3v);

  // sum modC*M + M*modC = dI5/dC
  LINALG::Matrix<3, 3> modCMsumMmodC(modCM);
  modCMsumMmodC.MultiplyNN(1.0, M, modC, 1.0);
  LINALG::Matrix<6, 1> modCMsumMmodCv(false);
  UTILS::VOIGT::Stresses::MatrixToVector(modCMsumMmodC, modCMsumMmodCv);

  // ficticious 2nd Piola-Kirchhoff stress tensor modS
  LINALG::Matrix<3, 3> modS(true);
  modS.Update(gamma1, Id3, 1.0);            // + gamma1*I
  modS.Update(gamma4, M, 1.0);              // + gamma4*M
  modS.Update(gamma5, modCMsumMmodC, 1.0);  // dyad(a0,modC*a0)+dyad(a0*C,a0) = (modC*M)' +modC'*M
  LINALG::Matrix<6, 1> modSv(false);
  UTILS::VOIGT::Stresses::MatrixToVector(modS, modSv);

  // isometirc 2nd Piola-Kirchhoff tensor S_iso from ficticious 2nd PK stress
  double traceCmodS = modSv(0) * rcg(0) + modSv(1) * rcg(1) + modSv(2) * rcg(2) +
                      (modSv(3) * rcg(3) + modSv(4) * rcg(4) + modSv(5) * rcg(5));
  LINALG::Matrix<6, 1> S_isov(modSv);
  S_isov.Update((-incJ / 3.0) * traceCmodS, icg, incJ);  // icg comes in stress-like notation

  // update 2nd Piola-Kirchhoff  tensor
  stress.Update(1.0, S_isov, 1.0);

  // scalar prefactors for computation of elasticity tensor
  double delta1 =
      (2. * G2 * modI4) / (std::pow(H2, 2.)) - (2. * B2 * G2 * H1 * modI4) / (std::pow(H2, 3.));
  double delta5 = (2. * B2 * G2) / (H2 * std::sqrt(modI4)) +
                  (4. * G2 * H3 * std::sqrt(modI4)) / (std::pow(H2, 2.)) -
                  (4. * B2 * G2 * H1 * H3 * std::sqrt(modI4)) / (std::pow(H2, 3.));
  double delta7 =
      (24. * G1 * modI5) / (std::pow(modI4, 4.)) +
      (8. * G2 * std::pow(H3, 2.)) / (std::pow(H2, 2.)) +
      (2. * B2 * G2 * (3. * H1 - 2. * modI1 * std::sqrt(modI4))) / (H2 * std::pow(modI4, 2.)) -
      (8. * B2 * G2 * std::pow(H3, 2.) * H1) / (std::pow(H2, 3.));
  delta7 += 2. * (-sigma_fiber_total / std::pow(modI4, 2.) + deriv_sigma_fiber_total / modI4);
  double delta8 = (2. * B2 * G2 * H1) / (std::pow(H2, 3.)) - (2. * G2) / (std::pow(H2, 2.));
  double delta10 =
      (2. * G2) / (std::pow(H2, 2) * modI4) - (2. * B2 * G2 * H1) / (modI4 * std::pow(H2, 3.));
  double delta11 = (2. * B2 * G2) / (H2 * std::pow(modI4, 3. / 2.)) -
                   (8. * G1) / (std::pow(modI4, 3.)) -
                   (4. * G2 * H3) / (std::pow(H2, 2.) * std::sqrt(modI4)) +
                   (4. * B2 * G2 * H3 * H1) / (std::sqrt(modI4) * std::pow(H2, 3.));
  double delta12 = (4. * G1) / (std::pow(modI4, 2.)) - (4. * B2 * G2) / (H2 * std::sqrt(modI4));

  // summand tensors for computation of elasticity tensor
  LINALG::Matrix<6, 6> IdId(false);
  IdId.MultiplyNT(Id3v, Id3v);  // summand1 = dyad(Id,Id) in Voigt-Notation

  LINALG::Matrix<6, 6> IdMsumMId(false);
  IdMsumMId.MultiplyNT(Id3v, Mv);
  IdMsumMId.MultiplyNT(
      1.0, Mv, Id3v, 1.0);  // summand5 = dyad(Id,Md) + dyad(M,Id) in Voigt-Notation

  LINALG::Matrix<6, 6> MM(false);
  MM.MultiplyNT(Mv, Mv);  // summand7 = dyad(M,M) in Voigt-Notation

  LINALG::Matrix<6, 6> IddI5sumdI5Id(false);
  IddI5sumdI5Id.MultiplyNT(Id3v, modCMsumMmodCv);
  IddI5sumdI5Id.MultiplyNT(
      1.0, modCMsumMmodCv, Id3v, 1.0);  // summand8 = dyad(Id,dI5) + dyad(dI5,Id)

  LINALG::Matrix<6, 6> dI5dI5(false);
  dI5dI5.MultiplyNT(modCMsumMmodCv, modCMsumMmodCv);  // summand10 = dyad(dI5,dI5)

  LINALG::Matrix<6, 6> MdI5sumdI5M(false);
  MdI5sumdI5M.MultiplyNT(Mv, modCMsumMmodCv);
  MdI5sumdI5M.MultiplyNT(1.0, modCMsumMmodCv, Mv, 1.0);  // summand11 = dyad(M,dI5) + dyad(dI5,M)

  // ficticious elasticiy tensor
  LINALG::Matrix<6, 6> modcmat(false);
  modcmat.Update(delta1, IdId);
  modcmat.Update(delta5, IdMsumMId, 1.0);
  modcmat.Update(delta7, MM, 1.0);
  modcmat.Update(delta8, IddI5sumdI5Id, 1.0);
  modcmat.Update(delta10, dI5dI5, 1.0);
  modcmat.Update(delta11, MdI5sumdI5M, 1.0);
  ElastSymTensorMultiply(modcmat, delta12, Id3, M, 1.0);  // summand 12 = ddI5/dC^2 = Id_ik*M_jl ...
  ElastSymTensorMultiply(modcmat, delta12, M, Id3, 1.0);  // ... + M_ik*Id_jl
  modcmat.Scale(std::pow(J, -4.0 / 3.0));

  // modified projection tensor Psl = Cinv o Cinv - 1/3 Cinv x Cinv
  LINALG::Matrix<6, 6> Psl(false);
  Psl.Clear();
  AddtoCmatHolzapfelProduct(Psl, icg, 1.0);
  Psl.MultiplyNT(-1.0 / 3.0, icg, icg, 1.0);

  // Right Cauchy-Green tensor in stress-like Voigt notation
  static LINALG::Matrix<6, 1> rcg_stress(false);
  UTILS::VOIGT::Strains::ToStressLike(rcg, rcg_stress);

  // compute the projection tensor P = II - 1/3 Cinv x C
  LINALG::Matrix<6, 6> P(false);
  UTILS::VOIGT::FourthOrderIdentityMatrix<UTILS::VOIGT::NotationType::stress,
      UTILS::VOIGT::NotationType::stress>(P);
  P.MultiplyNT(-1.0 / 3.0, icg, rcg_stress, 1.0);

  // compute the transpose of the projection tensor PT = II - 1/3 C x Cinv
  LINALG::Matrix<6, 6> PT(false);
  UTILS::VOIGT::FourthOrderIdentityMatrix<UTILS::VOIGT::NotationType::stress,
      UTILS::VOIGT::NotationType::stress>(PT);
  PT.MultiplyNT(-1.0 / 3.0, rcg_stress, icg, 1.0);

  // compute isochoric cmat from ficticious elasticiy tensor
  LINALG::Matrix<6, 6> cmatiso(false);
  cmatiso.MultiplyNN(P, modcmat);
  cmatiso.MultiplyNN(1.0, modcmat, PT, 1.0);
  cmatiso.MultiplyNT(-2.0 / 3.0, icg, S_isov, 1.0);
  cmatiso.Update(2.0 / 3.0 * incJ * traceCmodS, Psl, 1.0);
  cmatiso.MultiplyNT(-2.0 / 3.0, icg, S_isov, 1.0);
  cmatiso.MultiplyNT(-2.0 / 3.0, S_isov, icg, 1.0);

  // update cmat
  cmat.Update(1.0, cmatiso, 1.0);
}
double MAT::ELASTIC::IsoMuscleBlemker::EvaluateTimeDependentActivationLevel(
    double alpha, double beta, double t_tot)
{
  // get start of activation
  const double& t_act_start = params_->t_act_start_;

  // compute activation level
  double alpha_act = 0;
  if (t_tot >= t_act_start)
  {
    alpha_act = alpha * std::tanh(beta * (t_tot - t_act_start));
  }
  return alpha_act;
}

void MAT::ELASTIC::IsoMuscleBlemker::EvaluateTotalFiberCauchyStressAndDerivative(
    double modI4, double alpha_act, double& sigma_fiber_total, double& deriv_sigma_fiber_total)
{
  // get active material parameters
  const double& sigma_max = params_->sigma_max_;
  const double& lambda_star = params_->lambda_star_;
  const double& lambda_ofl = params_->lambda_ofl_;
  const double& P1 = params_->P1_;
  const double& P2 = params_->P2_;

  // compute fiber stretch and stretch ratio
  double lambda = std::sqrt(modI4);
  double ratio_lambda = lambda / lambda_ofl;

  // initialize fpassive and factive and the respective derivatives w.r.t. the fibre stretch
  double f_passive = 0;
  double deriv_f_passive = 0;
  double f_active = 0;
  double deriv_f_active = 0;

  // calculate constants for computation of passive fiber force
  double P3 = P1 * P2 * std::exp(P2 * (lambda_star / lambda_ofl - 1.0)) / lambda_ofl;
  double P4 = P1 * (std::exp(P2 * ((lambda_star / lambda_ofl) - 1.0)) - 1.0) -
              P3 * lambda_star / lambda_ofl;

  // compute normalized passive fiber force and derivative w.r.t. fibre stretch
  if (lambda <= lambda_ofl)
  {
    f_passive = 0.0;
    deriv_f_passive = 0.0;
  }
  else if (lambda < lambda_star)
  {
    f_passive = P1 * (std::exp(P2 * (ratio_lambda - 1.0)) - 1.0);
    deriv_f_passive = P1 * std::exp(P2 * (ratio_lambda - 1.0)) * P2 / lambda_ofl;
  }
  else
  {
    f_passive = P3 * ratio_lambda + P4;
    deriv_f_passive = P3 / lambda_ofl;
  }

  // compute normalized normalized active fiber force and derivative w.r.t. fibre stretch
  if (lambda <= 0.6 * lambda_ofl)
  {
    f_active = 9 * std::pow(ratio_lambda - 0.4, 2.0);
    deriv_f_active = 18 / lambda_ofl * (ratio_lambda - 0.4);
  }
  else if (lambda < 1.4 * lambda_ofl)
  {
    f_active = 1 - 4 * std::pow(1 - ratio_lambda, 2.0);
    deriv_f_active = 8 / lambda_ofl * (1 - ratio_lambda);
  }
  else
  {
    f_active = 9 * std::pow(ratio_lambda - 1.6, 2.0);
    deriv_f_active = 18 / lambda_ofl * (ratio_lambda - 1.6);
  }

  sigma_fiber_total = sigma_max * ratio_lambda * (f_passive + alpha_act * f_active);
  deriv_sigma_fiber_total =
      sigma_max * ratio_lambda * (deriv_f_passive + alpha_act * deriv_f_active) +
      sigma_max / lambda_ofl * sigma_fiber_total;
}