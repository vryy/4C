/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of a general remodel fiber constituent (basis for an explicit or implicit update
rule)

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/
#include "mixture_constituent_remodelfiber.H"
#include "../drt_lib/voigt_notation.H"
#include "../drt_mat/material_service.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_matelast/elast_coupanisoexpo.H"
#include "../drt_matelast/elast_coupanisoexpoactive.H"
#include "../drt_mat/anisotropy_extension.H"
#include "growth_evolution_linear_cauchy.H"
#include "../drt_matelast/elast_summand.H"
#include "../drt_matelast/elast_activesummand.H"


MIXTURE::PAR::MixtureConstituent_RemodelFiber::MixtureConstituent_RemodelFiber(
    const Teuchos::RCP<MAT::PAR::Material>& matdata, double ref_mass_fraction)
    : MixtureConstituent(matdata, ref_mass_fraction),
      matid_(matdata->GetInt("MATID")),
      poisson_decay_time_(matdata->GetDouble("DECAY_TIME")),
      growth_constant_(matdata->GetDouble("GROWTH_CONSTANT")),
      deposition_stretch_(matdata->GetDouble("DEPOSITION_STRETCH"))
{
}

MIXTURE::MixtureConstituent_RemodelFiber::MixtureConstituent_RemodelFiber(
    MIXTURE::PAR::MixtureConstituent_RemodelFiber* params)
    : MixtureConstituent(params), params_(params), orthogonalAnisotropyExtension_()
{
  summand_ = MAT::ELASTIC::Summand::Factory(params_->matid_);
  activeSummand_ = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::ActiveSummand>(summand_);

  if (summand_ == Teuchos::null)
  {
    dserror(
        "Failed to read strain energy summand of the remodel fiber. Use materials from the "
        "drt_matelast library.");
  }

  auto coupAnisoExpo = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpo>(summand_);
  auto coupAnisoExpoActive = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::CoupAnisoExpoActive>(summand_);
  // safety check: Currently only MAT::ELASTIC::CoupAnisoExpo is supported.
  if (coupAnisoExpo == Teuchos::null and coupAnisoExpoActive == Teuchos::null)
  {
    dserror(
        "Currently only the MAT::ELASTIC::CoupAnisoExpo and MAT::ELASTIC::CoupAnisoExpoActive "
        "summand are supported. You need to implement the derivatives with respect to the "
        "anisotropic invariants and extend the remodel fiber to handle materials that depent on "
        "the second anisotropic invariant.");
  }

  // Get anisotropy extension from CoupAnisoExpo
  auto anisotropyExtensionProvider =
      Teuchos::rcp_dynamic_cast<MAT::FiberAnisotropyExtensionProvider<1>>(summand_);
  if (coupAnisoExpo == Teuchos::null and coupAnisoExpoActive == Teuchos::null)
  {
    dserror(
        "The summand has to implement AnisotropyExtensionProvider so that the remodel fiber has "
        "access to the anisotropy information of the summand.");
  }

  fiberAnisotropyExtension_ =
      Teuchos::rcpFromRef(anisotropyExtensionProvider->GetFiberAnisotropyExtension());

  fiberAnisotropyExtension_->RegisterNeededTensors(
      MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR |
      MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR_STRESS);
  orthogonalAnisotropyExtension_.SetAnisotropyExtension(fiberAnisotropyExtension_);

  auto growthEvolution = Teuchos::rcp<CauchyLinearGrowthEvolution>(
      new CauchyLinearGrowthEvolution(params->growth_constant_));

  growth_evolution_ = Teuchos::rcp_static_cast<GrowthEvolution>(growthEvolution);
}

void MIXTURE::MixtureConstituent_RemodelFiber::PackConstituent(DRT::PackBuffer& data) const
{
  MixtureConstituent::PackConstituent(data);

  DRT::ParObject::AddtoPack(data, iFr_);

  DRT::ParObject::AddtoPack(data, cur_rho_);
  DRT::ParObject::AddtoPack(data, cur_lambdar_);
  DRT::ParObject::AddtoPack(data, lambda_pre_);

  DRT::ParObject::AddtoPack(data, sig_h_);

  summand_->PackSummand(data);
}

void MIXTURE::MixtureConstituent_RemodelFiber::UnpackConstituent(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MixtureConstituent::UnpackConstituent(position, data);

  DRT::ParObject::ExtractfromPack(position, data, iFr_);
  DRT::ParObject::ExtractfromPack(position, data, cur_rho_);
  DRT::ParObject::ExtractfromPack(position, data, cur_lambdar_);
  DRT::ParObject::ExtractfromPack(position, data, lambda_pre_);

  DRT::ParObject::ExtractfromPack(position, data, sig_h_);

  summand_->UnpackSummand(data, position);
}

void MIXTURE::MixtureConstituent_RemodelFiber::RegisterAnisotropyExtensions(
    MAT::Anisotropy& anisotropy)
{
  summand_->RegisterAnisotropyExtensions(anisotropy);
  anisotropy.RegisterAnisotropyExtension(orthogonalAnisotropyExtension_);
}

void MIXTURE::MixtureConstituent_RemodelFiber::ReadElement(
    int numgp, DRT::INPUT::LineDefinition* linedef)
{
  MixtureConstituent::ReadElement(numgp, linedef);

  summand_->Setup(numgp, linedef);

  // From here we know the number of Gauss-points
  // Resize all GP-vectors
  iFr_.resize(numgp);
  cur_rho_.resize(numgp, params_->RefMassFraction() * InitialRefDensity());
  cur_lambdar_.resize(numgp, 1.0);
  lambda_pre_.resize(numgp, params_->deposition_stretch_);
  sig_h_.resize(numgp, 0.0);
}

void MIXTURE::MixtureConstituent_RemodelFiber::EvaluateElasticPart(const LINALG::Matrix<3, 3>& FM,
    const LINALG::Matrix<3, 3>& iFextin, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  // read dt from ParameterList
  // ToDo: This is on GP level -> Find sth that is more efficient than reading from the
  //       ParameterList
  double dt = params.get<double>("delta time", -1.0);
  if (dt < 0)
  {
    dserror("Looks like the material does not write dt into the Parameter list");
  }

  // Newton Update of internal variables
  UpdateNewton(dt, gp);

  // compute inverse of the inelastic deformation gradient
  static LINALG::Matrix<3, 3> iFin(false);
  iFin.MultiplyNN(iFextin, iFr_[gp]);

  // Compute and add stresses and linearization
  AddStressCmat(FM, iFin, S_stress, cmat, params, gp, eleGID);
}

void MIXTURE::MixtureConstituent_RemodelFiber::Evaluate(const LINALG::Matrix<3, 3>& F,
    const LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  dserror(
      "The growth remodel fiber can currently only be used with an inelastic deformation. However, "
      "the extension is easily possible.");
}

double MIXTURE::MixtureConstituent_RemodelFiber::CurrentRefDensity(int gp) const
{
  return cur_rho_[gp];
}

void MIXTURE::MixtureConstituent_RemodelFiber::UpdateNewton(double dt, int gp)
{
  iFr_[gp].Update(lambda_pre_[gp] / cur_lambdar_[gp],
      fiberAnisotropyExtension_->GetStructuralTensor(gp, 0), 0.0);
  iFr_[gp].Update(1. / std::sqrt(lambda_pre_[gp] / cur_lambdar_[gp]),
      orthogonalAnisotropyExtension_.OrthogonalStructuralTensor(gp), 1.0);
}

void MIXTURE::MixtureConstituent_RemodelFiber::Setup(
    Teuchos::ParameterList& params, const int eleGID)
{
  MixtureConstituent::Setup(params, eleGID);

  UpdateSigH(eleGID);
}

void MIXTURE::MixtureConstituent_RemodelFiber::Update(LINALG::Matrix<3, 3> const& defgrd,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  MixtureConstituent::Update(defgrd, params, gp, eleGID);
}

void MIXTURE::MixtureConstituent_RemodelFiber::AddStressCmat(const LINALG::Matrix<3, 3>& F,
    const LINALG::Matrix<3, 3>& iFin, LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  static LINALG::Matrix<3, 3> C(false);
  C.MultiplyTN(F, F);

  // Compute elastic right Cauchy-Green deformation tensor
  static LINALG::Matrix<3, 3> Ce(false);
  EvaluateCe(C, iFin, Ce);

  // Evaluate the first/second derivative of the free-energy function w.r.t. C
  static LINALG::Matrix<3, 3> dPsidC(true);
  static LINALG::Matrix<6, 1> dPsidC_stress(true);
  static LINALG::Matrix<6, 6> ddPsidCdC_stress(true);

  EvaluatedPsidC(Ce, iFin, dPsidC, gp, eleGID);
  EvaluateddPsidCdC(Ce, iFin, ddPsidCdC_stress, gp, eleGID);

  // Add stress contributions
  UTILS::VOIGT::Stresses::MatrixToVector(dPsidC, dPsidC_stress);
  S_stress.Update(2.0 * CurrentRefDensity(gp), dPsidC_stress, 1.0);
  cmat.Update(4.0 * CurrentRefDensity(gp), ddPsidCdC_stress, 1.0);

  if (!Teuchos::is_null(activeSummand_))
  {
    static LINALG::Matrix<6, 1> activeStress(false);
    static LINALG::Matrix<6, 6> activeCMat(false);
    activeStress.Clear();
    activeCMat.Clear();

    activeSummand_->AddActiveStressCmatAniso(
        C, activeCMat, activeStress, gp, eleGID);  // FIXME: This should be Ce instead of C ???

    S_stress.Update(CurrentRefDensity(gp), activeStress, 1.0);
    cmat.Update(CurrentRefDensity(gp), activeCMat, 1.0);
  }
}

void MIXTURE::MixtureConstituent_RemodelFiber::UpdateGrowthAndRemodelingExpl(
    const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<3, 3>& iFg, double dt, int gp, int eleGID)
{
  // Update internal variables
  UpdateNewton(dt, gp);

  // Compute right Cauchy-Green tensor
  static LINALG::Matrix<3, 3> C(false);
  C.MultiplyTN(F, F);

  // Compute inverse inelastic deformation gradient
  static LINALG::Matrix<3, 3> iFin(false);
  iFin.MultiplyNN(iFg, iFr_[gp]);

  // Evaluate Cauchy fiber normal stress
  const double sig = EvaluateFiberCauchyStress(C, iFin, gp, eleGID);

  // Evaluate growth ode
  double delta_sig = (sig - sig_h_[gp]);
  const double drhodt = EvaluateGrowthEvolutionEquationdt(delta_sig, gp, eleGID);

  // Evaluate remodel ode
  const double dlambdardt = EvaluateRemodelEvolutionEquationdt(C, iFin, delta_sig, gp, eleGID);

  // Apply explicit update step
  cur_rho_[gp] += drhodt * dt;          // growth
  cur_lambdar_[gp] += dlambdardt * dt;  // remodel
}

double MIXTURE::MixtureConstituent_RemodelFiber::EvaluateFiberCauchyStress(
    const LINALG::Matrix<3, 3>& C, const LINALG::Matrix<3, 3>& iFin, int gp, int eleGID) const
{
  // Compute elastic right Cauchy-Green strain tensor
  static LINALG::Matrix<3, 3> Ce(false);
  EvaluateCe(C, iFin, Ce);


  static LINALG::Matrix<2, 1> dPIe(false);
  double dPIact = 0.0;

  summand_->EvaluateFirstDerivativesAniso(dPIe, Ce, gp, eleGID);

  if (!Teuchos::is_null(activeSummand_))
  {
    dPIact = activeSummand_->GetDerivativeAnisoActive();
  }

  return 2.0 * dPIe(0) * Ce.Dot(fiberAnisotropyExtension_->GetStructuralTensor(gp, 0)) + dPIact;
}

double MIXTURE::MixtureConstituent_RemodelFiber::EvaluateGrowthEvolutionEquationdt(
    double delta_sig, int gp, int eleGID) const
{
  // Evaluate right hand side of growth ode, which is
  return growth_evolution_->EvaluateRHS(delta_sig / sig_h_[gp]) * cur_rho_[gp];
}

double MIXTURE::MixtureConstituent_RemodelFiber::EvaluateRemodelEvolutionEquationdt(
    LINALG::Matrix<3, 3>& C, LINALG::Matrix<3, 3>& iFin, double delta_sig, int gp, int eleGID) const
{
  // Compute Ce = Cin^-T * C * C_in^-1
  static LINALG::Matrix<3, 3> Ce(false);
  EvaluateCe(C, iFin, Ce);

  static LINALG::Matrix<3, 3> dsigdCe(false);
  EvaluatedsigdCe(Ce, dsigdCe, gp, eleGID);

  // Compute time derivative of the growth deformation gradient, reduced with the time
  // derivative of the remodel fiber stretch dlambdardt
  static LINALG::Matrix<3, 3> Frdot_red;

  Frdot_red.Update(
      1.0 / lambda_pre_[gp], fiberAnisotropyExtension_->GetStructuralTensor(gp, 0), 0.0);
  Frdot_red.Update(-0.5 * std::pow(cur_lambdar_[gp], -3.0 / 2.0) * std::pow(lambda_pre_[gp], 0.5),
      orthogonalAnisotropyExtension_.OrthogonalStructuralTensor(gp), 1.0);

  // Compute velocity gradient, reduced with the time derivative of the remodel fiber stretch
  // dlambdardt
  static LINALG::Matrix<3, 3> Lr_red(false);
  Lr_red.MultiplyNN(Frdot_red, iFr_[gp]);

  static LINALG::Matrix<3, 3> CeLr_red(false);
  CeLr_red.MultiplyNN(Ce, Lr_red);

  return (growth_evolution_->EvaluateRHS(delta_sig / sig_h_[gp]) +
             1.0 / params_->poisson_decay_time_) *
         delta_sig / (2.0 * dsigdCe.Dot(CeLr_red));
}

void MIXTURE::MixtureConstituent_RemodelFiber::EvaluateCe(
    const LINALG::Matrix<3, 3>& C, const LINALG::Matrix<3, 3>& iFin, LINALG::Matrix<3, 3>& Ce)
{
  // Ce = iFin^T * C * iFin
  static LINALG::Matrix<3, 3> CiFin(false);
  CiFin.MultiplyNN(C, iFin);
  Ce.MultiplyTN(iFin, CiFin);
}

void MIXTURE::MixtureConstituent_RemodelFiber::EvaluatedPsidC(const LINALG::Matrix<3, 3>& Ce,
    const LINALG::Matrix<3, 3>& iFin, LINALG::Matrix<3, 3>& dPsidC, int gp, int eleGID) const
{
  // Evaluate derivatives of Free-Energy function
  static LINALG::Matrix<2, 1> dPI_aniso(false);
  summand_->EvaluateFirstDerivativesAniso(dPI_aniso, Ce, gp, eleGID);

  // Return derivative
  // The derivative is in general: 0.5*A_ij*( iFin_ik iFin_lj + iFin_il iFin_kj )
  // In case of a rotation free inelastic deformation gradient, this reduces to iFin * A * iFin
  static LINALG::Matrix<3, 3> iFinA(false);
  iFinA.MultiplyNN(iFin, fiberAnisotropyExtension_->GetStructuralTensor(gp, 0));
  dPsidC.MultiplyNN(dPI_aniso(0), iFinA, iFin);
}

void MIXTURE::MixtureConstituent_RemodelFiber::EvaluateddPsidCdC(const LINALG::Matrix<3, 3>& Ce,
    const LINALG::Matrix<3, 3>& iFin, LINALG::Matrix<6, 6>& ddPsidCdC_stress, int gp,
    int eleGID) const
{
  // Compute inverse inelastic right Cauchy-Green deformation tensor
  static LINALG::Matrix<3, 3> iCin(false);
  iCin.MultiplyNT(iFin, iFin);

  // Evaluate second derivative of the free-energy function with respect to the anisotropic
  // invariants
  static LINALG::Matrix<3, 1> ddPsiIIe_aniso(false);
  summand_->EvaluateSecondDerivativesAniso(ddPsiIIe_aniso, Ce, gp, eleGID);

  ddPsidCdC_stress.MultiplyNT(
      ddPsiIIe_aniso(0) *
          std::pow(iCin.Dot(fiberAnisotropyExtension_->GetStructuralTensor(gp, 0)), 2.0),
      fiberAnisotropyExtension_->GetStructuralTensor_stress(gp, 0),
      fiberAnisotropyExtension_->GetStructuralTensor_stress(gp, 0), 0.0);
}

void MIXTURE::MixtureConstituent_RemodelFiber::EvaluatedsigdCe(
    const LINALG::Matrix<3, 3>& Ce, LINALG::Matrix<3, 3>& dsigdCe, int gp, int eleGID) const
{
  // Compute dsigdCe (derivative of the Cauchy fiber stress w.r.t. the elastic RCG
  static LINALG::Matrix<2, 1> dPIe_aniso(false);
  static LINALG::Matrix<3, 1> ddPIIe_aniso(false);

  summand_->EvaluateFirstDerivativesAniso(dPIe_aniso, Ce, gp, eleGID);
  summand_->EvaluateSecondDerivativesAniso(ddPIIe_aniso, Ce, gp, eleGID);

  dsigdCe.Update(
      2.0 * (ddPIIe_aniso(0) * Ce.Dot(fiberAnisotropyExtension_->GetStructuralTensor(gp, 0)) +
                dPIe_aniso(0)),
      fiberAnisotropyExtension_->GetStructuralTensor(gp, 0));
}

void MIXTURE::MixtureConstituent_RemodelFiber::UpdateSigH(const int eleGID)
{
  LINALG::Matrix<3, 3> Id(false);
  MAT::IdentityMatrix(Id);

  LINALG::Matrix<3, 3> iFextin(false);

  for (unsigned gp = 0; gp < sig_h_.size(); ++gp)
  {
    iFextin.Update(lambda_pre_[gp], fiberAnisotropyExtension_->GetStructuralTensor(gp, 0));

    sig_h_[gp] = EvaluateFiberCauchyStress(Id, iFextin, gp, eleGID);
  }
}

MIXTURE::OrthogonalAnisotropyExtension::OrthogonalAnisotropyExtension()
    : FiberAnisotropyExtension(), other_(Teuchos::null), A_orth_(0)
{
}

void MIXTURE::OrthogonalAnisotropyExtension::SetAnisotropyExtension(
    Teuchos::RCP<MAT::FiberAnisotropyExtension<1>>& other)
{
  other_ = other;
  other_->RegisterNeededTensors(STRUCTURAL_TENSOR);
}

void MIXTURE::OrthogonalAnisotropyExtension::OnGlobalDataInitialized()
{
  FiberAnisotropyExtension::OnGlobalDataInitialized();

  LINALG::Matrix<3, 3> Id(false);
  MAT::IdentityMatrix(Id);

  A_orth_.resize(other_->GetFibersPerElement());
  for (int gp = 0; gp < other_->GetFibersPerElement(); ++gp)
  {
    A_orth_[gp].Update(1.0, Id, -1.0, other_->GetStructuralTensor(gp, 0));
  }
}
const LINALG::Matrix<3, 3>& MIXTURE::OrthogonalAnisotropyExtension::OrthogonalStructuralTensor(
    int gp) const
{
  if (gp > 0 and A_orth_.size() == 1) return A_orth_[0];
  return A_orth_[gp];
}

void MIXTURE::OrthogonalAnisotropyExtension::PackAnisotropy(DRT::PackBuffer& data) const
{
  FiberAnisotropyExtension::PackAnisotropy(data);
  DRT::ParObject::AddtoPack(data, A_orth_);
}

void MIXTURE::OrthogonalAnisotropyExtension::UnpackAnisotropy(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  FiberAnisotropyExtension::UnpackAnisotropy(data, position);
  DRT::ParObject::ExtractfromPack(position, data, A_orth_);
}