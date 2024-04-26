/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a hyperelastic constituent with a damage process and a 2D membrane material

\level 3

*/
/*----------------------------------------------------------------------*/

#include "4C_mixture_constituent_elasthyper_elastin_membrane.hpp"

#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix_generators.hpp"
#include "4C_mat_anisotropy_extension.hpp"
#include "4C_mat_multiplicative_split_defgrad_elasthyper_service.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_matelast_aniso_structuraltensor_strategy.hpp"
#include "4C_matelast_isoneohooke.hpp"
#include "4C_mixture_elastin_membrane_prestress_strategy.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN

MIXTURE::ElastinMembraneAnisotropyExtension::ElastinMembraneAnisotropyExtension(
    const Teuchos::RCP<MAT::ELASTIC::StructuralTensorStrategyBase>& structuralTensorStrategy)
    : FiberAnisotropyExtension<1>(structuralTensorStrategy)
{
  RegisterNeededTensors(
      FiberAnisotropyExtension::FIBER_VECTORS | FiberAnisotropyExtension::STRUCTURAL_TENSOR);
}

void MIXTURE::ElastinMembraneAnisotropyExtension::OnGlobalDataInitialized()
{
  if (GetAnisotropy()->HasGPCylinderCoordinateSystem())
  {
    for (int gp = 0; gp < GetAnisotropy()->GetNumberOfGaussPoints(); ++gp)
    {
      std::array<CORE::LINALG::Matrix<3, 1>, 1> fibers;

      fibers[0].Update(GetAnisotropy()->GetGPCylinderCoordinateSystem(gp).GetRad());
      MAT::FiberAnisotropyExtension<1>::SetFibers(gp, fibers);
    }

    orthogonal_structural_tensor_.resize(GetAnisotropy()->GetNumberOfGaussPoints());

    SetFiberLocation(MAT::FiberLocation::GPFibers);
  }
  else if (GetAnisotropy()->HasElementCylinderCoordinateSystem())
  {
    std::array<CORE::LINALG::Matrix<3, 1>, 1> fibers;

    fibers[0].Update(GetAnisotropy()->GetElementCylinderCoordinateSystem().GetRad());
    FiberAnisotropyExtension<1>::SetFibers(BaseAnisotropyExtension::GPDEFAULT, fibers);

    orthogonal_structural_tensor_.resize(1);

    SetFiberLocation(MAT::FiberLocation::ElementFibers);
  }
  else
  {
    FOUR_C_THROW(
        "Could not find either GP cylinder coordinates nor element cylinder coordinates. The "
        "elastin material needs at least one.");
  }

  for (unsigned gp = 0; gp < orthogonal_structural_tensor_.size(); ++gp)
  {
    orthogonal_structural_tensor_[gp] = CORE::LINALG::IdentityMatrix<3>();

    orthogonal_structural_tensor_[gp].Update(-1.0, GetStructuralTensor(gp, 0), 1.0);
  }
}

const CORE::LINALG::Matrix<3, 3>&
MIXTURE::ElastinMembraneAnisotropyExtension::GetOrthogonalStructuralTensor(int gp)
{
  if (orthogonal_structural_tensor_.empty())
  {
    FOUR_C_THROW("The coordinate system hast not been initialized yet.");

    // avoid compiler error
    return orthogonal_structural_tensor_[0];
  }

  if (orthogonal_structural_tensor_.size() == 1)
  {
    // using element fibers
    return orthogonal_structural_tensor_[0];
  }


  // using Gauss-point fibers
  return orthogonal_structural_tensor_[gp];
}

// Constructor for the parameter class
MIXTURE::PAR::MixtureConstituentElastHyperElastinMembrane::
    MixtureConstituentElastHyperElastinMembrane(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MixtureConstituentElastHyperBase(matdata),
      damage_function_id_(matdata->Get<int>("DAMAGE_FUNCT")),
      nummat_membrane_(matdata->Get<int>("MEMBRANENUMMAT")),
      matids_membrane_(matdata->Get<std::vector<int>>("MEMBRANEMATIDS"))
{
  if (nummat_membrane_ != (int)matids_membrane_.size())
  {
    FOUR_C_THROW(
        "number of membrane summands %d does not fit to the size of the membrane summands vector"
        " %d",
        nummat_membrane_, matids_membrane_.size());
  }
}

// Create an instance of MIXTURE::MixtureConstituentElastHyper from the parameters
std::unique_ptr<MIXTURE::MixtureConstituent>
MIXTURE::PAR::MixtureConstituentElastHyperElastinMembrane::CreateConstituent(int id)
{
  return std::unique_ptr<MIXTURE::MixtureConstituentElastHyperElastinMembrane>(
      new MIXTURE::MixtureConstituentElastHyperElastinMembrane(this, id));
}

// Constructor of the constituent holding the material parameters
MIXTURE::MixtureConstituentElastHyperElastinMembrane::MixtureConstituentElastHyperElastinMembrane(
    MIXTURE::PAR::MixtureConstituentElastHyperElastinMembrane* params, int id)
    : MixtureConstituentElastHyperBase(params, id),
      params_(params),
      anisotropy_extension_(Teuchos::rcp<MAT::ELASTIC::StructuralTensorStrategyBase>(
          new MAT::ELASTIC::StructuralTensorStrategyStandard(nullptr)))
{
  // Create summands
  for (const auto& matid : params_->matids_membrane_)
  {
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
    if (sum == Teuchos::null) FOUR_C_THROW("Failed to read elastic summand.");

    Teuchos::RCP<MAT::ELASTIC::IsoNeoHooke> neoHooke =
        Teuchos::rcp_dynamic_cast<MAT::ELASTIC::IsoNeoHooke>(sum);

    if (Teuchos::is_null(neoHooke))
    {
      FOUR_C_THROW(
          "Currently, only the an IsoNeoHooke material law is possible for use as an elastin "
          "membrane material");
    }

    potsum_membrane_.push_back(neoHooke);
  }
}

// Returns the material type
CORE::Materials::MaterialType MIXTURE::MixtureConstituentElastHyperElastinMembrane::MaterialType()
    const
{
  return CORE::Materials::mix_elasthyper_elastin_membrane;
}

// Pack the constituent
void MIXTURE::MixtureConstituentElastHyperElastinMembrane::PackConstituent(
    CORE::COMM::PackBuffer& data) const
{
  MixtureConstituentElastHyperBase::PackConstituent(data);

  CORE::COMM::ParObject::AddtoPack(data, current_reference_growth_);

  CORE::COMM::ParObject::AddtoPack(data, mue_frac_);

  anisotropy_extension_.PackAnisotropy(data);

  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (const auto& p : potsum_membrane_) p->PackSummand(data);
  }
}

// Unpack the constituent
void MIXTURE::MixtureConstituentElastHyperElastinMembrane::UnpackConstituent(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MixtureConstituentElastHyperBase::UnpackConstituent(position, data);

  CORE::COMM::ParObject::ExtractfromPack(position, data, current_reference_growth_);

  CORE::COMM::ParObject::ExtractfromPack(position, data, mue_frac_);

  anisotropy_extension_.UnpackAnisotropy(data, position);

  // loop map of associated potential summands
  for (auto& summand : potsum_membrane_) summand->UnpackSummand(data, position);
}

void MIXTURE::MixtureConstituentElastHyperElastinMembrane::RegisterAnisotropyExtensions(
    MAT::Anisotropy& anisotropy)
{
  MixtureConstituentElastHyperBase::RegisterAnisotropyExtensions(anisotropy);

  anisotropy.RegisterAnisotropyExtension(anisotropy_extension_);

  for (const auto& summand : potsum_membrane_) summand->RegisterAnisotropyExtensions(anisotropy);
}

// Reads the element from the input file
void MIXTURE::MixtureConstituentElastHyperElastinMembrane::ReadElement(
    int numgp, INPUT::LineDefinition* linedef)
{
  MixtureConstituentElastHyperBase::ReadElement(numgp, linedef);

  // Setup summands
  for (const auto& summand : potsum_membrane_) summand->Setup(numgp, linedef);

  current_reference_growth_.resize(numgp, 1.0);
  mue_frac_.resize(numgp, 1.0);
}

// Updates all summands
void MIXTURE::MixtureConstituentElastHyperElastinMembrane::Update(
    CORE::LINALG::Matrix<3, 3> const& defgrd, Teuchos::ParameterList& params, const int gp,
    const int eleGID)
{
  const auto& reference_coordinates = params.get<CORE::LINALG::Matrix<3, 1>>("gp_coords_ref");

  double totaltime = params.get<double>("total time", -1);
  if (totaltime < 0.0)
  {
    FOUR_C_THROW("Parameter 'total time' could not be read!");
  }

  current_reference_growth_[gp] =
      GLOBAL::Problem::Instance()
          ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(params_->damage_function_id_ - 1)
          .Evaluate(reference_coordinates.A(), totaltime, 0);

  MixtureConstituentElastHyperBase::Update(defgrd, params, gp, eleGID);

  // loop map of associated potential summands
  for (auto& summand : potsum_membrane_) summand->Update();
}

void MIXTURE::MixtureConstituentElastHyperElastinMembrane::PreEvaluate(
    MixtureRule& mixtureRule, Teuchos::ParameterList& params, int gp, int eleGID)
{
  MIXTURE::MixtureConstituentElastHyperBase::PreEvaluate(mixtureRule, params, gp, eleGID);

  // Evaluate mue frac
  std::shared_ptr<ElastinMembranePrestressStrategy> strategy =
      std::dynamic_pointer_cast<ElastinMembranePrestressStrategy>(PrestressStrategy());

  if (strategy == nullptr)
  {
    FOUR_C_THROW(
        "The used prestretch strategy is not compatible with elastin materials. It has to "
        "implement MIXTURE::ElastinMembranePrestressStrategy.");
  }

  mue_frac_[gp] = strategy->EvaluateMueFrac(mixtureRule,
      CylinderCoordinateSystemAnisotropyExtension().GetCoordinateSystemProvider(gp), *this, *this,
      params, gp, eleGID);
}

double MIXTURE::MixtureConstituentElastHyperElastinMembrane::GetGrowthScalar(int gp) const
{
  return current_reference_growth_[gp];
}

void MIXTURE::MixtureConstituentElastHyperElastinMembrane::Evaluate(
    const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<6, 1>& E_strain,
    Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
    CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  FOUR_C_THROW("This constituent does not support Evaluation without an elastic part.");
}

void MIXTURE::MixtureConstituentElastHyperElastinMembrane::EvaluateElasticPart(
    const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<3, 3>& iFextin,
    Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
    CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  // Compute total inelastic deformation gradient
  static CORE::LINALG::Matrix<3, 3> iFin(false);
  iFin.MultiplyNN(iFextin, PrestretchTensor(gp));

  // Evaluate 3D elastic part
  MAT::ElastHyperEvaluateElasticPart(
      F, iFin, S_stress, cmat, Summands(), SummandProperties(), gp, eleGID);

  // Evaluate Membrane
  static CORE::LINALG::Matrix<6, 1> Smembrane_stress(false);
  static CORE::LINALG::Matrix<6, 6> cmatmembrane(false);
  EvaluateStressCMatMembrane(F, iFin, params, Smembrane_stress, cmatmembrane, gp, eleGID);

  S_stress.Update(1.0, Smembrane_stress, 1.0);
  cmat.Update(1.0, cmatmembrane, 1.0);
}

void MIXTURE::MixtureConstituentElastHyperElastinMembrane::EvaluateMembraneStress(
    CORE::LINALG::Matrix<6, 1>& S, Teuchos::ParameterList& params, int gp, int eleGID)
{
  CORE::LINALG::Matrix<6, 6> cmat(false);
  const CORE::LINALG::Matrix<3, 3> Id = CORE::LINALG::IdentityMatrix<3>();
  CORE::LINALG::Matrix<3, 3> iFin(false);

  iFin.MultiplyNN(Id, PrestretchTensor(gp));

  EvaluateStressCMatMembrane(Id, iFin, params, S, cmat, gp, eleGID);
}

void MIXTURE::MixtureConstituentElastHyperElastinMembrane::EvaluateStressCMatMembrane(
    const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<3, 3>& iFin,
    Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
    CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID) const
{
  static CORE::LINALG::Matrix<3, 3> Ce(false);
  MAT::EvaluateCe(F, iFin, Ce);

  // Compute structural tensors in grown configuration
  static CORE::LINALG::Matrix<3, 3> Aradgr(false);
  static CORE::LINALG::Matrix<3, 3> Aorthgr(false);
  EvaluateStructuralTensorsInGrownConfiguration(Aradgr, Aorthgr, iFin, gp, eleGID);

  static CORE::LINALG::Matrix<3, 3> AorthgrCeAorthgrArad(false);
  EvaluateAorthgrCeAorthgrArad(AorthgrCeAorthgrArad, Aradgr, Aorthgr, Ce);
  double detX = AorthgrCeAorthgrArad.Determinant();

  static CORE::LINALG::Matrix<3, 3> iFinAorthgriFinT(false);
  EvaluateiFinAorthgriFinT(iFinAorthgriFinT, iFin, Aorthgr);

  // Z = F_{in}^{-T}*A_{gr}^{T}*X^{-1}*A_{gr}^{T}*F_{in}^{-1}
  static CORE::LINALG::Matrix<3, 3> iFinTAorthgrTiXTAorthgriFin_sym(false);
  static CORE::LINALG::Matrix<6, 1> iFinTAorthgrTiXTAorthgriFin_sym_stress(false);

  EvaluateiFinTAorthgrTiXTAorthgriFin(
      iFinTAorthgrTiXTAorthgriFin_sym, AorthgrCeAorthgrArad, iFin, Aorthgr);

  CORE::LINALG::VOIGT::Stresses::MatrixToVector(
      iFinTAorthgrTiXTAorthgriFin_sym, iFinTAorthgrTiXTAorthgriFin_sym_stress);

  // Get material parameter
  double mue = 0.0;
  for (const auto& summand : potsum_membrane_)
  {
    mue += summand->Mue();
  }

  // Compute membrane stress
  static CORE::LINALG::Matrix<3, 3> Smembrane;
  Smembrane.Update(mue * mue_frac_[gp], iFinAorthgriFinT, 0.0);
  Smembrane.Update(-mue * mue_frac_[gp] / detX, iFinTAorthgrTiXTAorthgriFin_sym, 1.0);

  CORE::LINALG::VOIGT::Stresses::MatrixToVector(Smembrane, S_stress);

  // Compute constitutive tensor
  static CORE::LINALG::Matrix<6, 6> dAradgriXAradgr_symdC(false);
  dAradgriXAradgr_symdC.Clear();

  MAT::AddtoCmatHolzapfelProduct(
      dAradgriXAradgr_symdC, iFinTAorthgrTiXTAorthgriFin_sym_stress, -2.0);

  cmat.MultiplyNT(2.0 * mue * mue_frac_[gp] / detX, iFinTAorthgrTiXTAorthgriFin_sym_stress,
      iFinTAorthgrTiXTAorthgriFin_sym_stress, 0.0);
  cmat.Update(-mue * mue_frac_[gp] / detX, dAradgriXAradgr_symdC, 1.0);
}

void MIXTURE::MixtureConstituentElastHyperElastinMembrane::
    EvaluateStructuralTensorsInGrownConfiguration(CORE::LINALG::Matrix<3, 3>& Aradgr,
        CORE::LINALG::Matrix<3, 3>& Aorthgr, const CORE::LINALG::Matrix<3, 3>& iFin, const int gp,
        const int eleGID) const
{
  // Compute inelastic right Cauchy-Green deformation gradient
  static CORE::LINALG::Matrix<3, 3> iCin(false);
  iCin.MultiplyNT(iFin, iFin);

  static CORE::LINALG::Matrix<3, 3> Fin(false);
  Fin.Invert(iFin);

  // Compute radial structural tensor in grown configuration
  static CORE::LINALG::Matrix<3, 3> FinArad(false);
  FinArad.MultiplyNN(Fin, anisotropy_extension_.GetStructuralTensor(gp, 0));
  Aradgr.MultiplyNT(iCin.Dot(anisotropy_extension_.GetStructuralTensor(gp, 0)), FinArad, Fin, 0.0);

  // Compute orthogonal (to radial) structural tensor in grown configuration
  Aorthgr = CORE::LINALG::IdentityMatrix<3>();
  Aorthgr.Update(-1.0, Aradgr, 1.0);
}

void MIXTURE::MixtureConstituentElastHyperElastinMembrane::EvaluateAorthgrCeAorthgrArad(
    CORE::LINALG::Matrix<3, 3>& AorthgrCeAorthgrArad, const CORE::LINALG::Matrix<3, 3>& Aradgr,
    const CORE::LINALG::Matrix<3, 3>& Aorthgr, const CORE::LINALG::Matrix<3, 3>& Ce)
{
  static CORE::LINALG::Matrix<3, 3> AorthgrCe(false);
  AorthgrCe.MultiplyNN(Aorthgr, Ce);
  AorthgrCeAorthgrArad.MultiplyNN(AorthgrCe, Aorthgr);
  AorthgrCeAorthgrArad.Update(1.0, Aradgr, 1.0);
}

void MIXTURE::MixtureConstituentElastHyperElastinMembrane::EvaluateiFinAorthgriFinT(
    CORE::LINALG::Matrix<3, 3>& iFinAorthgriFinT, const CORE::LINALG::Matrix<3, 3>& iFin,
    const CORE::LINALG::Matrix<3, 3>& Aorthgr)
{
  static CORE::LINALG::Matrix<3, 3> iFinAorthgr(false);
  iFinAorthgr.MultiplyNN(iFin, Aorthgr);
  iFinAorthgriFinT.MultiplyNT(iFinAorthgr, iFin);
}

void MIXTURE::MixtureConstituentElastHyperElastinMembrane::EvaluateiFinTAorthgrTiXTAorthgriFin(
    CORE::LINALG::Matrix<3, 3>& iFinTAorthgrTiXTAorthgriFin_sym,
    const CORE::LINALG::Matrix<3, 3>& AorthgrCeAorthgrArad, const CORE::LINALG::Matrix<3, 3>& iFin,
    const CORE::LINALG::Matrix<3, 3>& Aorthgr)
{
  static CORE::LINALG::Matrix<3, 3> iAorthgrCeAorthgrArad(false);
  iAorthgrCeAorthgrArad.Invert(AorthgrCeAorthgrArad);

  static CORE::LINALG::Matrix<3, 3> AorthgriFin(false);
  AorthgriFin.MultiplyNN(Aorthgr, iFin);

  static CORE::LINALG::Matrix<3, 3> iFinTAorthgrTiXT(false);
  iFinTAorthgrTiXT.MultiplyTT(AorthgriFin, iAorthgrCeAorthgrArad);

  static CORE::LINALG::Matrix<3, 3> iFinTAorthgrTiXTAorthgriFin(false);
  iFinTAorthgrTiXTAorthgriFin.MultiplyNN(iFinTAorthgrTiXT, AorthgriFin);

  iFinTAorthgrTiXTAorthgriFin_sym.Update(0.5, iFinTAorthgrTiXTAorthgriFin, 0.0);
  iFinTAorthgrTiXTAorthgriFin_sym.UpdateT(0.5, iFinTAorthgrTiXTAorthgriFin, 1.0);
}
FOUR_C_NAMESPACE_CLOSE
