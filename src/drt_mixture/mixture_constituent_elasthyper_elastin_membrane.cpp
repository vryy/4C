/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a hyperelastic constituent with a damage process and a 2D membrane material

\level 3

*/
/*----------------------------------------------------------------------*/

#include "mixture_constituent_elasthyper_elastin_membrane.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/multiplicative_split_defgrad_elasthyper_service.H"
#include "../drt_matelast/elast_aniso_structuraltensor_strategy.H"
#include "../drt_matelast/elast_isoneohooke.H"
#include "elastin_membrane_prestress_strategy.H"
#include "../drt_mat/anisotropy_extension.H"

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
      std::array<LINALG::Matrix<3, 1>, 1> fibers;

      fibers[0].Update(GetAnisotropy()->GetGPCylinderCoordinateSystem(gp).GetRad());
      MAT::FiberAnisotropyExtension<1>::SetFibers(gp, fibers);
    }

    orthogonalStructuralTensor_.resize(GetAnisotropy()->GetNumberOfGaussPoints());

    SetFiberLocation(MAT::FiberLocation::GPFibers);
  }
  else if (GetAnisotropy()->HasElementCylinderCoordinateSystem())
  {
    std::array<LINALG::Matrix<3, 1>, 1> fibers;

    fibers[0].Update(GetAnisotropy()->GetElementCylinderCoordinateSystem().GetRad());
    FiberAnisotropyExtension<1>::SetFibers(BaseAnisotropyExtension::GPDEFAULT, fibers);

    orthogonalStructuralTensor_.resize(1);

    SetFiberLocation(MAT::FiberLocation::ElementFibers);
  }
  else
  {
    dserror(
        "Could not find either GP cylinder coordinates nor element cylinder coordinates. The "
        "elastin material needs at least one.");
  }

  for (unsigned gp = 0; gp < orthogonalStructuralTensor_.size(); ++gp)
  {
    MAT::IdentityMatrix(orthogonalStructuralTensor_[gp]);
    orthogonalStructuralTensor_[gp].Update(-1.0, GetStructuralTensor(gp, 0), 1.0);
  }
}

const LINALG::Matrix<3, 3>&
MIXTURE::ElastinMembraneAnisotropyExtension::GetOrthogonalStructuralTensor(int gp)
{
  if (orthogonalStructuralTensor_.empty())
  {
    dserror("The coordinate system hast not been initialized yet.");

    // avoid compiler error
    return orthogonalStructuralTensor_[0];
  }

  if (orthogonalStructuralTensor_.size() == 1)
  {
    // using element fibers
    return orthogonalStructuralTensor_[0];
  }


  // using Gauss-point fibers
  return orthogonalStructuralTensor_[gp];
}

// Constructor for the parameter class
MIXTURE::PAR::MixtureConstituent_ElastHyperElastinMembrane::
    MixtureConstituent_ElastHyperElastinMembrane(
        const Teuchos::RCP<MAT::PAR::Material>& matdata, const double ref_mass_fraction)
    : MixtureConstituent_ElastHyperBase(matdata, ref_mass_fraction),
      damage_function_id_(matdata->GetInt("DAMAGE_FUNCT")),
      nummat_membrane_(matdata->GetInt("MEMBRANENUMMAT")),
      matids_membrane_(matdata->Get<std::vector<int>>("MEMBRANEMATIDS"))
{
  if (nummat_membrane_ != (int)matids_membrane_->size())
  {
    dserror(
        "number of membrane summands %d does not fit to the size of the membrane summands vector"
        " %d",
        nummat_membrane_, matids_membrane_->size());
  }
}

// Create an instance of MIXTURE::MixtureConstituent_ElastHyper from the parameters
Teuchos::RCP<MIXTURE::MixtureConstituent>
MIXTURE::PAR::MixtureConstituent_ElastHyperElastinMembrane::CreateConstituent()
{
  return Teuchos::rcp(new MIXTURE::MixtureConstituent_ElastHyperElastinMembrane(this));
}

// Constructor of the constituent holding the material parameters
MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::MixtureConstituent_ElastHyperElastinMembrane(
    MIXTURE::PAR::MixtureConstituent_ElastHyperElastinMembrane* params)
    : MixtureConstituent_ElastHyperBase(params),
      params_(params),
      anisotropyExtension_(Teuchos::rcp<MAT::ELASTIC::StructuralTensorStrategyBase>(
          new MAT::ELASTIC::StructuralTensorStrategyStandard(nullptr)))
{
  // Create summands
  for (const auto& matid : *params_->matids_membrane_)
  {
    Teuchos::RCP<MAT::ELASTIC::Summand> sum = MAT::ELASTIC::Summand::Factory(matid);
    if (sum == Teuchos::null) dserror("Failed to read elastic summand.");

    Teuchos::RCP<MAT::ELASTIC::IsoNeoHooke> neoHooke =
        Teuchos::rcp_dynamic_cast<MAT::ELASTIC::IsoNeoHooke>(sum);

    if (Teuchos::is_null(neoHooke))
    {
      dserror(
          "Currently, only the an IsoNeoHooke material law is possible for use as an elastin "
          "membrane material");
    }

    potsum_membrane_.push_back(neoHooke);
  }
}

// Returns the material type
INPAR::MAT::MaterialType MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::MaterialType() const
{
  return INPAR::MAT::mix_elasthyper_elastin_membrane;
}

// Pack the constituent
void MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::PackConstituent(
    DRT::PackBuffer& data) const
{
  MixtureConstituent_ElastHyperBase::PackConstituent(data);

  DRT::ParObject::AddtoPack(data, current_reference_density_);

  DRT::ParObject::AddtoPack(data, mue_frac_);

  anisotropyExtension_.PackAnisotropy(data);

  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (const auto& p : potsum_membrane_) p->PackSummand(data);
  }
}

// Unpack the constituent
void MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::UnpackConstituent(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MixtureConstituent_ElastHyperBase::UnpackConstituent(position, data);

  DRT::ParObject::ExtractfromPack(position, data, current_reference_density_);

  DRT::ParObject::ExtractfromPack(position, data, mue_frac_);

  anisotropyExtension_.UnpackAnisotropy(data, position);

  // loop map of associated potential summands
  for (auto& summand : potsum_membrane_) summand->UnpackSummand(data, position);
}

void MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::RegisterAnisotropyExtensions(
    MAT::Anisotropy& anisotropy)
{
  MixtureConstituent_ElastHyperBase::RegisterAnisotropyExtensions(anisotropy);

  anisotropy.RegisterAnisotropyExtension(anisotropyExtension_);

  for (const auto& summand : potsum_membrane_) summand->RegisterAnisotropyExtensions(anisotropy);
}

// Reads the element from the input file
void MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::ReadElement(
    int numgp, DRT::INPUT::LineDefinition* linedef)
{
  MixtureConstituent_ElastHyperBase::ReadElement(numgp, linedef);

  // Setup summands
  for (const auto& summand : potsum_membrane_) summand->Setup(numgp, linedef);

  current_reference_density_.resize(numgp, params_->RefMassFraction() * InitialRefDensity());
  mue_frac_.resize(numgp, 1.0);
}

// Updates all summands
void MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::Update(
    LINALG::Matrix<3, 3> const& defgrd, Teuchos::ParameterList& params, const int gp,
    const int eleGID)
{
  LINALG::Matrix<1, 3> gprefecoord(true);  // gp coordinates in reference configuration
  gprefecoord = params.get<LINALG::Matrix<1, 3>>("gprefecoord");

  double totaltime = params.get<double>("total time", 10);
  if (totaltime < 0.0)
  {
    dserror("Parameter 'total time' could not be read!");
  }

  current_reference_density_[gp] = params_->RefMassFraction() * InitialRefDensity() *
                                   DRT::Problem::Instance()
                                       ->Funct(params_->damage_function_id_ - 1)
                                       .Evaluate(0, gprefecoord.A(), totaltime);

  MixtureConstituent_ElastHyperBase::Update(defgrd, params, gp, eleGID);

  // loop map of associated potential summands
  for (auto& summand : potsum_membrane_) summand->Update();
}

// Add names for each summand for the quantities for post processing
void MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::VisNames(
    std::map<std::string, int>& names)
{
  MixtureConstituent_ElastHyperBase::VisNames(names);

  // loop map of associated potential summands
  for (auto& summand : potsum_membrane_) summand->VisNames(names);
}

// Add values for each summand of the quantities for post processing
bool MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleGID)
{
  // loop map of associated potential summands
  for (auto& summand : potsum_membrane_)
  {
    bool vis = summand->VisData(name, data, numgp, eleGID);
    if (vis) return true;
  }
  return MixtureConstituent_ElastHyperBase::VisData(name, data, numgp, eleGID);
}

// Returns the reference mass fraction of the constituent
double MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::CurrentRefDensity(int gp) const
{
  return current_reference_density_[gp];
}

void MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::PreEvaluate(
    MixtureRule& mixtureRule, Teuchos::ParameterList& params, int gp, int eleGID)
{
  MIXTURE::MixtureConstituent_ElastHyperBase::PreEvaluate(mixtureRule, params, gp, eleGID);

  // Evaluate mue frac
  Teuchos::RCP<ElastinMembranePrestressStrategy> strategy =
      Teuchos::rcp_dynamic_cast<ElastinMembranePrestressStrategy>(params_->PrestressStrategy());

  if (Teuchos::is_null(strategy))
  {
    dserror(
        "The used prestretch strategy is not compatible with elastin materials. It has to "
        "implement MIXTURE::ElastinMembranePrestressStrategy.");
  }

  mue_frac_[gp] = strategy->EvaluateMueFrac(mixtureRule,
      CylinderCoordinateSystemAnisotropyExtension().GetCylinderCoordinateSystem(gp), *this, *this,
      params, gp, eleGID);
}

void MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::Evaluate(const LINALG::Matrix<3, 3>& F,
    const LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  dserror("This constituent does not support Evaluation without an elastic part.");
}

void MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::EvaluateElasticPart(
    const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<3, 3>& iFextin,
    Teuchos::ParameterList& params, LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
    int gp, int eleGID)
{
  // Compute total inelastic deformation gradient
  static LINALG::Matrix<3, 3> iFin(false);
  iFin.MultiplyNN(iFextin, PrestretchTensor(gp));

  // Evaluate 3D elastic part
  MAT::ElastHyperEvaluateElasticPart(
      F, iFin, S_stress, cmat, Summands(), SummandProperties(), gp, eleGID);

  S_stress.Scale(CurrentRefDensity(gp));
  cmat.Scale(CurrentRefDensity(gp));

  // Evaluate Membrane
  static LINALG::Matrix<6, 1> Smembrane_stress(false);
  static LINALG::Matrix<6, 6> cmatmembrane(false);
  EvaluateStressCMatMembrane(F, iFin, params, Smembrane_stress, cmatmembrane, gp, eleGID);

  S_stress.Update(1.0, Smembrane_stress, 1.0);
  cmat.Update(1.0, cmatmembrane, 1.0);
}

void MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::EvaluateMembraneStress(
    LINALG::Matrix<6, 1>& S, Teuchos::ParameterList& params, int gp, int eleGID)
{
  LINALG::Matrix<6, 6> cmat(false);
  LINALG::Matrix<3, 3> Id(false);
  LINALG::Matrix<3, 3> iFin(false);

  MAT::IdentityMatrix(Id);
  iFin.MultiplyNN(Id, PrestretchTensor(gp));

  EvaluateStressCMatMembrane(Id, iFin, params, S, cmat, gp, eleGID);
}

void MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::EvaluateStressCMatMembrane(
    const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<3, 3>& iFin, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat, int gp, int eleGID) const
{
  static LINALG::Matrix<3, 3> Ce(false);
  MAT::EvaluateCe(F, iFin, Ce);

  // Compute structural tensors in grown configuration
  static LINALG::Matrix<3, 3> Aradgr(false);
  static LINALG::Matrix<3, 3> Aorthgr(false);
  EvaluateStructuralTensorsInGrownConfiguration(Aradgr, Aorthgr, iFin, gp, eleGID);

  static LINALG::Matrix<3, 3> AorthgrCeAorthgrArad(false);
  EvaluateAorthgrCeAorthgrArad(AorthgrCeAorthgrArad, Aradgr, Aorthgr, Ce);
  double detX = AorthgrCeAorthgrArad.Determinant();

  static LINALG::Matrix<3, 3> iFinAorthgriFinT(false);
  EvaluateiFinAorthgriFinT(iFinAorthgriFinT, iFin, Aorthgr);

  // Z = F_{in}^{-T}*A_{gr}^{T}*X^{-1}*A_{gr}^{T}*F_{in}^{-1}
  static LINALG::Matrix<3, 3> iFinTAorthgrTiXTAorthgriFin_sym(false);
  static LINALG::Matrix<6, 1> iFinTAorthgrTiXTAorthgriFin_sym_stress(false);

  EvaluateiFinTAorthgrTiXTAorthgriFin(
      iFinTAorthgrTiXTAorthgriFin_sym, AorthgrCeAorthgrArad, iFin, Aorthgr);

  UTILS::VOIGT::Stresses::MatrixToVector(
      iFinTAorthgrTiXTAorthgriFin_sym, iFinTAorthgrTiXTAorthgriFin_sym_stress);

  // Get material parameter
  double mue = 0.0;
  for (const auto& summand : potsum_membrane_)
  {
    mue += summand->Mue();
  }

  // Compute membrane stress
  static LINALG::Matrix<3, 3> Smembrane;
  Smembrane.Update(CurrentRefDensity(gp) * mue * mue_frac_[gp], iFinAorthgriFinT, 0.0);
  Smembrane.Update(
      -CurrentRefDensity(gp) * mue * mue_frac_[gp] / detX, iFinTAorthgrTiXTAorthgriFin_sym, 1.0);

  UTILS::VOIGT::Stresses::MatrixToVector(Smembrane, S_stress);

  // Compute constitutive tensor
  static LINALG::Matrix<6, 6> dAradgriXAradgr_symdC(false);
  dAradgriXAradgr_symdC.Clear();

  MAT::AddtoCmatHolzapfelProduct(
      dAradgriXAradgr_symdC, iFinTAorthgrTiXTAorthgriFin_sym_stress, -2.0);

  cmat.MultiplyNT(2.0 * CurrentRefDensity(gp) * mue * mue_frac_[gp] / detX,
      iFinTAorthgrTiXTAorthgriFin_sym_stress, iFinTAorthgrTiXTAorthgriFin_sym_stress, 0.0);
  cmat.Update(-CurrentRefDensity(gp) * mue * mue_frac_[gp] / detX, dAradgriXAradgr_symdC, 1.0);
}

void MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::
    EvaluateStructuralTensorsInGrownConfiguration(LINALG::Matrix<3, 3>& Aradgr,
        LINALG::Matrix<3, 3>& Aorthgr, const LINALG::Matrix<3, 3>& iFin, const int gp,
        const int eleGID) const
{
  // Compute inelastic right Cauchy-Green deformation gradient
  static LINALG::Matrix<3, 3> iCin(false);
  iCin.MultiplyNT(iFin, iFin);

  static LINALG::Matrix<3, 3> Fin(false);
  Fin.Invert(iFin);

  // Compute radial structural tensor in grown configuration
  static LINALG::Matrix<3, 3> FinArad(false);
  FinArad.MultiplyNN(Fin, anisotropyExtension_.GetStructuralTensor(gp, 0));
  Aradgr.MultiplyNT(iCin.Dot(anisotropyExtension_.GetStructuralTensor(gp, 0)), FinArad, Fin, 0.0);

  // Compute orthogonal (to radial) structural tensor in grown configuration
  MAT::IdentityMatrix(Aorthgr);
  Aorthgr.Update(-1.0, Aradgr, 1.0);
}

void MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::EvaluateAorthgrCeAorthgrArad(
    LINALG::Matrix<3, 3>& AorthgrCeAorthgrArad, const LINALG::Matrix<3, 3>& Aradgr,
    const LINALG::Matrix<3, 3>& Aorthgr, const LINALG::Matrix<3, 3>& Ce)
{
  static LINALG::Matrix<3, 3> AorthgrCe(false);
  AorthgrCe.MultiplyNN(Aorthgr, Ce);
  AorthgrCeAorthgrArad.MultiplyNN(AorthgrCe, Aorthgr);
  AorthgrCeAorthgrArad.Update(1.0, Aradgr, 1.0);
}

void MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::EvaluateiFinAorthgriFinT(
    LINALG::Matrix<3, 3>& iFinAorthgriFinT, const LINALG::Matrix<3, 3>& iFin,
    const LINALG::Matrix<3, 3>& Aorthgr)
{
  static LINALG::Matrix<3, 3> iFinAorthgr(false);
  iFinAorthgr.MultiplyNN(iFin, Aorthgr);
  iFinAorthgriFinT.MultiplyNT(iFinAorthgr, iFin);
}

void MIXTURE::MixtureConstituent_ElastHyperElastinMembrane::EvaluateiFinTAorthgrTiXTAorthgriFin(
    LINALG::Matrix<3, 3>& iFinTAorthgrTiXTAorthgriFin_sym,
    const LINALG::Matrix<3, 3>& AorthgrCeAorthgrArad, const LINALG::Matrix<3, 3>& iFin,
    const LINALG::Matrix<3, 3>& Aorthgr)
{
  static LINALG::Matrix<3, 3> iAorthgrCeAorthgrArad(false);
  iAorthgrCeAorthgrArad.Invert(AorthgrCeAorthgrArad);

  static LINALG::Matrix<3, 3> AorthgriFin(false);
  AorthgriFin.MultiplyNN(Aorthgr, iFin);

  static LINALG::Matrix<3, 3> iFinTAorthgrTiXT(false);
  iFinTAorthgrTiXT.MultiplyTT(AorthgriFin, iAorthgrCeAorthgrArad);

  static LINALG::Matrix<3, 3> iFinTAorthgrTiXTAorthgriFin(false);
  iFinTAorthgrTiXTAorthgriFin.MultiplyNN(iFinTAorthgrTiXT, AorthgriFin);

  iFinTAorthgrTiXTAorthgriFin_sym.Update(0.5, iFinTAorthgrTiXTAorthgriFin, 0.0);
  iFinTAorthgrTiXTAorthgriFin_sym.UpdateT(0.5, iFinTAorthgrTiXTAorthgriFin, 1.0);
}