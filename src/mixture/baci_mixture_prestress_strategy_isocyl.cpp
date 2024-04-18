/*----------------------------------------------------------------------*/
/*! \file

\brief Prestress strategy for isotropic materials used in a growth remodel simulation

\level 3


*/
/*----------------------------------------------------------------------*/
#include "baci_mixture_prestress_strategy_isocyl.hpp"

#include "baci_linalg_fixedsizematrix_generators.hpp"
#include "baci_mat_anisotropy.hpp"
#include "baci_mat_anisotropy_coordinate_system_provider.hpp"
#include "baci_mat_par_bundle.hpp"
#include "baci_mat_service.hpp"
#include "baci_matelast_isoneohooke.hpp"
#include "baci_matelast_volsussmanbathe.hpp"
#include "baci_mixture_constituent_elasthyper.hpp"
#include "baci_mixture_rule.hpp"
#include "baci_mixture_rule_growthremodel.hpp"

FOUR_C_NAMESPACE_OPEN

MIXTURE::PAR::IsotropicCylinderPrestressStrategy::IsotropicCylinderPrestressStrategy(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : PrestressStrategy(matdata),
      inner_radius_(*matdata->Get<double>("INNER_RADIUS")),
      wall_thickness_(*matdata->Get<double>("WALL_THICKNESS")),
      axial_prestretch_(*matdata->Get<double>("AXIAL_PRESTRETCH")),
      circumferential_prestretch_(*matdata->Get<double>("CIRCUMFERENTIAL_PRESTRETCH")),
      pressure_(*matdata->Get<double>("PRESSURE"))
{
}

std::unique_ptr<MIXTURE::PrestressStrategy>
MIXTURE::PAR::IsotropicCylinderPrestressStrategy::CreatePrestressStrategy()
{
  std::unique_ptr<MIXTURE::PrestressStrategy> prestressStrategy(
      new MIXTURE::IsotropicCylinderPrestressStrategy(this));
  return prestressStrategy;
}

MIXTURE::IsotropicCylinderPrestressStrategy::IsotropicCylinderPrestressStrategy(
    MIXTURE::PAR::IsotropicCylinderPrestressStrategy* params)
    : PrestressStrategy(params), params_(params)
{
}

void MIXTURE::IsotropicCylinderPrestressStrategy::Setup(
    MIXTURE::MixtureConstituent& constituent, Teuchos::ParameterList& params, int numgp, int eleGID)
{
  // nothing to do
}

void MIXTURE::IsotropicCylinderPrestressStrategy::EvaluatePrestress(const MixtureRule& mixtureRule,
    const Teuchos::RCP<const MAT::CoordinateSystemProvider> cosy,
    MIXTURE::MixtureConstituent& constituent, CORE::LINALG::Matrix<3, 3>& G,
    Teuchos::ParameterList& params, int gp, int eleGID)
{
  // We evaluate the stress in the reference configuration with a prestretch. Hence, the
  // deformation gradient is the identity matrix and the inverse inelastic deformation gradient ist
  // the prestretch. This results that the 2. Piola-Kirchhoff stress tensor is the same as the
  // Cauchy-stress tensor.

  // Currently, this prestressing technique implements only a certain material (IsoNeoHooke with
  // Sussman-Bathe penalty)

  auto& elhyper = dynamic_cast<MIXTURE::MixtureConstituentElastHyperBase&>(constituent);

  if (elhyper.Summands().size() != 2)
  {
    dserror(
        "Currrently, the prestressing technique is only implemented for an ElastHyper constituent "
        "using an IsoNeoHooke summand with Sussman-Bathe penalty (Hence, exactly 2 summands are "
        "needed)");
  }

  // Let's assume that for simplicity, the first index is the IsoNeoHooke material and the second
  // index is the Sussman-Bathe penalty parameter
  auto matiso = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::IsoNeoHooke>(elhyper.Summands()[0]);
  auto matvol = Teuchos::rcp_dynamic_cast<MAT::ELASTIC::VolSussmanBathe>(elhyper.Summands()[1]);

  if (Teuchos::is_null(matiso))
  {
    dserror(
        "The first summand of the constituent needs to be an IsoNeoHooke material law. This is a "
        "requirement from the prestressing technique.");
  }

  if (Teuchos::is_null(matvol))
  {
    dserror(
        "The second summand of the constituent needs to be a Sussman-Bathe penalty term. This is a "
        "requirement from the prestressing technique.");
  }

  // This prestress strategy is only valid for G&R simulations
  const auto& growth_remodel_rule =
      dynamic_cast<const MIXTURE::GrowthRemodelMixtureRule&>(mixtureRule);


  Teuchos::RCP<const MAT::CylinderCoordinateSystemProvider> cylinderCosy =
      cosy->GetCylinderCoordinateSystem();

  if (Teuchos::is_null(cylinderCosy))
  {
    dserror(
        "No cylinder coordinate system is defined but required by the cylinder prestress "
        "strategy!");
  }

  const auto& reference_coordinates = params.get<CORE::LINALG::Matrix<3, 1>>("gp_coords_ref");

  double r = 0;
  for (unsigned i = 0; i < 3; ++i)
  {
    r += cylinderCosy->GetRad()(i) * reference_coordinates(i);
  }

  double initial_constituent_reference_density =
      growth_remodel_rule.GetConstituentInitialReferenceMassDensity(constituent);

  double Res = 1.0;
  double dResdlamb_pre;
  double lamb_pre = 1. / (params_->circumferential_prestretch_ * params_->axial_prestretch_);
  while (std::abs(Res) > 1.0e-10)
  {
    Res =
        matiso->Mue() * initial_constituent_reference_density *
            std::pow(params_->circumferential_prestretch_ * params_->axial_prestretch_ * lamb_pre,
                -4. / 3.) *  // TODO: When deriving these equations by hand, I get -2.0 / 3.0. To be
                             //  compatible with the old implementation I decided for now to keep
                             //  this here. This has to be verified later.
            (lamb_pre * lamb_pre -
                (1. / 3.) *
                    (params_->circumferential_prestretch_ * params_->circumferential_prestretch_ +
                        params_->axial_prestretch_ * params_->axial_prestretch_ +
                        lamb_pre * lamb_pre)) +
        matvol->Kappa() * initial_constituent_reference_density *
            ((params_->circumferential_prestretch_ * params_->axial_prestretch_ * lamb_pre) *
                    (params_->circumferential_prestretch_ * params_->axial_prestretch_ * lamb_pre) -
                (params_->circumferential_prestretch_ * params_->axial_prestretch_ * lamb_pre)) +
        ((1.0 - (r - params_->inner_radius_) / params_->wall_thickness_) * params_->pressure_);

    dResdlamb_pre =
        matiso->Mue() * initial_constituent_reference_density *
            (-(4. / 3.) *
                std::pow(
                    params_->circumferential_prestretch_ * params_->axial_prestretch_ * lamb_pre,
                    -7. / 3.) *
                params_->circumferential_prestretch_ * params_->axial_prestretch_) *
            (lamb_pre * lamb_pre -
                (1. / 3.) *
                    (params_->circumferential_prestretch_ * params_->circumferential_prestretch_ +
                        params_->axial_prestretch_ * params_->axial_prestretch_ +
                        lamb_pre * lamb_pre)) +
        matiso->Mue() * initial_constituent_reference_density *
            std::pow(params_->circumferential_prestretch_ * params_->circumferential_prestretch_ *
                         params_->axial_prestretch_ * params_->axial_prestretch_ * lamb_pre *
                         lamb_pre,
                -2. / 3.) *
            (2.0 * lamb_pre - (1. / 3.) * (2.0 * lamb_pre)) +
        matvol->Kappa() * initial_constituent_reference_density *
            (2.0 * (params_->circumferential_prestretch_ * params_->axial_prestretch_ * lamb_pre) *
                    params_->circumferential_prestretch_ * params_->axial_prestretch_ -
                params_->circumferential_prestretch_ * params_->axial_prestretch_);

    lamb_pre = lamb_pre + (-Res / dResdlamb_pre);
  }

  // Build prestretch tensor
  G.MultiplyNT(lamb_pre, cylinderCosy->GetRad(), cylinderCosy->GetRad(), 0.0);
  G.MultiplyNT(params_->axial_prestretch_, cylinderCosy->GetAxi(), cylinderCosy->GetAxi(), 1.0);
  G.MultiplyNT(
      params_->circumferential_prestretch_, cylinderCosy->GetCir(), cylinderCosy->GetCir(), 1.0);
}

double MIXTURE::IsotropicCylinderPrestressStrategy::EvaluateMueFrac(MixtureRule& mixtureRule,
    const Teuchos::RCP<const MAT::CoordinateSystemProvider> cosy,
    MIXTURE::MixtureConstituent& constituent, ElastinMembraneEvaluation& membraneEvaluation,
    Teuchos::ParameterList& params, int gp, int eleGID) const
{
  Teuchos::RCP<const MAT::CylinderCoordinateSystemProvider> cylinderCosy =
      cosy->GetCylinderCoordinateSystem();

  if (Teuchos::is_null(cylinderCosy))
  {
    dserror(
        "No cylinder coordinate system is defined but required by the cylinder prestress "
        "strategy!");
  }

  CORE::LINALG::Matrix<3, 3> F = CORE::LINALG::IdentityMatrix<3>();
  CORE::LINALG::Matrix<6, 1> E_strain(true);
  CORE::LINALG::Matrix<6, 1> S_stress(true);
  CORE::LINALG::Matrix<6, 6> cmat(true);


  mixtureRule.Evaluate(F, E_strain, params, S_stress, cmat, gp, eleGID);

  CORE::LINALG::Matrix<6, 1> Acir(false);
  // Compute structural tensor
  for (int i = 0; i < 3; ++i) Acir(i) = cylinderCosy->GetCir()(i) * cylinderCosy->GetCir()(i);
  Acir(3) = 2.0 * cylinderCosy->GetCir()(0) * cylinderCosy->GetCir()(1);
  Acir(4) = 2.0 * cylinderCosy->GetCir()(1) * cylinderCosy->GetCir()(2);
  Acir(5) = 2.0 * cylinderCosy->GetCir()(0) * cylinderCosy->GetCir()(2);


  // This prestress strategy is only valid for G&R simulations
  const auto& growth_remodel_rule =
      dynamic_cast<const MIXTURE::GrowthRemodelMixtureRule&>(mixtureRule);
  double initial_constituent_reference_density =
      growth_remodel_rule.GetConstituentInitialReferenceMassDensity(constituent);

  CORE::LINALG::Matrix<6, 1> Smembrane(false);
  membraneEvaluation.EvaluateMembraneStress(Smembrane, params, gp, eleGID);
  Smembrane.Scale(initial_constituent_reference_density);

  double total_stress = S_stress.Dot(Acir);      // stress of all constituents in circular direction
  double membrane_stress = Smembrane.Dot(Acir);  // stress of the membrane in circular direction

  // Compute stress as a result of Barlow's formula ("Kesselformel")
  double target_stress = (params_->pressure_ * params_->inner_radius_) /
                         params_->wall_thickness_;  // stress that we need in circular direction

  return (target_stress - (total_stress - membrane_stress)) / membrane_stress;
}

void MIXTURE::IsotropicCylinderPrestressStrategy::Update(
    const Teuchos::RCP<const MAT::CoordinateSystemProvider> anisotropy,
    MIXTURE::MixtureConstituent& constituent, const CORE::LINALG::Matrix<3, 3>& F,
    CORE::LINALG::Matrix<3, 3>& G, Teuchos::ParameterList& params, int gp, int eleGID)
{
}
FOUR_C_NAMESPACE_CLOSE
