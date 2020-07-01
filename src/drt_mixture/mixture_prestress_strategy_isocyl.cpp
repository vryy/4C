/*----------------------------------------------------------------------*/
/*! \file

\brief Prestress strategy for isotropic materials used in a growth remodel simulation

\level 3


*/
/*----------------------------------------------------------------------*/
#include "mixture_prestress_strategy_isocyl.H"
#include "../drt_mat/matpar_bundle.H"
#include "mixture_constituent_elasthyper.H"
#include "../drt_matelast/elast_isoneohooke.H"
#include "../drt_matelast/elast_volsussmanbathe.H"
#include "../drt_mat/anisotropy.H"
#include "../drt_mat/material_service.H"
#include "mixture_rule.H"

MIXTURE::PAR::IsotropicCylinderPrestressStrategy::IsotropicCylinderPrestressStrategy(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : PrestressStrategy(matdata),
      inner_radius_(matdata->GetDouble("INNER_RADIUS")),
      wall_thickness_(matdata->GetDouble("WALL_THICKNESS")),
      axial_prestretch_(matdata->GetDouble("AXIAL_PRESTRETCH")),
      circumferential_prestretch_(matdata->GetDouble("CIRCUMFERENTIAL_PRESTRETCH")),
      pressure_(matdata->GetDouble("PRESSURE"))
{
}

Teuchos::RCP<MIXTURE::PrestressStrategy>
MIXTURE::PAR::IsotropicCylinderPrestressStrategy::CreatePrestressStrategy()
{
  return Teuchos::rcp(new MIXTURE::IsotropicCylinderPrestressStrategy(this));
}

MIXTURE::IsotropicCylinderPrestressStrategy::IsotropicCylinderPrestressStrategy(
    MIXTURE::PAR::IsotropicCylinderPrestressStrategy* params)
    : PrestressStrategy(params), params_(params)
{
}

void MIXTURE::IsotropicCylinderPrestressStrategy::EvaluatePrestress(
    const MAT::CylinderCoordinateSystemProvider& cosy, MIXTURE::MixtureConstituent& constituent,
    LINALG::Matrix<3, 3>& G, Teuchos::ParameterList& params, int gp, int eleGID)
{
  // We evaluate the stress in the reference configuration with a prestretch. Hence, the
  // deformation gradient is the identity matrix and the inverse inelastic deformation gradient ist
  // the prestretch. This results that the 2. Piola-Kirchhoff stress tensor is the same as the
  // Cauchy-stress tensor.

  // Currently, this prestressing technique implements only a certain material (IsoNeoHooke with
  // Sussman-Bathe penalty)

  auto& elhyper = dynamic_cast<MIXTURE::MixtureConstituent_ElastHyperBase&>(constituent);

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

  LINALG::Matrix<1, 3> gprefecoord(true);  // gp coordinates in reference configuration
  gprefecoord = params.get<LINALG::Matrix<1, 3>>("gprefecoord");

  double r = 0;
  for (unsigned i = 0; i < 3; ++i)
  {
    r += cosy.GetRad()(i) * gprefecoord(i);
  }

  double Res = 1.0;
  double dResdlamb_pre;
  double lamb_pre = 1. / (params_->circumferential_prestretch_ * params_->axial_prestretch_);
  while (std::abs(Res) > 1.0e-10)
  {
    Res =
        constituent.InitialConstituentRefDensity() * matiso->Mue() *
            std::pow(params_->circumferential_prestretch_ * params_->axial_prestretch_ * lamb_pre,
                -4. / 3.) *  // TODO: When deriving these equations by hand, I get -2.0 / 3.0. To be
                             //  compatible with the old implementation I decided for now to keep
                             //  this here. This has to be verified later.
            (lamb_pre * lamb_pre -
                (1. / 3.) *
                    (params_->circumferential_prestretch_ * params_->circumferential_prestretch_ +
                        params_->axial_prestretch_ * params_->axial_prestretch_ +
                        lamb_pre * lamb_pre)) +
        constituent.InitialConstituentRefDensity() * matvol->Kappa() *
            ((params_->circumferential_prestretch_ * params_->axial_prestretch_ * lamb_pre) *
                    (params_->circumferential_prestretch_ * params_->axial_prestretch_ * lamb_pre) -
                (params_->circumferential_prestretch_ * params_->axial_prestretch_ * lamb_pre)) +
        ((1.0 - (r - params_->inner_radius_) / params_->wall_thickness_) * params_->pressure_);

    dResdlamb_pre =
        matiso->Mue() * constituent.InitialConstituentRefDensity() *
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
        matiso->Mue() * constituent.InitialConstituentRefDensity() *
            std::pow(params_->circumferential_prestretch_ * params_->circumferential_prestretch_ *
                         params_->axial_prestretch_ * params_->axial_prestretch_ * lamb_pre *
                         lamb_pre,
                -2. / 3.) *
            (2.0 * lamb_pre - (1. / 3.) * (2.0 * lamb_pre)) +
        matvol->Kappa() * constituent.InitialConstituentRefDensity() *
            (2.0 * (params_->circumferential_prestretch_ * params_->axial_prestretch_ * lamb_pre) *
                    params_->circumferential_prestretch_ * params_->axial_prestretch_ -
                params_->circumferential_prestretch_ * params_->axial_prestretch_);

    lamb_pre = lamb_pre + (-Res / dResdlamb_pre);
  }

  // Build prestretch tensor
  G.MultiplyNT(lamb_pre, cosy.GetRad(), cosy.GetRad(), 0.0);
  G.MultiplyNT(params_->axial_prestretch_, cosy.GetAxi(), cosy.GetAxi(), 1.0);
  G.MultiplyNT(params_->circumferential_prestretch_, cosy.GetCir(), cosy.GetCir(), 1.0);
}

double MIXTURE::IsotropicCylinderPrestressStrategy::EvaluateMueFrac(MixtureRule& mixtureRule,
    const MAT::CylinderCoordinateSystemProvider& cosy, MIXTURE::MixtureConstituent& constituent,
    ElastinMembraneEvaluation& membraneEvaluation, Teuchos::ParameterList& params, int gp,
    int eleGID) const
{
  LINALG::Matrix<3, 3> F(false);
  LINALG::Matrix<6, 1> E_strain(true);
  LINALG::Matrix<6, 1> S_stress(true);
  LINALG::Matrix<6, 6> cmat(true);
  MAT::IdentityMatrix(F);


  mixtureRule.Evaluate(F, E_strain, params, S_stress, cmat, gp, eleGID);

  LINALG::Matrix<6, 1> Acir(false);
  // Compute structural tensor
  for (int i = 0; i < 3; ++i) Acir(i) = cosy.GetCir()(i) * cosy.GetCir()(i);
  Acir(3) = 2.0 * cosy.GetCir()(0) * cosy.GetCir()(1);
  Acir(4) = 2.0 * cosy.GetCir()(1) * cosy.GetCir()(2);
  Acir(5) = 2.0 * cosy.GetCir()(0) * cosy.GetCir()(2);

  LINALG::Matrix<6, 1> Smembrane(false);
  membraneEvaluation.EvaluateMembraneStress(Smembrane, params, gp, eleGID);

  double total_stress = S_stress.Dot(Acir);      // stress of all constituents in circular direction
  double membrane_stress = Smembrane.Dot(Acir);  // stress of the membrane in circular direction

  // Compute stress as a result of Barlow's formula ("Kesselformel")
  double target_stress = (params_->pressure_ * params_->inner_radius_) /
                         params_->wall_thickness_;  // stress that we need in circular direction

  return (target_stress - (total_stress - membrane_stress)) / membrane_stress;
}