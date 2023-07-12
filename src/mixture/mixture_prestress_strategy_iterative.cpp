/*----------------------------------------------------------------------*/
/*! \file

\brief Prestress strategy for isotropic materials used in a growth remodel simulation

\level 3

*/
/*----------------------------------------------------------------------*/
#include "mixture_prestress_strategy_iterative.H"
#include "mat_par_bundle.H"
#include "mixture_constituent_elasthyper.H"
#include "matelast_isoneohooke.H"
#include "matelast_volsussmanbathe.H"
#include "mat_anisotropy.H"
#include "mat_service.H"
#include "mixture_rule.H"
#include "linalg_utils_densematrix_svd.H"

MIXTURE::PAR::IterativePrestressStrategy::IterativePrestressStrategy(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : PrestressStrategy(matdata), isochoric_(static_cast<bool>(matdata->GetInt("ISOCHORIC")))
{
}

std::unique_ptr<MIXTURE::PrestressStrategy>
MIXTURE::PAR::IterativePrestressStrategy::CreatePrestressStrategy()
{
  std::unique_ptr<MIXTURE::PrestressStrategy> prestressStrategy(
      new MIXTURE::IterativePrestressStrategy(this));
  return prestressStrategy;
}

MIXTURE::IterativePrestressStrategy::IterativePrestressStrategy(
    MIXTURE::PAR::IterativePrestressStrategy* params)
    : PrestressStrategy(params), params_(params)
{
}

void MIXTURE::IterativePrestressStrategy::Setup(
    MIXTURE::MixtureConstituent& constituent, Teuchos::ParameterList& params, int numgp, int eleGID)
{
  // nothing to do
}

void MIXTURE::IterativePrestressStrategy::EvaluatePrestress(const MixtureRule& mixtureRule,
    const Teuchos::RCP<const MAT::CoordinateSystemProvider> cosy,
    MIXTURE::MixtureConstituent& constituent, CORE::LINALG::Matrix<3, 3>& G,
    Teuchos::ParameterList& params, int gp, int eleGID)
{
  // Start with zero prestretch
  MAT::IdentityMatrix(G);
}

void MIXTURE::IterativePrestressStrategy::UpdatePrestress(
    const Teuchos::RCP<const MAT::CoordinateSystemProvider> anisotropy,
    MIXTURE::MixtureConstituent& constituent, const CORE::LINALG::Matrix<3, 3>& F,
    CORE::LINALG::Matrix<3, 3>& G, Teuchos::ParameterList& params, int gp, int eleGID)
{
  // Compute isochoric part of the deformation
  CORE::LINALG::Matrix<3, 3> F_bar;
  if (params_->isochoric_)
  {
    F_bar.Update(std::pow(F.Determinant(), -1.0 / 3.0), F);
  }
  else
  {
    F_bar.Update(F);
  }

  // Compute new predeformation gradient
  CORE::LINALG::Matrix<3, 3> G_old(G);
  G.MultiplyNN(F_bar, G_old);


  // Compute polar decomposition of the prestretch deformation gradient

  // Singular value decomposition of F = RU
  CORE::LINALG::Matrix<3, 3> Q(true);
  CORE::LINALG::Matrix<3, 3> S(true);
  CORE::LINALG::Matrix<3, 3> VT(true);

  CORE::LINALG::SVD<3, 3>(G, Q, S, VT);

  // Compute stretch tensor G = U = V * S * VT
  CORE::LINALG::Matrix<3, 3> VS;

  VS.MultiplyTN(VT, S);
  G.MultiplyNN(VS, VT);
}