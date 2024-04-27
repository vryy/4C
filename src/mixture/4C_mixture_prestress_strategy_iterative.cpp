/*----------------------------------------------------------------------*/
/*! \file

\brief Prestress strategy for isotropic materials used in a growth remodel simulation

\level 3

*/
/*----------------------------------------------------------------------*/
#include "4C_mixture_prestress_strategy_iterative.hpp"

#include "4C_linalg_fixedsizematrix_generators.hpp"
#include "4C_linalg_utils_densematrix_svd.hpp"
#include "4C_mat_anisotropy.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"
#include "4C_matelast_isoneohooke.hpp"
#include "4C_matelast_volsussmanbathe.hpp"
#include "4C_mixture_constituent_elasthyper.hpp"
#include "4C_mixture_rule.hpp"

FOUR_C_NAMESPACE_OPEN

MIXTURE::PAR::IterativePrestressStrategy::IterativePrestressStrategy(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : PrestressStrategy(matdata),
      isochoric_(matdata->Get<bool>("ISOCHORIC")),
      is_active_(matdata->Get<bool>("ACTIVE"))
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
    const Teuchos::RCP<const MAT::CoordinateSystemProvider> anisotropy,
    MIXTURE::MixtureConstituent& constituent, CORE::LINALG::Matrix<3, 3>& G,
    Teuchos::ParameterList& params, int gp, int eleGID)
{
  // Start with zero prestretch
  G = CORE::LINALG::IdentityMatrix<3>();
}

void MIXTURE::IterativePrestressStrategy::Update(
    const Teuchos::RCP<const MAT::CoordinateSystemProvider> anisotropy,
    MIXTURE::MixtureConstituent& constituent, const CORE::LINALG::Matrix<3, 3>& F,
    CORE::LINALG::Matrix<3, 3>& G, Teuchos::ParameterList& params, int gp, int eleGID)
{
  // only update prestress if it is active
  if (!params_->is_active_) return;

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
FOUR_C_NAMESPACE_CLOSE
