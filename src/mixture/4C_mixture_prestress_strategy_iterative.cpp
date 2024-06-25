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
    const Core::Mat::PAR::Parameter::Data& matdata)
    : PrestressStrategy(matdata),
      isochoric_(matdata.parameters.get<bool>("ISOCHORIC")),
      is_active_(matdata.parameters.get<bool>("ACTIVE"))
{
}

std::unique_ptr<MIXTURE::PrestressStrategy>
MIXTURE::PAR::IterativePrestressStrategy::create_prestress_strategy()
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

void MIXTURE::IterativePrestressStrategy::setup(
    MIXTURE::MixtureConstituent& constituent, Teuchos::ParameterList& params, int numgp, int eleGID)
{
  // nothing to do
}

void MIXTURE::IterativePrestressStrategy::evaluate_prestress(const MixtureRule& mixtureRule,
    const Teuchos::RCP<const Mat::CoordinateSystemProvider> anisotropy,
    MIXTURE::MixtureConstituent& constituent, Core::LinAlg::Matrix<3, 3>& G,
    Teuchos::ParameterList& params, int gp, int eleGID)
{
  // Start with zero prestretch
  G = Core::LinAlg::IdentityMatrix<3>();
}

void MIXTURE::IterativePrestressStrategy::update(
    const Teuchos::RCP<const Mat::CoordinateSystemProvider> anisotropy,
    MIXTURE::MixtureConstituent& constituent, const Core::LinAlg::Matrix<3, 3>& F,
    Core::LinAlg::Matrix<3, 3>& G, Teuchos::ParameterList& params, int gp, int eleGID)
{
  // only update prestress if it is active
  if (!params_->is_active_) return;

  // Compute isochoric part of the deformation
  Core::LinAlg::Matrix<3, 3> F_bar;
  if (params_->isochoric_)
  {
    F_bar.update(std::pow(F.determinant(), -1.0 / 3.0), F);
  }
  else
  {
    F_bar.update(F);
  }

  // Compute new predeformation gradient
  Core::LinAlg::Matrix<3, 3> G_old(G);
  G.multiply_nn(F_bar, G_old);


  // Compute polar decomposition of the prestretch deformation gradient

  // Singular value decomposition of F = RU
  Core::LinAlg::Matrix<3, 3> Q(true);
  Core::LinAlg::Matrix<3, 3> S(true);
  Core::LinAlg::Matrix<3, 3> VT(true);

  Core::LinAlg::SVD<3, 3>(G, Q, S, VT);

  // Compute stretch tensor G = U = V * S * VT
  Core::LinAlg::Matrix<3, 3> VS;

  VS.multiply_tn(VT, S);
  G.multiply_nn(VS, VT);
}
FOUR_C_NAMESPACE_CLOSE
