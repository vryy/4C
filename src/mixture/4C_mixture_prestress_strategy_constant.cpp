/*----------------------------------------------------------------------*/
/*! \file

\brief Constant prestretch strategy

\level 3


*/
/*----------------------------------------------------------------------*/
#include "4C_mixture_prestress_strategy_constant.hpp"

#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_mat_anisotropy.hpp"
#include "4C_mat_anisotropy_coordinate_system_provider.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"
#include "4C_matelast_isoneohooke.hpp"
#include "4C_matelast_volsussmanbathe.hpp"
#include "4C_mixture_constituent_elasthyper.hpp"
#include "4C_mixture_rule.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

MIXTURE::PAR::ConstantPrestressStrategy::ConstantPrestressStrategy(
    const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata)
    : PrestressStrategy(matdata), prestretch_()
{
  std::copy_n(matdata->Get<std::vector<double>>("PRESTRETCH").begin(), 9, prestretch_.begin());
}

std::unique_ptr<MIXTURE::PrestressStrategy>
MIXTURE::PAR::ConstantPrestressStrategy::create_prestress_strategy()
{
  std::unique_ptr<MIXTURE::PrestressStrategy> prestressStrategy(
      new MIXTURE::ConstantPrestressStrategy(this));
  return prestressStrategy;
}

MIXTURE::ConstantPrestressStrategy::ConstantPrestressStrategy(
    MIXTURE::PAR::ConstantPrestressStrategy* params)
    : PrestressStrategy(params), params_(params)
{
}

void MIXTURE::ConstantPrestressStrategy::Setup(
    MIXTURE::MixtureConstituent& constituent, Teuchos::ParameterList& params, int numgp, int eleGID)
{
  // nothing to do
}

void MIXTURE::ConstantPrestressStrategy::EvaluatePrestress(const MixtureRule& mixtureRule,
    const Teuchos::RCP<const MAT::CoordinateSystemProvider> cosy,
    MIXTURE::MixtureConstituent& constituent, CORE::LINALG::Matrix<3, 3>& G,
    Teuchos::ParameterList& params, int gp, int eleGID)
{
  // setup prestretch
  const CORE::LINALG::Matrix<9, 1> prestretch_vector(params_->prestretch_.data(), true);

  CORE::LINALG::VOIGT::Matrix9x1to3x3(prestretch_vector, G);
}

void MIXTURE::ConstantPrestressStrategy::Update(
    const Teuchos::RCP<const MAT::CoordinateSystemProvider> anisotropy,
    MIXTURE::MixtureConstituent& constituent, const CORE::LINALG::Matrix<3, 3>& F,
    CORE::LINALG::Matrix<3, 3>& G, Teuchos::ParameterList& params, int gp, int eleGID)
{
}
FOUR_C_NAMESPACE_CLOSE
