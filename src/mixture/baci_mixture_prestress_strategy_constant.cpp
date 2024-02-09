/*----------------------------------------------------------------------*/
/*! \file

\brief Constant prestretch strategy

\level 3


*/
/*----------------------------------------------------------------------*/
#include "baci_mixture_prestress_strategy_constant.hpp"

#include "baci_linalg_fixedsizematrix_voigt_notation.hpp"
#include "baci_mat_anisotropy.hpp"
#include "baci_mat_anisotropy_coordinate_system_provider.hpp"
#include "baci_mat_par_bundle.hpp"
#include "baci_mat_service.hpp"
#include "baci_matelast_isoneohooke.hpp"
#include "baci_matelast_volsussmanbathe.hpp"
#include "baci_mixture_constituent_elasthyper.hpp"
#include "baci_mixture_rule.hpp"

#include <memory>

BACI_NAMESPACE_OPEN

MIXTURE::PAR::ConstantPrestressStrategy::ConstantPrestressStrategy(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : PrestressStrategy(matdata), prestretch_()
{
  std::copy_n(matdata->Get<std::vector<double>>("PRESTRETCH")->begin(), 9, prestretch_.begin());
}

std::unique_ptr<MIXTURE::PrestressStrategy>
MIXTURE::PAR::ConstantPrestressStrategy::CreatePrestressStrategy()
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
BACI_NAMESPACE_CLOSE
