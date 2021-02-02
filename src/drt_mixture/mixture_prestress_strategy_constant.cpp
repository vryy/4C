/*----------------------------------------------------------------------*/
/*! \file

\brief Constant prestretch strategy

\level 3


*/
/*----------------------------------------------------------------------*/
#include "mixture_prestress_strategy_constant.H"
#include "../drt_mat/matpar_bundle.H"
#include "mixture_constituent_elasthyper.H"
#include "../drt_matelast/elast_isoneohooke.H"
#include "../drt_matelast/elast_volsussmanbathe.H"
#include "../drt_mat/anisotropy.H"
#include "../drt_mat/material_service.H"
#include "../drt_mat/anisotropy_coordinate_system_provider.H"
#include "../drt_lib/voigt_notation.H"
#include "mixture_rule.H"

MIXTURE::PAR::ConstantPrestressStrategy::ConstantPrestressStrategy(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : PrestressStrategy(matdata), prestretch_()
{
  std::copy_n(matdata->Get<std::vector<double>>("PRESTRETCH")->begin(), 9, prestretch_.begin());
}

Teuchos::RCP<MIXTURE::PrestressStrategy>
MIXTURE::PAR::ConstantPrestressStrategy::CreatePrestressStrategy()
{
  return Teuchos::rcp(new MIXTURE::ConstantPrestressStrategy(this));
}

MIXTURE::ConstantPrestressStrategy::ConstantPrestressStrategy(
    MIXTURE::PAR::ConstantPrestressStrategy* params)
    : PrestressStrategy(params), params_(params)
{
}

void MIXTURE::ConstantPrestressStrategy::EvaluatePrestress(
    const Teuchos::RCP<const MAT::CoordinateSystemProvider> cosy,
    MIXTURE::MixtureConstituent& constituent, LINALG::Matrix<3, 3>& G,
    Teuchos::ParameterList& params, int gp, int eleGID)
{
  // setup prestretch
  const LINALG::Matrix<9, 1> prestretch_vector(params_->prestretch_.data(), true);

  UTILS::VOIGT::Matrix9x1to3x3(prestretch_vector, G);
}

void MIXTURE::ConstantPrestressStrategy::UpdatePrestress(
    const Teuchos::RCP<const MAT::CoordinateSystemProvider> anisotropy,
    MIXTURE::MixtureConstituent& constituent, const LINALG::Matrix<3, 3>& F,
    LINALG::Matrix<3, 3>& G, Teuchos::ParameterList& params, int gp, int eleGID)
{
  dserror(
      "The prestretching strategy that you have chosen does not need iterative prestretching. It "
      "ensures equilibrium during setup.");
}