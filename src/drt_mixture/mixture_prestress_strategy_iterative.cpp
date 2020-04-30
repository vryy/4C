/*----------------------------------------------------------------------*/
/*! \file

\brief Prestress strategy for isotropic materials used in a growth remodel simulation

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/
#include "mixture_prestress_strategy_iterative.H"
#include "../drt_mat/matpar_bundle.H"
#include "mixture_constituent_elasthyper.H"
#include "../drt_matelast/elast_isoneohooke.H"
#include "../drt_matelast/elast_volsussmanbathe.H"
#include "../drt_mat/anisotropy.H"
#include "../drt_mat/material_service.H"
#include "mixture_rule.H"
#include "../linalg/linalg_utils_densematrix_svd.H"

MIXTURE::PAR::IterativePrestressStrategy::IterativePrestressStrategy(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : PrestressStrategy(matdata)
{
}

Teuchos::RCP<MIXTURE::PrestressStrategy>
MIXTURE::PAR::IterativePrestressStrategy::CreatePrestressStrategy()
{
  return Teuchos::rcp(new MIXTURE::IterativePrestressStrategy(this));
}

MIXTURE::IterativePrestressStrategy::IterativePrestressStrategy(
    MIXTURE::PAR::IterativePrestressStrategy* params)
    : PrestressStrategy(params), params_(params)
{
}

void MIXTURE::IterativePrestressStrategy::EvaluatePrestress(
    const MAT::CylinderCoordinateSystemProvider& cosy, MIXTURE::MixtureConstituent& constituent,
    LINALG::Matrix<3, 3>& G, Teuchos::ParameterList& params, int gp, int eleGID)
{
  // Start with zero prestretch
  MAT::IdentityMatrix(G);
}

void MIXTURE::IterativePrestressStrategy::UpdatePrestress(
    const MAT::CylinderCoordinateSystemProvider& anisotropy,
    MIXTURE::MixtureConstituent& constituent, const LINALG::Matrix<3, 3>& F,
    LINALG::Matrix<3, 3>& G, Teuchos::ParameterList& params, int gp, int eleGID)
{
  // 1. Computate polar decomposition

  // 1.1 Singular value decomposition of F = RU
  LINALG::Matrix<3, 3> Q(true);
  LINALG::Matrix<3, 3> S(true);
  LINALG::Matrix<3, 3> VT(true);

  LINALG::SVD<3, 3>(F, Q, S, VT);

  // Compute stretch tensor U = V * S * VT
  LINALG::Matrix<3, 3> VS;
  LINALG::Matrix<3, 3> U;

  VS.MultiplyTN(VT, S);
  U.MultiplyNN(VS, VT);

  // 2. Update G with inverse of stretch
  LINALG::Matrix<3, 3> Gold(G);
  G.MultiplyNN(U, Gold);
}