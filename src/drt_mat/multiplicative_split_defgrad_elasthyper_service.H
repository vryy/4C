/*----------------------------------------------------------------------*/
/*! \file

\brief Common functionalities for materials that base on a multiplicative split of the deformation
gradient.

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_MULTIPLICATIVE_SPLIT_DEFGRAD_ELASTHYPER_SERVICE_H
#define BACI_MULTIPLICATIVE_SPLIT_DEFGRAD_ELASTHYPER_SERVICE_H
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/voigt_notation.H"
#include "../drt_mat/material_service.H"
#include "elasthyper_service.H"

namespace MAT
{
  inline void EvaluateCe(
      const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<3, 3>& iFin, LINALG::Matrix<3, 3>& Ce)
  {
    static LINALG::Matrix<3, 3> FiFin(false);
    FiFin.MultiplyNN(F, iFin);
    Ce.MultiplyTN(FiFin, FiFin);
  }

  inline void EvaluateiCinCiCin(const LINALG::Matrix<3, 3>& C, const LINALG::Matrix<3, 3>& iCin,
      LINALG::Matrix<3, 3>& iCinCiCin)
  {
    static LINALG::Matrix<3, 3> CiCin(false);
    CiCin.MultiplyNN(C, iCin);
    iCinCiCin.MultiplyNN(iCin, CiCin);
  }

  inline void ElastHyperEvaluateElasticPart(const LINALG::Matrix<3, 3>& F,
      const LINALG::Matrix<3, 3>& iFin, LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
      const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum,
      MAT::SummandProperties summandProperties, const int gp, const int eleGID)
  {
    if (summandProperties.anisomod or summandProperties.anisoprinc)
    {
      dserror("An additional inelastic part is not yet implemented for anisotropic materials.");
    }

    S_stress.Clear();
    cmat.Clear();

    // Variables needed for the computation of the stress resultants
    static LINALG::Matrix<3, 3> C(true);
    static LINALG::Matrix<3, 3> Ce(true);
    static LINALG::Matrix<3, 3> iC(true);
    static LINALG::Matrix<3, 3> iCin(true);
    static LINALG::Matrix<3, 3> iCinCiCin(true);

    static LINALG::Matrix<6, 1> iCinv(true);
    static LINALG::Matrix<6, 1> iCinCiCinv(true);
    static LINALG::Matrix<6, 1> iCv(true);
    static LINALG::Matrix<3, 1> principleInvariantsCe(true);

    // Compute right Cauchy-Green tensor C=F^TF
    C.MultiplyTN(F, F);

    // Compute inverse right Cauchy-Green tensor C^-1
    iC.Invert(C);

    // Compute inverse inelastic right Cauchy-Green Tensor
    iCin.MultiplyNT(iFin, iFin);

    // Compute iCin * C * iCin
    MAT::EvaluateiCinCiCin(C, iCin, iCinCiCin);

    // Compute Ce
    MAT::EvaluateCe(F, iFin, Ce);

    // Compute principal invariants
    MAT::InvariantsPrincipal(principleInvariantsCe, Ce);

    LINALG::Matrix<3, 1> dPIe(true);
    LINALG::Matrix<6, 1> ddPIIe(true);

    MAT::ElastHyperEvaluateInvariantDerivatives(
        principleInvariantsCe, dPIe, ddPIIe, potsum, summandProperties, gp, eleGID);

    // 2nd Piola Kirchhoff stress factors (according to Holzapfel-Nonlinear Solid Mechanics p. 216)
    static LINALG::Matrix<3, 1> gamma(true);
    // constitutive tensor factors (according to Holzapfel-Nonlinear Solid Mechanics p. 261)
    static LINALG::Matrix<8, 1> delta(true);

    MAT::CalculateGammaDelta(gamma, delta, principleInvariantsCe, dPIe, ddPIIe);

    // Convert necessary tensors to stress-like Voigt-Notation
    UTILS::VOIGT::Stresses::MatrixToVector(iCin, iCinv);
    UTILS::VOIGT::Stresses::MatrixToVector(iCinCiCin, iCinCiCinv);
    UTILS::VOIGT::Stresses::MatrixToVector(iC, iCv);

    // Contribution to 2nd Piola-Kirchhoff stress tensor
    S_stress.Update(gamma(0), iCinv, 1.0);
    S_stress.Update(gamma(1), iCinCiCinv, 1.0);
    S_stress.Update(gamma(2), iCv, 1.0);

    // Contribution to the linearization
    cmat.MultiplyNT(delta(0), iCinv, iCinv, 1.);
    cmat.MultiplyNT(delta(1), iCinCiCinv, iCinv, 1.);
    cmat.MultiplyNT(delta(1), iCinv, iCinCiCinv, 1.);
    cmat.MultiplyNT(delta(2), iCinv, iCv, 1.);
    cmat.MultiplyNT(delta(2), iCv, iCinv, 1.);
    cmat.MultiplyNT(delta(3), iCinCiCinv, iCinCiCinv, 1.);
    cmat.MultiplyNT(delta(4), iCinCiCinv, iCv, 1.);
    cmat.MultiplyNT(delta(4), iCv, iCinCiCinv, 1.);
    cmat.MultiplyNT(delta(5), iCv, iCv, 1.);
    MAT::AddtoCmatHolzapfelProduct(cmat, iCv, delta(6));
    MAT::AddtoCmatHolzapfelProduct(cmat, iCinv, delta(7));
  }

}  // namespace MAT
#endif  // BACI_MULTIPLICATIVE_SPLIT_DEFGRAD_ELASTHYPER_SERVICE_H