/*----------------------------------------------------------------------*/
/*! \file

\brief Common functionalities for materials that base on a multiplicative split of the deformation
gradient.

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_MULTIPLICATIVE_SPLIT_DEFGRAD_ELASTHYPER_SERVICE_HPP
#define FOUR_C_MAT_MULTIPLICATIVE_SPLIT_DEFGRAD_ELASTHYPER_SERVICE_HPP
#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_mat_elasthyper_service.hpp"
#include "4C_mat_service.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  inline void EvaluateCe(const Core::LinAlg::Matrix<3, 3>& F,
      const Core::LinAlg::Matrix<3, 3>& iFin, Core::LinAlg::Matrix<3, 3>& Ce)
  {
    static Core::LinAlg::Matrix<3, 3> FiFin(false);
    FiFin.MultiplyNN(F, iFin);
    Ce.MultiplyTN(FiFin, FiFin);
  }

  inline void EvaluateiCinCiCin(const Core::LinAlg::Matrix<3, 3>& C,
      const Core::LinAlg::Matrix<3, 3>& iCin, Core::LinAlg::Matrix<3, 3>& iCinCiCin)
  {
    static Core::LinAlg::Matrix<3, 3> CiCin(false);
    CiCin.MultiplyNN(C, iCin);
    iCinCiCin.MultiplyNN(iCin, CiCin);
  }

  inline void ElastHyperEvaluateElasticPart(const Core::LinAlg::Matrix<3, 3>& F,
      const Core::LinAlg::Matrix<3, 3>& iFin, Core::LinAlg::Matrix<6, 1>& S_stress,
      Core::LinAlg::Matrix<6, 6>& cmat,
      const std::vector<Teuchos::RCP<Mat::Elastic::Summand>>& potsum,
      Mat::SummandProperties summandProperties, const int gp, const int eleGID)
  {
    if (summandProperties.anisomod or summandProperties.anisoprinc)
    {
      FOUR_C_THROW(
          "An additional inelastic part is not yet implemented for anisotropic materials.");
    }

    S_stress.Clear();
    cmat.Clear();

    // Variables needed for the computation of the stress resultants
    static Core::LinAlg::Matrix<3, 3> C(true);
    static Core::LinAlg::Matrix<3, 3> Ce(true);
    static Core::LinAlg::Matrix<3, 3> iC(true);
    static Core::LinAlg::Matrix<3, 3> iCin(true);
    static Core::LinAlg::Matrix<3, 3> iCinCiCin(true);

    static Core::LinAlg::Matrix<6, 1> iCinv(true);
    static Core::LinAlg::Matrix<6, 1> iCinCiCinv(true);
    static Core::LinAlg::Matrix<6, 1> iCv(true);
    static Core::LinAlg::Matrix<3, 1> principleInvariantsCe(true);

    // Compute right Cauchy-Green tensor C=F^TF
    C.MultiplyTN(F, F);

    // Compute inverse right Cauchy-Green tensor C^-1
    iC.Invert(C);

    // Compute inverse inelastic right Cauchy-Green Tensor
    iCin.MultiplyNT(iFin, iFin);

    // Compute iCin * C * iCin
    Mat::EvaluateiCinCiCin(C, iCin, iCinCiCin);

    // Compute Ce
    Mat::EvaluateCe(F, iFin, Ce);

    // Compute principal invariants
    Mat::invariants_principal(principleInvariantsCe, Ce);

    Core::LinAlg::Matrix<3, 1> dPIe(true);
    Core::LinAlg::Matrix<6, 1> ddPIIe(true);

    Mat::ElastHyperEvaluateInvariantDerivatives(
        principleInvariantsCe, dPIe, ddPIIe, potsum, summandProperties, gp, eleGID);

    // 2nd Piola Kirchhoff stress factors (according to Holzapfel-Nonlinear Solid Mechanics p. 216)
    static Core::LinAlg::Matrix<3, 1> gamma(true);
    // constitutive tensor factors (according to Holzapfel-Nonlinear Solid Mechanics p. 261)
    static Core::LinAlg::Matrix<8, 1> delta(true);

    Mat::CalculateGammaDelta(gamma, delta, principleInvariantsCe, dPIe, ddPIIe);

    // Convert necessary tensors to stress-like Voigt-Notation
    Core::LinAlg::Voigt::Stresses::matrix_to_vector(iCin, iCinv);
    Core::LinAlg::Voigt::Stresses::matrix_to_vector(iCinCiCin, iCinCiCinv);
    Core::LinAlg::Voigt::Stresses::matrix_to_vector(iC, iCv);

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
    Mat::add_holzapfel_product(cmat, iCv, delta(6));
    Mat::add_holzapfel_product(cmat, iCinv, delta(7));
  }

}  // namespace Mat
FOUR_C_NAMESPACE_CLOSE

#endif