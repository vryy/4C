/*----------------------------------------------------------------------*/
/*! \file

\brief Contains the declaration of service functions for hyperelastic membrane materials

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef MEMBRANE_ELASTHYPER_SERVICE_H
#define MEMBRANE_ELASTHYPER_SERVICE_H

#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_matelast/elast_summand.H"
namespace MAT
{
  // Forward declaration
  class SummandProperties;

  /*!
   * @brief hyperelastic stress response plus elasticity tensor for membranes assuming
   * incompressibility and plane stress
   *
   * @param cauchygreen (in) : right Cauchy-Green tensor
   * @param params (in) : Container for additional information
   * @param Q_trafo (in) : Trafo from local membrane orthonormal coordinates to global coordinates
   * @param stress (out) : 2nd Piola-Kirchhoff stresses in stress-like voigt notation
   * @param cmat (out) : Constitutive matrix
   * @param gp (in) : Gauss point
   * @param eleGID (in) : Global element id
   * @param potsum (in) : Vector of summands of the strain energy function
   * @param properties (in) : Properties of the summands of the strain energy function
   */
  void MembraneElastHyperEvaluateIsotropicStressCmat(const LINALG::Matrix<3, 3>& cauchygreen,
      Teuchos::ParameterList& params, const LINALG::Matrix<3, 3>& Q_trafo,
      LINALG::Matrix<3, 1>& stress, LINALG::Matrix<3, 3>& cmat, int gp, int eleGID,
      const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum,
      const SummandProperties& properties);


  /*!
   * @brief calculates the kinematic quantities and tensors used afterwards
   *
   * @param cauchygreen (in) : right Cauchy-Green tensor
   * @param id2 (out) : cartesian identity 2-tensor I_{AB}
   * @param id4sharp (out) : cartesian identity 4-tensor
   * @param rcg (out) : right Cauchy-Green in stress-like 3-Voigt notation
   * @param rcg33 (out) : principal stretch in thickness direction
   * @param icg (out) : inverse right Cauchy-Green in stress-like 3-Voigt notation
   */
  void MembraneElastHyperEvaluateKinQuant(const LINALG::Matrix<3, 3>& cauchygreen,
      LINALG::Matrix<3, 1>& id2, LINALG::Matrix<3, 3>& id4sharp, LINALG::Matrix<3, 1>& rcg,
      double& rcg33, LINALG::Matrix<3, 1>& icg);

  /*!
   * @brief Computes the isotropic stress response and the linearization
   *
   * @param stress_iso (out) : Isotropic stress response in local membrane coordinate system (plane
   * stress)
   * @param cmat_iso (out) : Linearization of stress response
   * @param id2 (in) : 2D identity tensor
   * @param id4sharp (in) : Fourth order identity tensor
   * @param rcg (in) : Right Cauchy-Green tensor in local membrane coordinate system (strain like)
   * @param rcg33 (in) : Right Cauchy-Green tensor in membrane thickness direction
   * @param icg (in) : Inverse Right Cauchy-Green tensor in local membrane coordinate system (strain
   * like)
   * @param gp (in) : Guass point
   * @param eleGID (in) : global element id
   * @param potsum (in) : Vector of summands of the strain energy function
   * @param properties (in) : Properties of the summands of the strain energy function
   */
  void MembraneElastHyperEvaluateIsotropicStressCmat(LINALG::Matrix<3, 1>& stress_iso,
      LINALG::Matrix<3, 3>& cmat_iso, const LINALG::Matrix<3, 1>& id2,
      const LINALG::Matrix<3, 3>& id4sharp, const LINALG::Matrix<3, 1>& rcg, const double& rcg33,
      const LINALG::Matrix<3, 1>& icg, int gp, int eleGID,
      const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum,
      const SummandProperties& properties);

  /*!
   * @brief calculate principal invariants
   *
   * @param prinv (out) : principal invariants
   * @param rcg (in) : right cauchy-green in stress-like 3-Voigt notation
   * @param rcg33 (in) :principal stretch in thickness direction
   */
  void MembraneElastHyperInvariantsPrincipal(
      LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<3, 1>& rcg, const double& rcg33);

}  // namespace MAT

#endif  // MEMBRANE_ELASTHYPER_SERVICE_H