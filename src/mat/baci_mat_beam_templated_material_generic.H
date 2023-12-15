/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief templated interface for constitutive relations for beam cross-section resultants
\level 1
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BACI_MAT_BEAM_TEMPLATED_MATERIAL_GENERIC_H
#define BACI_MAT_BEAM_TEMPLATED_MATERIAL_GENERIC_H

#include "baci_mat_beam_material_generic.H"

BACI_NAMESPACE_OPEN


// forward declaration
namespace DRT
{
  class ParObject;
}

namespace MAT
{
  // forward declaration
  namespace PAR
  {
    class BeamElastHyperMaterialParameterGeneric;
  }

  /*---------------------------------------------------------------------------------------------*/
  /// constitutive relations for beam cross-section resultants (hyperelastic stored energy function)
  template <typename T>
  class BeamMaterialTemplated : public BeamMaterial
  {
   public:
    /*
     * \brief Compute axial stress contributions
     *
     *\param[out] stressM axial stress
     *
     *\param[in] CM constitutive matrix
     *
     *\param[in] Cur curvature
     */
    virtual void EvaluateMomentContributionsToStress(CORE::LINALG::Matrix<3, 1, T>& stressM,
        const CORE::LINALG::Matrix<3, 3, T>& CM, const CORE::LINALG::Matrix<3, 1, T>& Cur,
        const unsigned int gp) = 0;

    /*
     * \brief Compute axial stress contributions
     *
     *\param[out] stressN axial stress
     *
     *\param[in] CN constitutive matrix
     *
     *\param[in] Gamma triad
     */

    virtual void EvaluateForceContributionsToStress(CORE::LINALG::Matrix<3, 1, T>& stressN,
        const CORE::LINALG::Matrix<3, 3, T>& CN, const CORE::LINALG::Matrix<3, 1, T>& Gamma,
        const unsigned int gp) = 0;

    /*
     * \brief Update material-dependent variables
     */
    virtual void ComputeConstitutiveParameter(
        CORE::LINALG::Matrix<3, 3, T>& C_N, CORE::LINALG::Matrix<3, 3, T>& C_M) = 0;

    /** \brief get constitutive matrix relating stress force resultants and translational strain
     *         measures, expressed w.r.t. material frame
     *
     */
    virtual void GetConstitutiveMatrixOfForcesMaterialFrame(
        CORE::LINALG::Matrix<3, 3, T>& C_N) const = 0;

    /** \brief get constitutive matrix relating stress moment resultants and rotational strain
     *         measures, expressed w.r.t. material frame
     *
     */
    virtual void GetConstitutiveMatrixOfMomentsMaterialFrame(
        CORE::LINALG::Matrix<3, 3, T>& C_M) const = 0;

    /** \brief get linearization of the constitutive law relating stress moment resultants and
     * rotational strain measures, expressed w.r.t. material frame
     *
     *\param[in] C_M constitutive matrix
     *\param[out] stiffness_matrix
     */
    virtual void GetStiffnessMatrixOfMoments(CORE::LINALG::Matrix<3, 3, T>& stiffness_matrix,
        const CORE::LINALG::Matrix<3, 3, T>& C_M, const int gp) = 0;

    /** \brief get linearization of the constitutive law relating stress force resultants and
     * translational strain measures, expressed w.r.t. material frame
     *
     *\param[in] C_N constitutive matrix
     *\param[out] stiffness_matrix
     */
    virtual void GetStiffnessMatrixOfForces(CORE::LINALG::Matrix<3, 3, T>& stiffness_matrix,
        const CORE::LINALG::Matrix<3, 3, T>& C_N, const int gp) = 0;
  };
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif
