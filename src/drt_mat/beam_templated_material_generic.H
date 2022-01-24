/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief templated interface for constitutive relations for beam cross-section resultants
\level 1
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef MAT_BEAM_TEMPLATED_H_
#define MAT_BEAM_TEMPLATED_H_

#include "beam_material_generic.H"


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
     *\param[in] K curvature
     */
    virtual void EvaluateMomentContributionsToStress(LINALG::Matrix<3, 1, T>& stressM,
        const LINALG::Matrix<3, 3, T>& CM, const LINALG::Matrix<3, 1, T>& K) const = 0;

    /*
     * \brief Compute axial stress contributions
     *
     *\param[out] stressN axial stress
     *
     *\param[in] CN constitutive matrix
     *
     *\param[in] Gamma triad
     */

    virtual void EvaluateForceContributionsToStress(LINALG::Matrix<3, 1, T>& stressN,
        const LINALG::Matrix<3, 3, T>& CN, const LINALG::Matrix<3, 1, T>& Gamma) const = 0;

    /*
     * \brief Update material-dependent variables
     */
    virtual void ComputeConstitutiveParameter(
        LINALG::Matrix<3, 3, T>& C_N, LINALG::Matrix<3, 3, T>& C_M) const = 0;

    /** \brief get constitutive matrix relating stress force resultants and translational strain
     *         measures, expressed w.r.t. material frame
     *
     */
    virtual void GetConstitutiveMatrixOfForcesMaterialFrame(LINALG::Matrix<3, 3, T>& C_N) const = 0;

    /** \brief get constitutive matrix relating stress moment resultants and rotational strain
     *         measures, expressed w.r.t. material frame
     *
     */
    virtual void GetConstitutiveMatrixOfMomentsMaterialFrame(
        LINALG::Matrix<3, 3, T>& C_M) const = 0;

    virtual void GetStiffnessMatrixOfMoments(
        LINALG::Matrix<3, 3, T>& stiffness_matrix, const LINALG::Matrix<3, 3, T>& C_M) const = 0;

    virtual void GetStiffnessMatrixOfForces(
        LINALG::Matrix<3, 3, T>& stiffness_matrix, const LINALG::Matrix<3, 3, T>& C_N) const = 0;
  };
}  // namespace MAT

#endif
