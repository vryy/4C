/*-----------------------------------------------------------*/
/*! \file

\brief Input parameter for beam-to-solid interaction.


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef INPAR_BEAM_TO_SOLID_H
#define INPAR_BEAM_TO_SOLID_H


#include "inpar_parameterlist_utils.H"

#include "../drt_lib/drt_condition.H"


// Forward declaration.
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT
namespace INPAR
{
  namespace BEAMINTERACTION
  {
    enum class BeamInteractionConditions;
  }
}  // namespace INPAR


namespace INPAR
{
  namespace BEAMTOSOLID
  {
    /**
     * \brief Get the condition strings for the beam-to-solid interactions.
     */
    void BeamToSolidInteractionGetString(
        const INPAR::BEAMINTERACTION::BeamInteractionConditions& interaction,
        std::array<std::string, 2>& condition_names);

    /**
     * \brief Constraint enforcement for beam to solid interactions.
     */
    enum class BeamToSolidConstraintEnforcement
    {
      //! Default value.
      none,
      //! Penalty method.
      penalty
    };

    /**
     * \brief Contact discretization for beam to solid interactions.
     */
    enum class BeamToSolidContactDiscretization
    {
      //! Default value.
      none,
      //! Gauss point to segment approach.
      gauss_point_to_segment,
      //! Gauss point on cross section method.
      gauss_point_cross_section,
      //! Mortar type coupling.
      mortar
    };

    /**
     * \brief Contact discretization for beam to solid interactions.
     */
    enum class BeamToSolidRotationCoupling
    {
      //! Default value.
      none,
      //! Polar decomposition of the deformation gradient.
      polar_decomposition_2d,
      //! Deformed y basis vector.
      deformation_gradient_y_2d,
      //! Deformed z basis vector.
      deformation_gradient_z_2d,
      //! Average of first and second basis vector of Deformation gradient.
      deformation_gradient_average_2d,
      //! Fix the beam crosssection triad to the solid triad.
      fix_triad_2d,
      //! 3D case where the triad construction starts with the local e1 vector.
      deformation_gradient_3d_local_1,
      //! 3D case where the triad construction starts with the local e2 vector.
      deformation_gradient_3d_local_2,
      //! 3D case where the triad construction starts with the local e3 vector.
      deformation_gradient_3d_local_3
    };

    /**
     * \brief Shape function for the mortar Lagrange-multiplicators
     */
    enum class BeamToSolidMortarShapefunctions
    {
      //! Default value.
      none,
      //! Linear.
      line2,
      //! Quadratic.
      line3,
      //! Cubic.
      line4
    };

    /**
     * \brief How the coupling between a beam and a surface should be evaluated.
     */
    enum class BeamToSolidSurfaceCoupling
    {
      //! Default value.
      none,
      //! Couple the reference configurations of the beam and solid, but force the normal distance
      //! to 0. This will lead to angular momentum conservation, but initial strains.
      reference_configuration_forced_to_zero,
      //! Same as above, but evaluated with FAD types.
      reference_configuration_forced_to_zero_fad,
      //! Couple the displacements, if no mesh initialization is performed, this will lead to
      //! non-conserving angular momentum.
      displacement,
      //! Same as above, but evaluated with FAD types.
      displacement_fad,
      //! Consistent coupling defined via the normal vector. This requires derivatives of the
      //! averaged normals.
      consistent_fad
    };

    /**
    \brief Set the beam-to-solid parameters
    */
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /**
     * \brief Set beam beam-to-solid specific conditions.
     */
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);
  }  // namespace BEAMTOSOLID
}  // namespace INPAR

#endif
