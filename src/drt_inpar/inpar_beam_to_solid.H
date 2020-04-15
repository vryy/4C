/*-----------------------------------------------------------*/
/*! \file

\brief Input parameter for beam-to-solid interaction.

\maintainer Ivo Steinbrecher

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
      configurations_forced_to_zero,
      //! Couple the displacements, if no mesh initialization is performed, this will lead to
      //! non-conserving angular momentum.
      displacements
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
