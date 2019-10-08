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

#include "../drt_fem_general/drt_utils_integration.H"
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
  namespace BEAMTOSOLID
  {
    /**
     * \brief Types of beam-to-solid interaction
     */
    enum class BeamToSolidInteraction
    {
      //! Default value.
      none,
      //! Beam-to-solid volume meshtying.
      beam_to_solid_volume_meshtying,
      //! Beam-to-solid surface meshtying.
      beam_to_solid_surface_meshtying
    };

    /**
     * \brief Get all available beam-to-solid interactions, excluding the default value.
     */
    void BeamToSolidInteractionGetAll(std::vector<BeamToSolidInteraction>& interactions);

    /**
     * \brief Get the condition strings for the beam-to-solid interactions.
     */
    void BeamToSolidInteractionGetString(
        const BeamToSolidInteraction& interaction, std::array<std::string, 2>& condition_names);

    /**
     * \brief Constraint enforcement for beam to solid volume mesh tying.
     */
    enum class BeamToSolidVolumeConstraintEnforcement
    {
      //! Default value.
      none,
      //! Penalty method.
      penalty
    };

    /**
     * \brief Contact discretization for beam to solid volume mesh tying.
     */
    enum class BeamToSolidVolumeContactDiscretization
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
    enum class BeamToSolidVolumeMortarShapefunctions
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
     * \brief Map number of gauss points to 1D gauss rule
     */
    DRT::UTILS::GaussRule1D IntToGaussRule1D(const int n_gauss_points);

    /**
     * \brief Set the beam-to-solid parameters
     */
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /**
     * \brief Set beam beam-to-solid specific conditions.
     */
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);
  }  // namespace BEAMTOSOLID
}  // namespace INPAR

#endif
