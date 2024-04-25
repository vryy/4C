/*----------------------------------------------------------------------*/
/*! \file

\brief Object to store the beam to fluid volume meshtying output (visualization) parameters.

\level 3

*/



#ifndef FOUR_C_FBI_BEAM_TO_FLUID_MESHTYING_OUTPUT_PARAMS_HPP
#define FOUR_C_FBI_BEAM_TO_FLUID_MESHTYING_OUTPUT_PARAMS_HPP

#include "4C_config.hpp"

#include "4C_beaminteraction_beam_to_solid_volume_meshtying_visualization_output_params.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FBI
{
  class BeamToFluidMeshtyingVtkOutputParams
      : public BEAMINTERACTION::BeamToSolidVolumeMeshtyingVisualizationOutputParams
  {
   public:
    /**
     * \brief Constructor.
     */
    BeamToFluidMeshtyingVtkOutputParams();

    /**
     * \brief Setup member variables.
     */
    void Setup();

    // returns flag enabling/disabling output of the constraint violation
    bool GetConstraintViolationOutputFlag() const { return constraint_violation_; };

   private:
    // bool enabling/disabling output of the constraint_enforcer violation
    bool constraint_violation_;
  };
}  // namespace FBI

FOUR_C_NAMESPACE_CLOSE

#endif
