/*----------------------------------------------------------------------*/
/*! \file

\brief Object to store the beam to fluid volume meshtying output (visualization) parameters.

\level 3

*/



#ifndef BACI_FBI_BEAM_TO_FLUID_MESHTYING_VTK_OUTPUT_PARAMS_H
#define BACI_FBI_BEAM_TO_FLUID_MESHTYING_VTK_OUTPUT_PARAMS_H

#include "baci_beaminteraction_beam_to_solid_volume_meshtying_vtk_output_params.H"


namespace FBI
{
  class BeamToFluidMeshtyingVtkOutputParams
      : public BEAMINTERACTION::BeamToSolidVolumeMeshtyingVtkOutputParams
  {
   public:
    /**
     * \brief Constructor.
     */
    BeamToFluidMeshtyingVtkOutputParams();

    /**
     * \brief Destructor.
     */
    ~BeamToFluidMeshtyingVtkOutputParams() override = default;

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

#endif
