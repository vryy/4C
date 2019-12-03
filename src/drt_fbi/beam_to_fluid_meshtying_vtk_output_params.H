/*----------------------------------------------------------------------*/
/*! \file

\brief Object to store the beam to fluid volume meshtying output (visualization) parameters.

\level 3

\maintainer Nora Hagmeyer
*/



#ifndef BEAM_TO_FLUID_MESHTYING_VTK_OUTPUT_PARAMS_H
#define BEAM_TO_FLUID_MESHTYING_VTK_OUTPUT_PARAMS_H

#include "../drt_beaminteraction/beam_to_solid_volume_meshtying_vtk_output_params.H"


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
    virtual ~BeamToFluidMeshtyingVtkOutputParams(){};

    /**
     * \brief Setup member variables.
     */
    void Setup();
  };
}  // namespace FBI

#endif
