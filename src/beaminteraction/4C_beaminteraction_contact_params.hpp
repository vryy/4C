/*----------------------------------------------------------------------------*/
/*! \file

\brief data container holding pointers to all subcontainers that in turn hold
       all input parameters specific to their problem type

\level 3

*/
/*----------------------------------------------------------------------------*/


#ifndef FOUR_C_BEAMINTERACTION_CONTACT_PARAMS_HPP
#define FOUR_C_BEAMINTERACTION_CONTACT_PARAMS_HPP

#include "4C_config.hpp"

#include "4C_inpar_beamcontact.hpp"

FOUR_C_NAMESPACE_OPEN


namespace BEAMINTERACTION
{
  class BeamToBeamContactParams;
  class BeamToSphereContactParams;
  class BeamToSolidVolumeMeshtyingParams;
  class BeamToSolidSurfaceMeshtyingParams;
  class BeamToSolidSurfaceContactParams;
  class BeamContactRuntimeVisualizationOutputParams;

  /*!
   *  */
  class BeamContactParams
  {
   public:
    //! constructor
    BeamContactParams();

    //! destructor
    virtual ~BeamContactParams() = default;

    //! builds a new beam_to_beam_contact_params object
    void build_beam_to_beam_contact_params();

    //! builds a new beam_to_sphere_contact_params object
    void build_beam_to_sphere_contact_params();

    //! builds a new beam_to_solid_volume_meshtying_params object
    void build_beam_to_solid_volume_meshtying_params();

    //! builds a new beam_to_solid_surface_meshtying_params object
    void build_beam_to_solid_surface_meshtying_params();

    //! builds a new beam_to_solid_surface_contact_params object
    void build_beam_to_solid_surface_contact_params();

    //! builds a new BeamContactRuntimeOutputParams object
    void build_beam_contact_runtime_output_params(double restart_time);


    inline Teuchos::RCP<BEAMINTERACTION::BeamToBeamContactParams> beam_to_beam_contact_params()
        const
    {
      return beam_to_beam_contact_params_;
    }

    inline Teuchos::RCP<BEAMINTERACTION::BeamToSphereContactParams> beam_to_sphere_contact_params()
        const
    {
      return beam_to_sphere_contact_params_;
    }

    inline Teuchos::RCP<BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams>
    beam_to_solid_volume_meshtying_params() const
    {
      return beam_to_solid_volume_meshtying_params_;
    }

    inline Teuchos::RCP<BEAMINTERACTION::BeamToSolidSurfaceMeshtyingParams>
    beam_to_solid_surface_meshtying_params() const
    {
      return beam_to_solid_surface_meshtying_params_;
    }

    inline Teuchos::RCP<BEAMINTERACTION::BeamToSolidSurfaceContactParams>
    beam_to_solid_surface_contact_params() const
    {
      return beam_to_solid_surface_contact_params_;
    }

    inline Teuchos::RCP<BEAMINTERACTION::BeamContactRuntimeVisualizationOutputParams>
    beam_contact_runtime_visualization_output_params() const
    {
      return beam_contact_runtime_output_params_;
    }


   private:
    //! pointer to the parameter class of beam-to-beam contact
    Teuchos::RCP<BEAMINTERACTION::BeamToBeamContactParams> beam_to_beam_contact_params_;

    //! pointer to the parameter class of beam-to-sphere contact
    Teuchos::RCP<BEAMINTERACTION::BeamToSphereContactParams> beam_to_sphere_contact_params_;

    //! pointer to the parameter class of beam-to-solid-volume contact
    Teuchos::RCP<BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams>
        beam_to_solid_volume_meshtying_params_;

    //! pointer to the parameter class of beam-to-solid-surface mesh tying
    Teuchos::RCP<BEAMINTERACTION::BeamToSolidSurfaceMeshtyingParams>
        beam_to_solid_surface_meshtying_params_;

    //! pointer to the parameter class of beam-to-solid-surface contact
    Teuchos::RCP<BEAMINTERACTION::BeamToSolidSurfaceContactParams>
        beam_to_solid_surface_contact_params_;

    //! pointer to the parameter class of beam contact visualization output
    Teuchos::RCP<BEAMINTERACTION::BeamContactRuntimeVisualizationOutputParams>
        beam_contact_runtime_output_params_;
  };

}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
