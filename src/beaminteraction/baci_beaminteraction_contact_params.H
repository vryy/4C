/*----------------------------------------------------------------------------*/
/*! \file

\brief data container holding pointers to all subcontainers that in turn hold
       all input parameters specific to their problem type

\level 3

*/
/*----------------------------------------------------------------------------*/


#ifndef BACI_BEAMINTERACTION_CONTACT_PARAMS_H
#define BACI_BEAMINTERACTION_CONTACT_PARAMS_H

#include "baci_config.H"

#include "baci_inpar_beamcontact.H"

BACI_NAMESPACE_OPEN


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

    //! builds a new BeamToBeamContactParams object
    void BuildBeamToBeamContactParams();

    //! builds a new BeamToSphereContactParams object
    void BuildBeamToSphereContactParams();

    //! builds a new BeamToSolidVolumeMeshtyingParams object
    void BuildBeamToSolidVolumeMeshtyingParams();

    //! builds a new BeamToSolidSurfaceMeshtyingParams object
    void BuildBeamToSolidSurfaceMeshtyingParams();

    //! builds a new BeamToSolidSurfaceContactParams object
    void BuildBeamToSolidSurfaceContactParams();

    //! builds a new BeamContactRuntimeOutputParams object
    void BuildBeamContactRuntimeOutputParams(double restart_time);


    inline Teuchos::RCP<BEAMINTERACTION::BeamToBeamContactParams> BeamToBeamContactParams() const
    {
      return beam_to_beam_contact_params_;
    }

    inline Teuchos::RCP<BEAMINTERACTION::BeamToSphereContactParams> BeamToSphereContactParams()
        const
    {
      return beam_to_sphere_contact_params_;
    }

    inline Teuchos::RCP<BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams>
    BeamToSolidVolumeMeshtyingParams() const
    {
      return beam_to_solid_volume_meshtying_params_;
    }

    inline Teuchos::RCP<BEAMINTERACTION::BeamToSolidSurfaceMeshtyingParams>
    BeamToSolidSurfaceMeshtyingParams() const
    {
      return beam_to_solid_surface_meshtying_params_;
    }

    inline Teuchos::RCP<BEAMINTERACTION::BeamToSolidSurfaceContactParams>
    BeamToSolidSurfaceContactParams() const
    {
      return beam_to_solid_surface_contact_params_;
    }

    inline Teuchos::RCP<BEAMINTERACTION::BeamContactRuntimeVisualizationOutputParams>
    BeamContactRuntimeVisualizationOutputParams() const
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

BACI_NAMESPACE_CLOSE

#endif
