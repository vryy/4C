/*----------------------------------------------------------------------------*/
/*! \file

\brief data container holding pointers to all subcontainers that in turn hold
       all input parameters specific to their problem type

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/


#ifndef BEAMCONTACTPARAMS_H
#define BEAMCONTACTPARAMS_H

#include "../drt_inpar/inpar_beamcontact.H"
//#include "beam_to_beam_contact_params.H"
//#include "beam_to_sphere_contact_params.H"
//#include "beam_to_solid_volume_meshtying_params.H"
//#include "beam_contact_runtime_vtk_output_params.H"


namespace BEAMINTERACTION
{
  class BeamToBeamContactParams;
  class BeamToSphereContactParams;
  class BeamToSolidVolumeMeshtyingParams;
  class BeamToSolidSurfaceMeshtyingParams;
  class BeamContactRuntimeVtkOutputParams;

  /*!
   *  */
  class BeamContactParams
  {
   public:
    //! constructor
    BeamContactParams();

    //! destructor
    virtual ~BeamContactParams(){};

    //! builds a new BeamToBeamContactParams object
    void BuildBeamToBeamContactParams();

    //! builds a new BeamToSphereContactParams object
    void BuildBeamToSphereContactParams();

    //! builds a new BeamToSolidVolumeMeshtyingParams object
    void BuildBeamToSolidVolumeMeshtyingParams();

    //! builds a new BeamToSolidSurfaceMeshtyingParams object
    void BuildBeamToSolidSurfaceMeshtyingParams();

    //! builds a new BeamContactRuntimeVtkParams object
    void BuildBeamContactRuntimeVtkOutputParams();


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

    inline Teuchos::RCP<BEAMINTERACTION::BeamContactRuntimeVtkOutputParams>
    BeamContactRuntimeVtkOutputParams() const
    {
      return beam_contact_runtime_vtk_output_params_;
    }


   private:
    //! pointer to the parameter class of beam-to-beam contact
    Teuchos::RCP<BEAMINTERACTION::BeamToBeamContactParams> beam_to_beam_contact_params_;

    //! pointer to the parameter class of beam-to-sphere contact
    Teuchos::RCP<BEAMINTERACTION::BeamToSphereContactParams> beam_to_sphere_contact_params_;

    //! pointer to the parameter class of beam-to-solid-volume contact
    Teuchos::RCP<BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams>
        beam_to_solid_volume_meshtying_params_;

    //! pointer to the parameter class of beam-to-solid-surface contact
    Teuchos::RCP<BEAMINTERACTION::BeamToSolidSurfaceMeshtyingParams>
        beam_to_solid_surface_meshtying_params_;

    //! pointer to the parameter class of beam contact vtk output
    Teuchos::RCP<BEAMINTERACTION::BeamContactRuntimeVtkOutputParams>
        beam_contact_runtime_vtk_output_params_;
  };

}  // namespace BEAMINTERACTION

#endif
