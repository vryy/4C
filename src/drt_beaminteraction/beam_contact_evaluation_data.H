/*----------------------------------------------------------------------------*/
/*!
\file beam_contact_evaluation_data.H

\brief data container holding pointers to all subcontainers that in turn hold
       all evaluation data specific to their problem type,
       that can not be stored pairwise.

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/


#ifndef BEAMCONTACTEVALUATIONDATA_H
#define BEAMCONTACTEVALUATIONDATA_H

#include <Teuchos_RCP.hpp>

namespace BEAMINTERACTION
{
  // forward declaration
  class BeamToSolidVolumeMeshtyingEvaluationData;

  /*!
   *  */
  class BeamContactEvaluationData
  {
   public:
    //! constructor
    BeamContactEvaluationData();

    //! destructor
    virtual ~BeamContactEvaluationData(){};

    //! builds a new BeamToSolidVolumeMeshtyingEvaluationData object
    void BuildBeamToSolidVolumeMeshtyingEvaluationData();

    //! access to the evaluation data object
    inline Teuchos::RCP<BEAMINTERACTION::BeamToSolidVolumeMeshtyingEvaluationData>
    BeamToSolidVolumeMeshtyingEvaluationData() const
    {
      return beam_to_solid_volume_meshtying_evaluation_data_;
    }


   private:
    //! pointer to the evaluation data class of beam-to-solid-volume contact
    Teuchos::RCP<BEAMINTERACTION::BeamToSolidVolumeMeshtyingEvaluationData>
        beam_to_solid_volume_meshtying_evaluation_data_;
  };

}  // namespace BEAMINTERACTION

#endif
