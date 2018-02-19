/*----------------------------------------------------------------------------*/
/*!
\file beam_contact_evaluation_data.cpp

\brief data container holding pointers to all subcontainers that in turn hold
       all evaluation data specific to their problem type,
       that can not be stored pairwise.

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/


#include "beam_contact_evaluation_data.H"

#include "beam_to_solid_volume_meshtying_evaluation_data.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamContactEvaluationData::BeamContactEvaluationData():
beam_to_solid_volume_meshtying_evaluation_data_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactEvaluationData::
    BuildBeamToSolidVolumeMeshtyingEvaluationData()
{
  beam_to_solid_volume_meshtying_evaluation_data_ = Teuchos::rcp( new BEAMINTERACTION::
      BeamToSolidVolumeMeshtyingEvaluationData() );

  beam_to_solid_volume_meshtying_evaluation_data_->Init();
  beam_to_solid_volume_meshtying_evaluation_data_->Setup();
}
