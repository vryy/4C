/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_to_solid_volume_meshtying_evaluation_data.cpp

\brief data container holding all beam to solid volume meshtying
       evaluation data, that can not be stored pairwise.

\level 3

\maintainer Alexander Popp
*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam_to_solid_volume_meshtying_evaluation_data.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamToSolidVolumeMeshtyingEvaluationData::
    BeamToSolidVolumeMeshtyingEvaluationData()
    : isinit_(false),
      issetup_(false),
      BTSVOLMT_GaussPointProjectionTracker_ptr_(new std::map<int, std::vector<bool>>),
      BTSVOLMT_GaussPointEvaluateTracker_ptr_(new std::map<int, std::vector<bool>>)
{
  // Empty Constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingEvaluationData::Init()
{
  // empty for now

  issetup_ = false;

  isinit_ = true;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingEvaluationData::Setup()
{
  CheckInit();

  // empty for now

  issetup_ = true;
}
