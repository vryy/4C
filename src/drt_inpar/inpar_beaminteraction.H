/*-----------------------------------------------------------*/
/*!
\file inpar_beaminteraction.H

\brief input parameter for beaminteraction

\maintainer Jonas Eichinger, Maximilian Grill

\level 2

*/
/*-----------------------------------------------------------*/
#ifndef INPAR_BEAMINTERACTION_H
#define INPAR_BEAMINTERACTION_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters used in statistical mechanics
  namespace BEAMINTERACTION
  {

  /// Type of the used submodel for beaminteraction
  enum SubModelType
  {
    submodel_crosslinking,            ///< evaluate the structural model
    submodel_beamcontact,              ///< evaluate the contact model
    submodel_potential,            ///< evaluate the meshtying model
    submodel_vague                 ///< undefined model type
  };

    /// set the beaminteraction parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }

}

/*----------------------------------------------------------------------*/
#endif
