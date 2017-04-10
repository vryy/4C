/*-----------------------------------------------------------*/
/*!
\file inpar_browniandyn.H

\brief input parameter for Brownian dynamics simulation

\maintainer Jonas Eichinger, Maximilian Grill

\level 2

*/
/*-----------------------------------------------------------*/
#ifndef INPAR_BROWNIANDYN_H
#define INPAR_BROWNIANDYN_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters used in statistical mechanics
  namespace BROWNIANDYN
  {
    /// the way how damping coefficient values for beams are specified
    enum BeamDampingCoefficientSpecificationType
    {
      cylinder_geometry_approx,
      input_file,
      vague
    };

    /// set the brownian dynamic parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }

}

/*----------------------------------------------------------------------*/
#endif
