
/*-----------------------------------------------------------*/
/*!
\file inpar_crosslinking.H

\brief input parameter for statistical mechanic problem

\maintainer Jonas Eichinger, Maximilian Grill

\level 2

*/
/*-----------------------------------------------------------*/
#ifndef INPAR_CROSSLINKING_H_
#define INPAR_CROSSLINKING_H_

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters used in statistical mechanics
  namespace CROSSLINKING
  {

    /// set the crosslinking parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }

}

/*----------------------------------------------------------------------*/
#endif


