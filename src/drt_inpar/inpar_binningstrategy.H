/*-----------------------------------------------------------*/
/*!
\file inpar_binningstrategy.H

\brief input parameter for binning strategy

\maintainer Jonas Eichinger

\level 2

*/
/*-----------------------------------------------------------*/

#ifndef INPAR_BINNINGSTRATEGY_H
#define INPAR_BINNINGSTRATEGY_H

#include "inpar_parameterlist_utils.H"

namespace INPAR
{

/// Input parameters for Meshfree
namespace BINSTRATEGY
{

  //! type of bin specification
  enum xaabbspectype {
    input,
    dynamic
  };

  //! type of bin specification
  enum binspectype {
    cutoff,
    binsperdir,
    largestele
  };

  //! type of gap function/projection
  enum writebins {
    none,
    rows,
    cols
  };

  /// set the meshfree parameters
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

}

}

/*----------------------------------------------------------------------*/
#endif
