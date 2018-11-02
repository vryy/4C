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
  namespace BINSTRATEGY
  {
    //! type of gap function/projection
    enum writebins
    {
      none,
      rows,
      cols
    };

    /// set the meshfree parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace BINSTRATEGY

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif
