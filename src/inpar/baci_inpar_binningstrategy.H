/*-----------------------------------------------------------*/
/*! \file

\brief input parameter for binning strategy


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef BACI_INPAR_BINNINGSTRATEGY_H
#define BACI_INPAR_BINNINGSTRATEGY_H

#include "baci_inpar_parameterlist_utils.H"

namespace INPAR
{
  namespace BINSTRATEGY
  {
    /*!
     * \ brief write either no, row or column bins for visualization
     */
    enum writebins
    {
      none,
      rows,
      cols
    };

    /// set the binning strategy parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace BINSTRATEGY

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif
