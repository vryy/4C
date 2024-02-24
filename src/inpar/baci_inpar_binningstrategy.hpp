/*-----------------------------------------------------------*/
/*! \file

\brief input parameter for binning strategy


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef BACI_INPAR_BINNINGSTRATEGY_HPP
#define BACI_INPAR_BINNINGSTRATEGY_HPP

#include "baci_config.hpp"

#include "baci_utils_parameter_list.hpp"

BACI_NAMESPACE_OPEN

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
BACI_NAMESPACE_CLOSE

#endif
