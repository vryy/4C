/*-----------------------------------------------------------*/
/*! \file

\brief input parameter for binning strategy


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_INPAR_BINNINGSTRATEGY_HPP
#define FOUR_C_INPAR_BINNINGSTRATEGY_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Inpar
{
  namespace BINSTRATEGY
  {

    /// set the binning strategy parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace BINSTRATEGY

}  // namespace Inpar

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
