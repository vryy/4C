/*-----------------------------------------------------------*/
/*! \file

\brief input parameter for rebalancing the discretization

\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_INPAR_REBALANCE_HPP
#define FOUR_C_INPAR_REBALANCE_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Inpar::Rebalance
{
  //! set the parameters for the geometric search strategy
  void set_valid_parameters(Teuchos::ParameterList& list);
}  // namespace Inpar::Rebalance

FOUR_C_NAMESPACE_CLOSE

#endif
