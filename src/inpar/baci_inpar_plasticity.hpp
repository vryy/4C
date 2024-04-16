/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for plasticity
\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_INPAR_PLASTICITY_HPP
#define FOUR_C_INPAR_PLASTICITY_HPP

#include "baci_config.hpp"

#include "baci_utils_parameter_list.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace PLASTICITY
  {
    /// set the plasticity parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace PLASTICITY

}  // namespace INPAR

/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
