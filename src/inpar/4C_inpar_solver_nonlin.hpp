/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for nonlinear solvers


\level 1
*/

/*----------------------------------------------------------------------*/

#ifndef FOUR_C_INPAR_SOLVER_NONLIN_HPP
#define FOUR_C_INPAR_SOLVER_NONLIN_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"



FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace Inpar
{
  namespace NlnSol
  {
    /// set the nonlinear solver parameters
    void set_valid_parameters(Teuchos::ParameterList& list);

  }  // namespace NlnSol
}  // namespace Inpar

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
