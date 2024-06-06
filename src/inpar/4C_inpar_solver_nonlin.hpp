/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for nonlinear solvers


\level 1
*/

/*----------------------------------------------------------------------*/

#ifndef FOUR_C_INPAR_SOLVER_NONLIN_HPP
#define FOUR_C_INPAR_SOLVER_NONLIN_HPP

#include "4C_config.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCPDecl.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace Inpar
{
  namespace NlnSol
  {
    /// set the nonlinear solver parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace NlnSol
}  // namespace Inpar

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
