/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for nonlinear solvers


\level 1
*/

/*----------------------------------------------------------------------*/

#ifndef BACI_INPAR_SOLVER_NONLIN_H
#define BACI_INPAR_SOLVER_NONLIN_H

#include "baci_config.H"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCPDecl.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace NLNSOL
  {
    /// set the nonlinear solver parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace NLNSOL
}  // namespace INPAR

/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif  // INPAR_SOLVER_NONLIN_H
