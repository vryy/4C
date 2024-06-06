/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameters for linear solvers

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_INPAR_SOLVER_HPP
#define FOUR_C_INPAR_SOLVER_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Inpar::SOLVER
{
  //! set the parameters for the linear solver
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

}  // namespace Inpar::SOLVER

FOUR_C_NAMESPACE_CLOSE

#endif
