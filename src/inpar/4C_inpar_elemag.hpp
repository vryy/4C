/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for electromagnetic simulations

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_INPAR_ELEMAG_HPP
#define FOUR_C_INPAR_ELEMAG_HPP

#include "4C_config.hpp"

#include "4C_discretization_condition_definition.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration
/*----------------------------------------------------------------------*/
namespace Inpar
{
  namespace EleMag
  {
    /// Type of time integrator
    enum DynamicType
    {
      /// one-step-theta time integration
      elemag_ost,
      /// explicit euler method
      elemag_explicit_euler,
      /// implicit euler method
      elemag_bdf1,
      /// BDF2
      elemag_bdf2,
      /// BDF4
      elemag_bdf4,
      /// Generalized-Alpha method
      elemag_genAlpha,
      /// runge-kutta method
      elemag_rk,
      /// crank-nicolson method
      elemag_cn
    };

    /// Initial field for electromagnetic problems.
    enum InitialField
    {
      /// Initialize a zero field on all the components
      initfield_zero_field,
      /// Initialize the components as specified by the function
      initfield_field_by_function,
      /// Initialize the electric field with a CG scatra solution
      initfield_scatra,
      /// Initialize the electric field with a HDG scatra solution
      initfield_scatra_hdg
    };

    /// Define all valid parameters for electromagnetic problem.
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// Set specific electromagnetic conditions.
    void SetValidConditions(
        std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist);

  }  // namespace EleMag
}  // namespace Inpar


FOUR_C_NAMESPACE_CLOSE

#endif
