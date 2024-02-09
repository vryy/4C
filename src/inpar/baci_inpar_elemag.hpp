/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for electromagnetic simulations

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_INPAR_ELEMAG_HPP
#define BACI_INPAR_ELEMAG_HPP

#include "baci_config.hpp"

#include "baci_inpar_parameterlist_utils.hpp"

BACI_NAMESPACE_OPEN

// forward declaration
namespace INPUT
{
  class ConditionDefinition;
}

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace ELEMAG
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
    void SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist);

  }  // namespace ELEMAG
}  // namespace INPAR


BACI_NAMESPACE_CLOSE

#endif  // INPAR_ELEMAG_H
