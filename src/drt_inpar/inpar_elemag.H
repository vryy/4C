/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for electromagnetic simulations

\level 3

\maintainer Luca Berardocco
*/
/*----------------------------------------------------------------------*/

#ifndef INPAR_ELEMAG_H
#define INPAR_ELEMAG_H

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

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
      elemag_implicit_euler,
      /// BDF2
      elemag_bdf,
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
      initfield_field_by_function
    };

    /// Define all valid parameters for electromagnetic problem.
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// Set specific electromagnetic conditions.
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

  }  // namespace ELEMAG
}  // namespace INPAR


#endif  // #ifndef INPAR_ELEMAG_H
