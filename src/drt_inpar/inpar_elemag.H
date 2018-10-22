/*----------------------------------------------------------------------*/
/*!
\file inpar_elemag.H

\brief Input parameters for electromagnetic simulations

<pre>
\level 3

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            089 - 289-15244
</pre>
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
  // Input parameters for electromagnetic simulations
  namespace ELEMAG
  {
    // Type of time integrator
    enum DynamicType
    {
      elemag_ost,             // one-step-theta time integration
      elemag_explicit_euler,  // explicit euler method
      elemag_implicit_euler,  // implicit euler method
      elemag_rk,              // runge-kutta method
      elemag_cn               // crank-nicolson method
    };

    /// initial field for electromagnetic problems
    enum InitialField
    {
      initfield_zero_field,
      initfield_field_by_function
    };

    /// define all valid parameters for electromagnetic problem
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific electromagnetic conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

  }  // namespace ELEMAG
}  // namespace INPAR


#endif  // #ifndef INPAR_ELEMAG_H
