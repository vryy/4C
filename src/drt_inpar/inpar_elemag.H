/*----------------------------------------------------------------------*/
/*!
\file inpar_elemag.H

\brief Input parameters for electromagnetic simulations

<pre>
\level 3

\maintainer Volker Gravemeier
            gravemeier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
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
  namespace ELEMAG
  {
    // Type of time integrator
    enum DynamicType
    {
      elemag_ost  // one-step-theta time integration
    };

    /// define all valid parameters for electromagnetic problem
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific electromagnetic conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

  }  // namespace ELEMAG
}  // namespace INPAR


#endif  // #ifndef INPAR_ELEMAG_H
