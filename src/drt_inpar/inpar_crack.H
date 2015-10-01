/*----------------------------------------------------------------------*/
/*!
\file inpar_crack.H

\brief Input parameters for crack propagation

<pre>
Maintainer: sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_CRACK_H
#define INPAR_CRACK_H

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT { namespace INPUT { class ConditionDefinition; } }

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for Inverse Analysis
  namespace CRACK
  {

    enum crackModel
    {
      crack_none,
      crack_lefm,
      crack_cohesive
    };

    /// Type of cohesive crack modeling
    enum CohesiveModel
    {
      cohesive_none,          // no carack propagation
      cohesive_dczm,          // Discrete cohesive zone modeling
      cohesive_ddzm           // Discrete damage zone modeling
    };

    enum tractionSeparation
    {
      linear,
      trapezoidal,
      exponential,
      sinusoidal,
      ppr
    };

    enum propagationCriterion
    {
      displacementCorrelation,
      J_Integral
    };

    /// set the crack parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific crack conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

  }  // namespace CRACK

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif
