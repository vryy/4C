/*----------------------------------------------------------------------*/
/*!
\file inpar_acou.H

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_ACOU_H
#define INPAR_ACOU_H

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT { namespace INPUT { class ConditionDefinition; } }

/*----------------------------------------------------------------------*/
namespace INPAR
{
  // Input parameters for acoustical simulations
  namespace ACOU
  {
    // Type of time integrator
    enum DynamicType
    {
      acou_impleuler,   // implicit Euler
      acou_expleuler,   // explicit Euler
      acou_classrk4,    // classical explicit Runge-Kutta with 4 stages
      acou_lsrk45reg2,  // low storage Runge-Kutta
      acou_lsrk33reg2,  // low storage Runge-Kutta
      acou_lsrk45reg3,  // low storage Runge-Kutta
      acou_ssprk        // strong stability preserving Runge-Kutta
    };

    // physical type
    enum PhysicalType
    {
      acou_lossless,
      acou_solid
    };

    /// Type of Inverse Analysis
    enum InvAnalysisType
    {
      inv_none,          // no inverse analysis
      inv_pat            // photoacoustic inverse analysis
    };

    /// Type of Algorithm
    enum OptimizationType
    {
      inv_lbfgs,        // limited memory BFGS
      inv_gd            // gradient descent
    };

    /// Defines all valid parameters for acou problem
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific acou conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

  } // namespace ACOU
} // namespace INPAR


#endif // #ifndef INPAR_ACOU_H
