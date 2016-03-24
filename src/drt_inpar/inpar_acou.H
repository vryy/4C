/*----------------------------------------------------------------------*/
/*!
\file inpar_acou.H

<pre>
\maintainer Svenja Schoeder
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
      acou_ssprk,       // strong stability preserving Runge-Kutta
      acou_ader         // ader time integration
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
      pat_none,       // no inverse analysis
      pat_opti,       // only for absorption
      pat_optisplit,  // for absorption and diffusion
      pat_optisplitacousplit, // for absorption diffusion, sos and density
      pat_optisplitacouident // use optical parameters to identify acoustical ones (what julia fischer did in her master's thesis)
    };

    /// Type of Algorithm
    enum OptimizationType
    {
      inv_lbfgs,        // limited memory BFGS
      inv_gd            // gradient descent
    };

    /// Type of Regularization
    enum RegulaType
    {
      pat_regula_none,
      pat_regula_tikh,        // tikhonov
      pat_regula_tvd,         // total variation
      pat_regula_tikhtvd      // both
    };


    /// Defines all valid parameters for acou problem
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific acou conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

  } // namespace ACOU
} // namespace INPAR


#endif // #ifndef INPAR_ACOU_H
