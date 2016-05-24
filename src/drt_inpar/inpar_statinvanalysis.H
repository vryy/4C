/*----------------------------------------------------------------------*/
/*!
\file inpar_statinvanalysis.H

\brief enums for the inverse analysis

<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@lnm.mw.tum.de
            089-289-10361
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_STATINVANALYSIS_H
#define INPAR_STATINVANALYSIS_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for Inverse Analysis
  namespace INVANA
  {
    enum StatInvRegularization
    {
      stat_inv_reg_none,              // nothing defined
      stat_inv_reg_tikhonov,          // thikonov regularization
      stat_inv_reg_totalvariation     // total variation diminishing regularization
    };

    // Parametrization of material parameters
    enum StatInvMatParametrization
    {
      stat_inv_mp_none,                  // nothing defined
      stat_inv_mp_smoothkernel,          // smooth kernel approximation
      stat_inv_mp_elementwise,            // parameters per element
      stat_inv_mp_uniform                // uniform i.e spatially constant parameters
    };

    /// Type of Inverse Analysis
    enum StatInvAnalysisType
    {
      stat_inv_none,          // no
      stat_inv_mc,            // statistical inverse analysis using montecarlo
      stat_inv_graddesc,      // deterministic inverse analysis using gradient descent
      stat_inv_lbfgs          // deterministic inverse analysis using limited memory bfgs
    };

    /// Type of the objective function
    enum StatInvObjFunctType
    {
      stat_inv_obj_none,         // none
      stat_inv_obj_disp,         // based on per dof displacement measurements
      stat_inv_obj_surfcurr      // based on a current representation of surfaces
    };

    /// Type of the metaparameterization
    enum StatInvMetaParamsType
    {
      stat_inv_meta_none,        // no meta parameterization
      stat_inv_meta_quad,        // parabolic meta parametrization
      stat_inv_meta_arctan       // arcus tangens meta parametrization
    };

  }  // namespace INVANA

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_STATINVANALYSIS_H
