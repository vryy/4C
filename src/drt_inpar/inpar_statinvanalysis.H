/*----------------------------------------------------------------------*/
/*!
\file inpar_statinvanalysis.H

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
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
  namespace STR
  {
    enum StatInvRegularization
    {
      stat_inv_reg_none,              // nothing defined
      stat_inv_reg_tikhonov,          // thikonov regularization
      stat_inv_reg_totalvariation,
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

  }  // namespace STR

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_STATINVANALYSIS_H
