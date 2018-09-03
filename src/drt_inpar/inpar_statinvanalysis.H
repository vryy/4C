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
      stat_inv_reg_none,           // nothing defined
      stat_inv_reg_tikhonov,       // thikonov regularization
      stat_inv_reg_totalvariation  // total variation diminishing regularization
    };

    // Parametrization of material parameters
    enum StatInvMatParametrization
    {
      stat_inv_mp_none,         // nothing defined
      stat_inv_mp_patchwise,    // patchwise parameter approximation
      stat_inv_mp_tvsvd,        // svd of total variation approximation
      stat_inv_mp_elementwise,  // parameters per element
      stat_inv_mp_uniform       // uniform i.e spatially constant parameters
    };

    /// Type of Inverse Analysis
    enum StatInvAnalysisType
    {
      stat_inv_none,        // no
      stat_inv_smc,         // sequential monte carlo
      stat_inv_mh,          // metropolis hastings
      stat_inv_lbfgs,       // deterministic inverse analysis using limited memory bfgs
      stat_inv_bruteforce,  // brute force sampling
      stat_inv_prediction   // particle based prediction
    };

    /// Type of the objective function
    enum StatInvObjFunctType
    {
      stat_inv_obj_none,     // none
      stat_inv_obj_disp,     // based on per dof displacement measurements
      stat_inv_obj_surfcurr  // based on a current representation of surfaces
    };

    /// Type of the metaparameterization
    enum StatInvMetaParamsType
    {
      stat_inv_meta_none,   // no meta parameterization
      stat_inv_meta_quad,   // parabolic meta parametrization
      stat_inv_meta_exp,    // log-normal distribution of material parameters
      stat_inv_meta_arctan  // arcus tangens meta parametrization
    };

    /// Type of the initial guess
    enum StatInvInitType
    {
      stat_inv_init_dat,  // initialization from the datfile
      stat_inv_init_map   // initialization from some MAP approximation
    };

    /// Type of graph weights in the elementwise parametrization
    enum StatInvGraphWeight
    {
      stat_inv_graph_area,  // weights from area of faces dividing elements
      stat_inv_graph_unity  // unit weights between nodes of the graph
    };

  }  // namespace INVANA

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_STATINVANALYSIS_H
