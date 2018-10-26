/*----------------------------------------------------------------------*/
/*!
\file inpar_acou.H
\level 2
<pre>
\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_ACOU_H
#define INPAR_ACOU_H

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
      acou_ader,        // ader time integration
      acou_ader_lts,    // ader with local time stepping
      acou_ader_tritet
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
      pat_none,                // no inverse analysis
      pat_opti,                // only for absorption
      pat_optisplit,           // for absorption and diffusion
      pat_optisplitacousplit,  // for absorption diffusion, sos and density
      pat_optisplitacouident,  // use optical parameters to identify acoustical ones (what julia
                               // fischer did in her master's thesis)
      pat_optisplitacouidentsmart,  // use optical parameters to identify acoustical ones but with
                                    // advanced decisions for ambigious materials
      pat_reduction  // to reduce the computational domain on a ring (e.g. a mouse completely
                     // surrounded by lots of water)
    };

    ///  list of all tomograph types (important for interpretation of measurement data)
    enum TomographType
    {
      pat_circle,  // detectors on circle
      pat_unknown  // to be filled in
    };

    /// Type of Algorithm
    enum OptimizationType
    {
      inv_lbfgs,  // limited memory BFGS
      inv_gd      // gradient descent
    };

    /// Type of Regularization
    enum RegulaType
    {
      pat_regula_none,
      pat_regula_tikh,    // tikhonov
      pat_regula_tvd,     // total variation
      pat_regula_tikhtvd  // both
    };

    /// Patch Type
    enum PatchType
    {
      pat_patch_none,
      pat_patch_self,
      pat_patch_reacgrad,
      pat_patch_reacvals,
      pat_patch_mixed
    };

    /// Defines all valid parameters for acou problem
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific acou conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

  }  // namespace ACOU
}  // namespace INPAR


#endif  // #ifndef INPAR_ACOU_H
