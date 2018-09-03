/*----------------------------------------------------------------------*/
/*!
\file inpar_invanalysis.H

\brief Input parameters for inverse analysis

\level 2

<pre>
\maintainer Harald Willmann
            willmann@lnm.mw.tum.de
            089-289-15255
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_INVANALYSIS_H
#define INPAR_INVANALYSIS_H

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
  /// Input parameters for Inverse Analysis
  namespace STR
  {
    /// Type of Inverse Analysis
    enum InvAnalysisType
    {
      inv_none,        // no inverse analysis
      inv_generalized  // generalized inverse analysis
    };

    /// Type of Special Inverse Analysis
    enum SpecInvAnalysisType
    {
      spec_inv_none,  // no special inverse analysis
      spec_inv_coup,  // coupled inverse analysis: use information of two different experiments,
                      // here implemented for pressure-volume-change experiment and tensile test
                      // (just applicable for lung)
      spec_inv_mult   // multible inverse analysis: fit more inverse analysis (different dat-files)
                      // in one identification procedure
    };

    /// Update strategy for regularization parameter mu
    enum RegStratUpdate
    {
      reg_update_grad,  // Update based on gradient of objective function
      reg_update_res,   //  Update based on objective function
    };

  }  // namespace STR

  namespace INVANA
  {
    /// set the invana parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific inverse analysis conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

  }  // namespace INVANA

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_INVANALYSIS_H
