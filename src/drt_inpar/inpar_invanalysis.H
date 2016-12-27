/*----------------------------------------------------------------------*/
/*!
\file inpar_invanalysis.H

\brief Input parameters for inverse analysis

\level 3

<pre>
\maintainer Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_INVANALYSIS_H
#define INPAR_INVANALYSIS_H

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT { namespace INPUT { class ConditionDefinition; } }

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for Inverse Analysis
  namespace STR
  {

    /// Type of Inverse Analysis
    enum InvAnalysisType
    {
      inv_none,          // no inverse analysis
      inv_generalized    // generalized inverse analysis
    };
    /// Update strategy for regularization parameter mu
    enum RegStratUpdate
      {
       reg_update_grad,          // Update based on gradient of objective function
       reg_update_res,          //  Update based on objective function
      };

  }  // namespace STR

  namespace INVANA
  {
    /// set the invana parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific inverse analysis conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

  }

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_INVANALYSIS_H
