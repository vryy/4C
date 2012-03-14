/*!-----------------------------------------------------------------------------------------------*
\file combust_reinitialization_pde.H

\brief
<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef COMBUST_REINITIALIZATION_PDE_H
#define COMBUST_REINITIALIZATION_PDE_H

#include "../drt_scatra/scatra_timint_implicit.H"


namespace ADAPTER
{
class ScaTraBaseAlgorithm;
}

namespace COMBUST
{

class ReinitializationPDE
{
public:

  ReinitializationPDE();

  virtual ~ReinitializationPDE();

  bool CallReinitialization();

private:

  void TimeLoop_Reinit();

  void ReinitializeInfo();

 /* void CheckSteadyState(
      double&    Gradient_Error_old,
      bool&      STOP);
*/
//  double EvaluateGradientNormError();

  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> reinit_;

  Teuchos::RCP<Teuchos::ParameterList> combustdynreinit_;

  /// processor id
  int                        myrank_;

  int stepmax_;

}; // class ReinitializationPDE

} // namespace COMBUST

#endif //COMBUST_REINITIALIZATION_PDE_H
