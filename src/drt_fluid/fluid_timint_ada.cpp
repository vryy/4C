/*!
\file fluid_timint_ada.cpp

\brief Fluid wrapper for adaptive time stepping

<pre>
Maintainer: Raffaela Kruse
            kruse@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_timint_ada.H"

#include "fluidimplicitintegration.H"

void FLD::TimIntAda::Indicate
(
  bool& accepted,
  double& dt_new
)
{
  //Todo: estimation based on local error evaluation not yet implemented
  accepted = estimator_->TimeStepAccepted();
  dt_new = estimator_->CalculateDt();
}


double FLD::TimIntAda::EstimatorCFL::CalculateDt(
  const double norm///< current norm of local discretization error
)
{
  return fluidimplicitintegration_->EvaluateDtViaCflIfApplicable();
}
