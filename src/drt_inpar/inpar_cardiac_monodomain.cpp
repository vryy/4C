/*----------------------------------------------------------------------*/
/*!
\file inpar_cardiac_monodomain.cpp

\brief Input parameters for cardiac monodomain

\level 2

<pre>
\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15264
</pre>
*/
/*----------------------------------------------------------------------*/
#include "inpar_cardiac_monodomain.H"

#include "drt_validparameters.H"

#include "../drt_lib/drt_conditiondefinition.H"

void INPAR::EP::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& epcontrol = list->sublist("CARDIAC MONODOMAIN CONTROL", false,
      "control parameters for cardiac electrophysiology problems\n");


  // Parameters for reaction-diffusion systems (for example cardiac electrophysiology)
  IntParameter("WRITEMAXINTSTATE", 0,
      "number of maximal internal state variables to be postprocessed", &epcontrol);
  IntParameter("WRITEMAXIONICCURRENTS", 0, "number of maximal ionic currents to be postprocessed",
      &epcontrol);

  DoubleParameter("ACTTHRES", 1.0,
      "threshold for the potential for computing and postprocessing activation time ", &epcontrol);
}


void INPAR::EP::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
}
