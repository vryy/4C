/*----------------------------------------------------------------------*/
/*! \file
\brief entry point for (in)stationary heat conduction
\level 1
*/

/*----------------------------------------------------------------------*
 | definitions                                                gjb 01/08 |
 *----------------------------------------------------------------------*/

#include "baci_thermo_dyn.hpp"

#include "baci_adapter_thermo.hpp"
#include "baci_global_data.hpp"
#include "baci_inpar_validparameters.hpp"
#include "baci_thermo_resulttest.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | main control routine for (in)stationary heat conduction              |
 *----------------------------------------------------------------------*/
void thr_dyn_drt()
{
  // access the discretization
  Teuchos::RCP<DRT::Discretization> thermodis = Teuchos::null;
  thermodis = GLOBAL::Problem::Instance()->GetDis("thermo");

  // set degrees of freedom in the discretization
  if (not thermodis->Filled()) thermodis->FillComplete();

  const Teuchos::ParameterList& tdyn = GLOBAL::Problem::Instance()->ThermalDynamicParams();

  // create instance of thermo basis algorithm (no structure discretization)
  Teuchos::RCP<ADAPTER::ThermoBaseAlgorithm> thermoonly =
      Teuchos::rcp(new ADAPTER::ThermoBaseAlgorithm(tdyn, thermodis));

  // do restart if demanded from input file
  const int restart = GLOBAL::Problem::Instance()->Restart();
  if (restart)
  {
    thermoonly->ThermoField().ReadRestart(restart);
  }

  // enter time loop to solve problem
  (thermoonly->ThermoField()).Integrate();

  // perform the result test if required
  GLOBAL::Problem::Instance()->AddFieldTest(thermoonly->ThermoField().CreateFieldTest());
  GLOBAL::Problem::Instance()->TestAll(thermodis->Comm());

  // done
  return;

}  // thr_dyn_drt()


/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
