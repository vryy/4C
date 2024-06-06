/*----------------------------------------------------------------------*/
/*! \file
\brief entry point for (in)stationary heat conduction
\level 1
*/

/*----------------------------------------------------------------------*
 | definitions                                                gjb 01/08 |
 *----------------------------------------------------------------------*/

#include "4C_thermo_dyn.hpp"

#include "4C_adapter_thermo.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_thermo_resulttest.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | main control routine for (in)stationary heat conduction              |
 *----------------------------------------------------------------------*/
void thr_dyn_drt()
{
  // access the discretization
  Teuchos::RCP<Discret::Discretization> thermodis = Teuchos::null;
  thermodis = Global::Problem::Instance()->GetDis("thermo");

  // set degrees of freedom in the discretization
  if (not thermodis->Filled()) thermodis->fill_complete();

  const Teuchos::ParameterList& tdyn = Global::Problem::Instance()->thermal_dynamic_params();

  // create instance of thermo basis algorithm (no structure discretization)
  Teuchos::RCP<Adapter::ThermoBaseAlgorithm> thermoonly =
      Teuchos::rcp(new Adapter::ThermoBaseAlgorithm(tdyn, thermodis));

  // do restart if demanded from input file
  const int restart = Global::Problem::Instance()->Restart();
  if (restart)
  {
    thermoonly->ThermoField().read_restart(restart);
  }

  // enter time loop to solve problem
  (thermoonly->ThermoField()).Integrate();

  // perform the result test if required
  Global::Problem::Instance()->AddFieldTest(thermoonly->ThermoField().CreateFieldTest());
  Global::Problem::Instance()->TestAll(thermodis->Comm());

  // done
  return;

}  // thr_dyn_drt()


/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
