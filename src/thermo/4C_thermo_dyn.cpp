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
  Teuchos::RCP<Core::FE::Discretization> thermodis = Teuchos::null;
  thermodis = Global::Problem::instance()->get_dis("thermo");

  // set degrees of freedom in the discretization
  if (not thermodis->filled()) thermodis->fill_complete();

  const Teuchos::ParameterList& tdyn = Global::Problem::instance()->thermal_dynamic_params();

  // create instance of thermo basis algorithm (no structure discretization)
  Teuchos::RCP<Adapter::ThermoBaseAlgorithm> thermoonly =
      Teuchos::rcp(new Adapter::ThermoBaseAlgorithm(tdyn, thermodis));

  // do restart if demanded from input file
  const int restart = Global::Problem::instance()->restart();
  if (restart)
  {
    thermoonly->thermo_field().read_restart(restart);
  }

  // enter time loop to solve problem
  (thermoonly->thermo_field()).integrate();

  // perform the result test if required
  Global::Problem::instance()->add_field_test(thermoonly->thermo_field().create_field_test());
  Global::Problem::instance()->test_all(thermodis->get_comm());

  // done
  return;

}  // thr_dyn_drt()


/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
