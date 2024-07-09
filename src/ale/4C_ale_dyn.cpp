/*----------------------------------------------------------------------------*/
/*! \file

\brief Entry routine for pure ALE problems

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "4C_ale_dyn.hpp"

#include "4C_adapter_ale.hpp"
#include "4C_ale_resulttest.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void dyn_ale_drt()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<Core::FE::Discretization> actdis = Global::Problem::instance()->get_dis("ale");

  // -------------------------------------------------------------------
  // ask ALE::AleBaseAlgorithm for the ale time integrator
  // -------------------------------------------------------------------
  Teuchos::RCP<Adapter::AleBaseAlgorithm> ale = Teuchos::rcp(
      new Adapter::AleBaseAlgorithm(Global::Problem::instance()->ale_dynamic_params(), actdis));
  Teuchos::RCP<Adapter::Ale> aletimint = ale->ale_field();

  // -------------------------------------------------------------------
  // read the restart information, set vectors and variables if necessary
  // -------------------------------------------------------------------
  const int restart = Global::Problem::instance()->restart();
  if (restart) aletimint->read_restart(restart);

  // -------------------------------------------------------------------
  // call time loop
  // -------------------------------------------------------------------
  aletimint->create_system_matrix();
  aletimint->integrate();

  // -------------------------------------------------------------------
  // do the result test
  // -------------------------------------------------------------------
  // test results
  Global::Problem::instance()->add_field_test(aletimint->create_field_test());
  Global::Problem::instance()->test_all(actdis->get_comm());

  return;
}

FOUR_C_NAMESPACE_CLOSE
