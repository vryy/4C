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
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void dyn_ale_drt()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> actdis = GLOBAL::Problem::Instance()->GetDis("ale");

  // -------------------------------------------------------------------
  // ask ALE::AleBaseAlgorithm for the ale time integrator
  // -------------------------------------------------------------------
  Teuchos::RCP<ADAPTER::AleBaseAlgorithm> ale = Teuchos::rcp(
      new ADAPTER::AleBaseAlgorithm(GLOBAL::Problem::Instance()->AleDynamicParams(), actdis));
  Teuchos::RCP<ADAPTER::Ale> aletimint = ale->ale_field();

  // -------------------------------------------------------------------
  // read the restart information, set vectors and variables if necessary
  // -------------------------------------------------------------------
  const int restart = GLOBAL::Problem::Instance()->Restart();
  if (restart) aletimint->read_restart(restart);

  // -------------------------------------------------------------------
  // call time loop
  // -------------------------------------------------------------------
  aletimint->create_system_matrix();
  aletimint->Integrate();

  // -------------------------------------------------------------------
  // do the result test
  // -------------------------------------------------------------------
  // test results
  GLOBAL::Problem::Instance()->AddFieldTest(aletimint->CreateFieldTest());
  GLOBAL::Problem::Instance()->TestAll(actdis->Comm());

  return;
}

FOUR_C_NAMESPACE_CLOSE
