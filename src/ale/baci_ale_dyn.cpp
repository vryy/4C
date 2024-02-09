/*----------------------------------------------------------------------------*/
/*! \file

\brief Entry routine for pure ALE problems

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "baci_ale_dyn.hpp"

#include "baci_adapter_ale.hpp"
#include "baci_ale_resulttest.hpp"
#include "baci_global_data.hpp"
#include "baci_lib_discret.hpp"

#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

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
  Teuchos::RCP<ADAPTER::Ale> aletimint = ale->AleField();

  // -------------------------------------------------------------------
  // read the restart information, set vectors and variables if necessary
  // -------------------------------------------------------------------
  const int restart = GLOBAL::Problem::Instance()->Restart();
  if (restart) aletimint->ReadRestart(restart);

  // -------------------------------------------------------------------
  // call time loop
  // -------------------------------------------------------------------
  aletimint->CreateSystemMatrix();
  aletimint->Integrate();

  // -------------------------------------------------------------------
  // do the result test
  // -------------------------------------------------------------------
  // test results
  GLOBAL::Problem::Instance()->AddFieldTest(aletimint->CreateFieldTest());
  GLOBAL::Problem::Instance()->TestAll(actdis->Comm());

  return;
}

BACI_NAMESPACE_CLOSE
