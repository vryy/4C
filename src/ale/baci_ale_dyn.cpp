/*----------------------------------------------------------------------------*/
/*! \file

\brief Entry routine for pure ALE problems

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "baci_ale_dyn.H"

#include "baci_adapter_ale.H"
#include "baci_ale_resulttest.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"

#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void dyn_ale_drt()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> actdis = DRT::Problem::Instance()->GetDis("ale");

  // -------------------------------------------------------------------
  // ask ALE::AleBaseAlgorithm for the ale time integrator
  // -------------------------------------------------------------------
  Teuchos::RCP<::ADAPTER::AleBaseAlgorithm> ale = Teuchos::rcp(
      new ::ADAPTER::AleBaseAlgorithm(DRT::Problem::Instance()->AleDynamicParams(), actdis));
  Teuchos::RCP<::ADAPTER::Ale> aletimint = ale->AleField();

  // -------------------------------------------------------------------
  // read the restart information, set vectors and variables if necessary
  // -------------------------------------------------------------------
  const int restart = DRT::Problem::Instance()->Restart();
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
  DRT::Problem::Instance()->AddFieldTest(aletimint->CreateFieldTest());
  DRT::Problem::Instance()->TestAll(actdis->Comm());

  return;
}
