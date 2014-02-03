/*----------------------------------------------------------------------*/
/*!
\file ale_dyn.cpp

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/


#include "ale_dyn.H"
#include "ale.H"
#include "ale_resulttest.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include <Teuchos_RCP.hpp>


void dyn_ale_drt()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> actdis = DRT::Problem::Instance()->GetDis("ale");

  // -------------------------------------------------------------------
  // ask ALE::AleBaseAlgorithm for the ale time integrator
  // -------------------------------------------------------------------
  Teuchos::RCP<ALE::AleBaseAlgorithm> ale
    = Teuchos::rcp(new ALE::AleBaseAlgorithm(DRT::Problem::Instance()->AleDynamicParams(), actdis));
  Teuchos::RCP<ALE::Ale> aletimint = ale->AleFieldrcp();

  // -------------------------------------------------------------------
  // read the restart information, set vectors and variables if necessary
  // -------------------------------------------------------------------
  const int restart = DRT::Problem::Instance()->Restart();
  if (restart)
    aletimint->ReadRestart(restart);

  // -------------------------------------------------------------------
  // call time loop
  // -------------------------------------------------------------------
  aletimint->BuildSystemMatrix();
  aletimint->Integrate();

  // -------------------------------------------------------------------
  // do the result test
  // -------------------------------------------------------------------
  DRT::Problem::Instance()->AddFieldTest(Teuchos::rcp(new ALE::AleResultTest(*aletimint)));
  DRT::Problem::Instance()->TestAll(actdis->Comm());
}


