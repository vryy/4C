/*----------------------------------------------------------------------*/
/*!
\file thr_dyn.cpp
\brief entry point for (in)stationary heat conduction
\level 1
<pre>
\maintainer Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                                gjb 01/08 |
 *----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 | headers                                                   gjb 01/08 |
 *----------------------------------------------------------------------*/
#include "thr_dyn.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_adapter/adapter_thermo.H"
#include "thr_resulttest.H"

/*----------------------------------------------------------------------*
 | main control routine for (in)stationary heat conduction              |
 *----------------------------------------------------------------------*/
void thr_dyn_drt()
{
  // access the discretization
  Teuchos::RCP<DRT::Discretization> thermodis = Teuchos::null;
  thermodis = DRT::Problem::Instance()->GetDis("thermo");

  // set degrees of freedom in the discretization
  if (not thermodis->Filled()) thermodis->FillComplete();

  const Teuchos::ParameterList& tdyn = DRT::Problem::Instance()->ThermalDynamicParams();

  // create instance of thermo basis algorithm (no structure discretization)
  Teuchos::RCP<ADAPTER::ThermoBaseAlgorithm> thermoonly =
      Teuchos::rcp(new ADAPTER::ThermoBaseAlgorithm(tdyn, thermodis));

  // do restart if demanded from input file
  const int restart = DRT::Problem::Instance()->Restart();
  if (restart)
  {
    thermoonly->ThermoField().ReadRestart(restart);
  }

  // enter time loop to solve problem
  (thermoonly->ThermoField()).Integrate();

  // perform the result test if required
  DRT::Problem::Instance()->AddFieldTest(thermoonly->ThermoField().CreateFieldTest());
  DRT::Problem::Instance()->TestAll(thermodis->Comm());

  // done
  return;

}  // thr_dyn_drt()


/*----------------------------------------------------------------------*/
