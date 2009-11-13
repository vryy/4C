/*----------------------------------------------------------------------*/
/*!
\file thr_dyn.cpp
\brief entry point for (in)stationary heat conduction

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                               gjb 01/08 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <Teuchos_TimeMonitor.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

/*----------------------------------------------------------------------*
 |  headers                                                   gjb 01/08 |
 *----------------------------------------------------------------------*/
#include "thr_dyn.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_adapter/adapter_thermo.H"
#include "../drt_adapter/adapter_thermo_timint.H"
#include "thr_resulttest.H"

/*----------------------------------------------------------------------*
 | general problem data                                     m.gee 06/01 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 | Main control routine for (in)stationary heat conduction              |
 *----------------------------------------------------------------------*/
void thr_dyn_drt()
{
  // access the discretization
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numtf, 0);

  // set degrees of freedom in the discretization
  if (not actdis->Filled()) actdis->FillComplete();

  // context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter> output
    = Teuchos::rcp(new IO::DiscretizationWriter(actdis));

  // get input parameter lists
  //const Teuchos::ParameterList& probtype
  //  = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& ioflags
    = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& tdyn
    = DRT::Problem::Instance()->ThermalDynamicParams();

  // show default parameters
  if (actdis->Comm().MyPID() == 0)
    DRT::INPUT::PrintDefaultParameters(std::cout, tdyn);

  // add extra parameters (a kind of work-around)
  Teuchos::ParameterList xparams;
  xparams.set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());

  // create a solver
  Teuchos::RCP<LINALG::Solver> solver
    = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->ThermalSolverParams(),
                                      actdis->Comm(),
                                      DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // create marching time integrator
  Teuchos::RCP<ADAPTER::Thermo> atti
    = Teuchos::rcp(new ADAPTER::ThermoTimInt(Teuchos::rcp(new Teuchos::ParameterList(ioflags)),
                                             Teuchos::rcp(new Teuchos::ParameterList(tdyn)),
                                             Teuchos::rcp(new Teuchos::ParameterList(xparams)),
                                             actdis, solver, output));
  if (atti == Teuchos::null) dserror("Failed in creating integrator.");

  // do restart if demanded from input file
  if (genprob.restart)
  {
    atti->ReadRestart(genprob.restart);
  }

  // write mesh always at beginning of calc or restart
  {
    int step = atti->GetTimeStep();
    double time = atti->GetTime();
    output->WriteMesh(step, time);
  }

  // integrate in time
  atti->Integrate();

  // test results
  DRT::Problem::Instance()->AddFieldTest(atti->CreateFieldTest());
  DRT::Problem::Instance()->TestAll(actdis->Comm());

  // done
  return;

} // end of thr_dyn_drt()

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
