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


#include "../drt_io/io_control.H"
#include "ale_dyn.H"
#include "ale.H"
#include "ale_lin.H"
#include "ale_resulttest.H"

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Epetra_Time.h>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_io/io.H"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;



using namespace std;
using namespace Teuchos;



void dyn_ale_drt()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numaf,0);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  RefCountPtr<IO::DiscretizationWriter> output =
    rcp(new IO::DiscretizationWriter(actdis));
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  //const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& adyn     = DRT::Problem::Instance()->AleDynamicParams();

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // get the solver number used for ALE problems
  const int linsolvernumber = adyn.get<int>("LINEAR_SOLVER");
  // check if the TSI solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for ALE problems. Please set LINEAR_SOLVER in ALE DYNAMIC to a valid number!");

  RefCountPtr<LINALG::Solver> solver =
    rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                           actdis->Comm(),
                           DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  RefCountPtr<ParameterList> params = rcp(new ParameterList());
  params->set<int>("numstep", adyn.get<int>("NUMSTEP"));
  params->set<double>("maxtime", adyn.get<double>("MAXTIME"));
  params->set<double>("dt", adyn.get<double>("TIMESTEP"));

  //int aletype = DRT::INPUT::IntegralValue<int>(adyn,"ALE_TYPE");
  ALE::AleLinear ale(actdis, solver, params, output, false, true);

  const int restart = DRT::Problem::Instance()->Restart();
  if (restart)
  {
    // read the restart information, set vectors and variables
    ale.ReadRestart(restart);
  }

  ale.BuildSystemMatrix();
  ale.Integrate();

  // do the result test
  DRT::Problem::Instance()->AddFieldTest(rcp(new ALE::AleResultTest(ale)));
  DRT::Problem::Instance()->TestAll(actdis->Comm());
}


