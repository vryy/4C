/*----------------------------------------------------------------------*/
/*!
\file

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "adapter_ale.H"
#include "adapter_ale_lin.H"
#include "adapter_ale_springs.H"

// further includes for AleBaseAlgorithm:
#include "../drt_lib/drt_globalproblem.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

using namespace std;
using namespace Teuchos;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::AleBaseAlgorithm::AleBaseAlgorithm()
{
  SetupAle();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::AleBaseAlgorithm::~AleBaseAlgorithm()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::AleBaseAlgorithm::SetupAle()
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("ADAPTER::AleBaseAlgorithm::SetupAle");
  Teuchos::TimeMonitor monitor(*t);

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numaf,0);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  RCP<IO::DiscretizationWriter> output =
    rcp(new IO::DiscretizationWriter(actdis));
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR        *actsolv  = &solv[genprob.numaf];

  const Teuchos::ParameterList& adyn     = DRT::Problem::Instance()->AleDynamicParams();
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RCP<ParameterList> solveparams = rcp(new ParameterList());
  RCP<LINALG::Solver> solver =
    rcp(new LINALG::Solver(solveparams,actdis->Comm(),allfiles.out_err));
  solver->TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  RCP<ParameterList> params = rcp(new ParameterList());
  params->set<int>("numstep",    fsidyn.get<int>("NUMSTEP"));
  params->set<double>("maxtime", fsidyn.get<double>("MAXTIME"));
  params->set<double>("dt",      fsidyn.get<double>("TIMESTEP"));

  // ----------------------------------------------- restart and output
  // restart
  params->set<int>("write restart every", fsidyn.get<int>("RESTARTEVRY"));

  int aletype = Teuchos::getIntegralValue<int>(adyn,"ALE_TYPE");
  if (aletype==ALE_DYNAMIC::classic_lin)
    ale_ = rcp(new AleLinear(actdis, solver, params, output));
  else if (aletype==ALE_DYNAMIC::springs)
    ale_ = rcp(new AleSprings(actdis, solver, params, output));
  else
    dserror("ale type '%s' unsupported",adyn.get<std::string>("ALE_TYPE").c_str());
}


#endif
