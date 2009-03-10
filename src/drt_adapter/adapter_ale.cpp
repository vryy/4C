/*----------------------------------------------------------------------*/
/*!
 * \file adapter_ale.cpp
 * 
\brief ALE base implementation

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
#include "adapter_ale_laplace.H"
#include "adapter_ale_springs.H"
#include "adapter_ale_springs_fixed_ref.H"

// further includes for AleBaseAlgorithm:
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "../drt_inpar/inpar_fsi.H"
#include "../drt_fluid/drt_periodicbc.H"

using namespace std;
using namespace Teuchos;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


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
  // connect degrees of freedom for coupled nodes
  // -------------------------------------------------------------------
  PeriodicBoundaryConditions::PeriodicBoundaryConditions pbc(actdis);
  pbc.UpdateDofsForPeriodicBoundaryConditions();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  RCP<IO::DiscretizationWriter> output =
    rcp(new IO::DiscretizationWriter(actdis));
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  const Teuchos::ParameterList& adyn     = DRT::Problem::Instance()->AleDynamicParams();
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RCP<LINALG::Solver> solver =
    rcp(new LINALG::Solver(DRT::Problem::Instance()->AleSolverParams(),
                           actdis->Comm(),
                           DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  RCP<ParameterList> params = rcp(new ParameterList());
  params->set<int>("numstep",    fsidyn.get<int>("NUMSTEP"));
  params->set<double>("maxtime", fsidyn.get<double>("MAXTIME"));
  params->set<double>("dt",      fsidyn.get<double>("TIMESTEP"));

  // ----------------------------------------------- restart and output
  // restart
  params->set<int>("write restart every", fsidyn.get<int>("RESTARTEVRY"));

  params->set<int>("ALE_TYPE",      Teuchos::getIntegralValue<int>(adyn,"ALE_TYPE"));

  
  bool dirichletcond = true;
  if (genprob.probtyp == prb_fsi)
  {
    // FSI input parameters
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    int coupling = Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO");
    if (coupling == fsi_iter_monolithic or
        coupling == fsi_iter_monolithiclagrange or
        coupling == fsi_iter_monolithicstructuresplit)
    {
      // partitioned MFSI solvers require Dirichlet conditions
      INPAR::FSI::LinearBlockSolver linearsolverstrategy =
        Teuchos::getIntegralValue<INPAR::FSI::LinearBlockSolver>(fsidyn,"LINEARBLOCKSOLVER");
      if (linearsolverstrategy==INPAR::FSI::PartitionedAitken or
          linearsolverstrategy==INPAR::FSI::PartitionedVectorExtrapolation or
          linearsolverstrategy==INPAR::FSI::PartitionedJacobianFreeNewtonKrylov)
        dirichletcond = true;
      else
        dirichletcond = false;
    }
  }

  int aletype = Teuchos::getIntegralValue<int>(adyn,"ALE_TYPE");
  if (aletype==ALE_DYNAMIC::classic_lin)
    ale_ = rcp(new AleLinear(actdis, solver, params, output, false, dirichletcond));
  else if (aletype==ALE_DYNAMIC::incr_lin)
    ale_ = rcp(new AleLinear(actdis, solver, params, output, true , dirichletcond));
  else if (aletype==ALE_DYNAMIC::laplace)
    ale_ = rcp(new AleLaplace(actdis, solver, params, output, true, dirichletcond));
  else if (aletype==ALE_DYNAMIC::springs)
    ale_ = rcp(new AleSprings(actdis, solver, params, output, dirichletcond));
  else if (aletype==ALE_DYNAMIC::springs_fixed_ref)
    ale_ = rcp(new AleSpringsFixedRef(actdis, solver, params, output, true, dirichletcond));
  else
    dserror("ale type '%s' unsupported",adyn.get<std::string>("ALE_TYPE").c_str());
}


#endif
