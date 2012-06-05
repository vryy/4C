/*----------------------------------------------------------------------*/
/*!
 * \file ale.cpp
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



#include "ale.H"
#include "ale_lin.H"
#include "ale_laplace.H"
#include "ale_springs.H"
#include "ale_springs_fixed_ref.H"

// further includes for AleBaseAlgorithm:
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "../drt_inpar/inpar_ale.H"
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
ALE::AleBaseAlgorithm::AleBaseAlgorithm(const Teuchos::ParameterList& prbdyn, int disnum)
{
  SetupAle(prbdyn,disnum);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ALE::AleBaseAlgorithm::~AleBaseAlgorithm()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ALE::AleBaseAlgorithm::SetupAle(const Teuchos::ParameterList& prbdyn, int disnum)
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("ALE::AleBaseAlgorithm::SetupAle");
  Teuchos::TimeMonitor monitor(*t);

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numaf,disnum);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // connect degrees of freedom for coupled nodes
  // -------------------------------------------------------------------
  PeriodicBoundaryConditions pbc(actdis);
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

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // get the solver number used for ALE problems
    const int linsolvernumber = adyn.get<int>("LINEAR_SOLVER");
    // check if the TSI solver has a valid solver number
    if (linsolvernumber == (-1))
      dserror("no linear solver defined for ALE problems. Please set LINEAR_SOLVER in ALE DYNAMIC to a valid number!");

  RCP<LINALG::Solver> solver =
    rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                           actdis->Comm(),
                           DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  RCP<ParameterList> params = rcp(new ParameterList());
  params->set<int>("numstep",    prbdyn.get<int>("NUMSTEP"));
  params->set<double>("maxtime", prbdyn.get<double>("MAXTIME"));
  params->set<double>("dt",      prbdyn.get<double>("TIMESTEP"));

  // ----------------------------------------------- restart and output
  // restart
  params->set<int>("write restart every", prbdyn.get<int>("RESTARTEVRY"));

  params->set<int>("ALE_TYPE",DRT::INPUT::IntegralValue<int>(adyn,"ALE_TYPE"));


  bool dirichletcond = true;
  // what's the current problem type?
  PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();
  if (probtype == prb_fsi or
      probtype == prb_fsi_lung or
      probtype == prb_gas_fsi or
      probtype == prb_thermo_fsi or
      probtype == prb_biofilm_fsi or
      probtype == prb_fluid_fluid_fsi)
  {
    // FSI input parameters
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
    if (coupling == fsi_iter_monolithicfluidsplit or
        coupling == fsi_iter_monolithicstructuresplit or
        coupling == fsi_iter_constr_monolithicfluidsplit or
        coupling == fsi_iter_constr_monolithicstructuresplit or
        coupling == fsi_iter_lung_monolithicfluidsplit or
        coupling == fsi_iter_lung_monolithicstructuresplit or
        coupling == fsi_iter_mortar_monolithicstructuresplit or
        coupling == fsi_iter_mortar_monolithicfluidsplit or
        coupling == fsi_iter_fluidfluid_monolithicstructuresplit)
    {
        dirichletcond = false;
    }
  }

  if (probtype == prb_freesurf)
  {
    // FSI input parameters
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
    if (coupling == fsi_iter_monolithicfluidsplit or
         coupling == fsi_iter_monolithicstructuresplit or
         coupling == fsi_iter_constr_monolithicfluidsplit or
         coupling == fsi_iter_constr_monolithicstructuresplit or
         coupling == fsi_iter_lung_monolithicfluidsplit or
         coupling == fsi_iter_lung_monolithicstructuresplit or
         coupling == fsi_iter_mortar_monolithicstructuresplit or
         coupling == fsi_iter_mortar_monolithicfluidsplit)
    {
      dirichletcond = false;
    }
  }

  int aletype = DRT::INPUT::IntegralValue<int>(adyn,"ALE_TYPE");
  if (aletype==INPAR::ALE::classic_lin)
    ale_ = rcp(new AleLinear(actdis, solver, params, output, false, dirichletcond));
  else if (aletype==INPAR::ALE::incr_lin)
    ale_ = rcp(new AleLinear(actdis, solver, params, output, true , dirichletcond));
  else if (aletype==INPAR::ALE::laplace)
    ale_ = rcp(new AleLaplace(actdis, solver, params, output, true, dirichletcond));
  else if (aletype==INPAR::ALE::springs)
    ale_ = rcp(new AleSprings(actdis, solver, params, output, dirichletcond));
  else if (aletype==INPAR::ALE::springs_fixed_ref)
    ale_ = rcp(new AleSpringsFixedRef(actdis, solver, params, output, true, dirichletcond));
  else
    dserror("ale type '%s' unsupported",adyn.get<std::string>("ALE_TYPE").c_str());
}


