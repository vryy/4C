/*--------------------------------------------------------------------------*/
/*! \file

\brief Lubrication field base algorithm

\level 3

\maintainer Mostafa Faraji
*/
/*--------------------------------------------------------------------------*/

#include "adapter_lubrication.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io.H"

#include "../linalg/linalg_solver.H"

#include "../drt_lubrication/lubrication_resulttest.H"

#include "../drt_lubrication/lubrication_timint_stat.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::LubricationBaseAlgorithm::Setup(
    const Teuchos::ParameterList& prbdyn,          ///< parameter list for global problem
    const Teuchos::ParameterList& lubricationdyn,  ///< parameter list for Lubrication subproblem
    const Teuchos::ParameterList& solverparams,    ///< parameter list for Lubrication solver
    const std::string& disname,                    ///< name of Lubrication discretization
    const bool isale                               ///< ALE flag
)
{
  // setup Lubrication algorithm (overriding some dynamic parameters
  // with values specified in given problem-dependent ParameterList prbdyn)

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->GetDis(disname);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<IO::DiscretizationWriter> output = actdis->Writer();
  output->WriteMesh(0, 0.0);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // TODO: TAW use of solverparams??? change input parameter to solver number instead of parameter
  // list? -> no default paramter possible any more
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(new LINALG::Solver(
      solverparams, actdis->Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  // make a copy (inside an Teuchos::rcp) containing also all sublists
  Teuchos::RCP<Teuchos::ParameterList> lubricationtimeparams =
      Teuchos::rcp(new Teuchos::ParameterList(lubricationdyn));

  // -------------------------------------------------------------------
  // overrule certain parameters for coupled problems
  // -------------------------------------------------------------------
  // the default time step size
  lubricationtimeparams->set<double>("TIMESTEP", prbdyn.get<double>("TIMESTEP"));
  // maximum simulation time
  lubricationtimeparams->set<double>("MAXTIME", prbdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  lubricationtimeparams->set<int>("NUMSTEP", prbdyn.get<int>("NUMSTEP"));
  // restart
  lubricationtimeparams->set("RESTARTEVRY", prbdyn.get<int>("RESTARTEVRY"));
  // solution output
  lubricationtimeparams->set("RESULTSEVRY", prbdyn.get<int>("RESULTSEVRY"));

  // -------------------------------------------------------------------
  // list for extra parameters
  // (put here everything that is not available in lubricationdyn or its sublists)
  // -------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> extraparams = Teuchos::rcp(new Teuchos::ParameterList());

  // ------------------------------pointer to the error file (for output)
  extraparams->set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());

  // ----------------Eulerian or ALE formulation of transport equation(s)
  extraparams->set<bool>("isale", isale);

  // create instance of time integration class (call the constructor)
  lubrication_ = Teuchos::rcp(new LUBRICATION::TimIntStationary(
      actdis, solver, lubricationtimeparams, extraparams, output));

  lubrication_->Init();
  // initialize algorithm for specific time-integration scheme

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::LubricationBaseAlgorithm::CreateLubricationFieldTest()
{
  return Teuchos::rcp(new LUBRICATION::ResultTest(lubrication_));
}

Teuchos::RCP<IO::DiscretizationWriter> ADAPTER::LubricationBaseAlgorithm::DiscWriter()
{
  return lubrication_->DiscWriter();
}
