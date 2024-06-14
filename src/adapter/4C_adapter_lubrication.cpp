/*--------------------------------------------------------------------------*/
/*! \file

\brief Lubrication field base algorithm

\level 3

*/
/*--------------------------------------------------------------------------*/

#include "4C_adapter_lubrication.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_lubrication_resulttest.hpp"
#include "4C_lubrication_timint_stat.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::LubricationBaseAlgorithm::Setup(
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
  Teuchos::RCP<Core::FE::Discretization> actdis = Teuchos::null;
  actdis = Global::Problem::Instance()->GetDis(disname);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->fill_complete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<Core::IO::DiscretizationWriter> output = actdis->Writer();
  output->WriteMesh(0, 0.0);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // TODO: TAW use of solverparams??? change input parameter to solver number instead of parameter
  // list? -> no default paramter possible any more
  Teuchos::RCP<Core::LinAlg::Solver> solver = Teuchos::rcp(new Core::LinAlg::Solver(solverparams,
      actdis->Comm(), Global::Problem::Instance()->solver_params_callback(),
      Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::Instance()->IOParams(), "VERBOSITY")));
  actdis->compute_null_space_if_necessary(solver->Params());

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
Teuchos::RCP<Core::UTILS::ResultTest>
Adapter::LubricationBaseAlgorithm::create_lubrication_field_test()
{
  return Teuchos::rcp(new LUBRICATION::ResultTest(lubrication_));
}

Teuchos::RCP<Core::IO::DiscretizationWriter> Adapter::LubricationBaseAlgorithm::DiscWriter()
{
  return lubrication_->DiscWriter();
}

FOUR_C_NAMESPACE_CLOSE
