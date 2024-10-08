/*--------------------------------------------------------------------------*/
/*! \file

\brief Lubrication field base algorithm

\level 3

*/
/*--------------------------------------------------------------------------*/

#include "4C_lubrication_adapter.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_lubrication_resulttest.hpp"
#include "4C_lubrication_timint_stat.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LUBRICATION::LubricationBaseAlgorithm::setup(
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
  actdis = Global::Problem::instance()->get_dis(disname);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->filled()) actdis->fill_complete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<Core::IO::DiscretizationWriter> output = actdis->writer();
  output->write_mesh(0, 0.0);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // TODO: TAW use of solverparams??? change input parameter to solver number instead of parameter
  // list? -> no default paramter possible any more
  Teuchos::RCP<Core::LinAlg::Solver> solver = Teuchos::make_rcp<Core::LinAlg::Solver>(solverparams,
      actdis->get_comm(), Global::Problem::instance()->solver_params_callback(),
      Teuchos::getIntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::instance()->io_params(), "VERBOSITY"));
  actdis->compute_null_space_if_necessary(solver->params());

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  // make a copy (inside an Teuchos::rcp) containing also all sublists
  Teuchos::RCP<Teuchos::ParameterList> lubricationtimeparams =
      Teuchos::make_rcp<Teuchos::ParameterList>(lubricationdyn);

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
  Teuchos::RCP<Teuchos::ParameterList> extraparams = Teuchos::make_rcp<Teuchos::ParameterList>();

  // ----------------Eulerian or ALE formulation of transport equation(s)
  extraparams->set<bool>("isale", isale);

  // create instance of time integration class (call the constructor)
  lubrication_ = Teuchos::make_rcp<LUBRICATION::TimIntStationary>(
      actdis, solver, lubricationtimeparams, extraparams, output);

  lubrication_->init();
  // initialize algorithm for specific time-integration scheme

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::UTILS::ResultTest>
LUBRICATION::LubricationBaseAlgorithm::create_lubrication_field_test()
{
  return Teuchos::make_rcp<LUBRICATION::ResultTest>(lubrication_);
}

Teuchos::RCP<Core::IO::DiscretizationWriter> LUBRICATION::LubricationBaseAlgorithm::disc_writer()
{
  return lubrication_->disc_writer();
}

FOUR_C_NAMESPACE_CLOSE
