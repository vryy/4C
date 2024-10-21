// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_thermo_adapter.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_thermo.hpp"
#include "4C_io_pstream.hpp"
#include "4C_thermo_timint_expleuler.hpp"
#include "4C_thermo_timint_genalpha.hpp"
#include "4C_thermo_timint_ost.hpp"
#include "4C_thermo_timint_statics.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |                                                          bborn 08/09 |
 *----------------------------------------------------------------------*/
Thermo::BaseAlgorithm::BaseAlgorithm(
    const Teuchos::ParameterList& prbdyn, Teuchos::RCP<Core::FE::Discretization> actdis)
{
  setup_thermo(prbdyn, actdis);
}



/*----------------------------------------------------------------------*
 |                                                          bborn 08/09 |
 *----------------------------------------------------------------------*/
void Thermo::BaseAlgorithm::setup_thermo(
    const Teuchos::ParameterList& prbdyn, Teuchos::RCP<Core::FE::Discretization> actdis)
{
  const Teuchos::ParameterList& tdyn = Global::Problem::instance()->thermal_dynamic_params();

  // major switch to different time integrators
  auto timinttype = Teuchos::getIntegralValue<Inpar::Thermo::DynamicType>(tdyn, "DYNAMICTYP");
  switch (timinttype)
  {
    case Inpar::Thermo::dyna_statics:
    case Inpar::Thermo::dyna_onesteptheta:
    case Inpar::Thermo::dyna_genalpha:
    case Inpar::Thermo::dyna_expleuler:
      setup_tim_int(prbdyn, timinttype, actdis);  // <-- here is the show
      break;
    default:
      FOUR_C_THROW("unknown time integration scheme '%s'",
          Teuchos::getStringValue<Inpar::Thermo::DynamicType>(tdyn, "DYNAMICTYP").c_str());
  }

}  // setup_thermo()


/*----------------------------------------------------------------------*
 | setup of thermal time integration                        bborn 08/09 |
 *----------------------------------------------------------------------*/
void Thermo::BaseAlgorithm::setup_tim_int(const Teuchos::ParameterList& prbdyn,
    Inpar::Thermo::DynamicType timinttype, Teuchos::RCP<Core::FE::Discretization> actdis)
{
  // this is not exactly a one hundred meter race, but we need timing
  Teuchos::RCP<Teuchos::Time> t =
      Teuchos::TimeMonitor::getNewTimer("Thermo::BaseAlgorithm::setup_thermo");
  Teuchos::TimeMonitor monitor(*t);

  // set degrees of freedom in the discretization
  if (not actdis->filled()) actdis->fill_complete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<Core::IO::DiscretizationWriter> output = actdis->writer();
  output->write_mesh(0, 0.0);

  //  // get input parameter lists and copy them, because a few parameters are overwritten
  Teuchos::ParameterList ioflags(Global::Problem::instance()->io_params());
  const Teuchos::RCP<Teuchos::ParameterList> tdyn = Teuchos::make_rcp<Teuchos::ParameterList>(
      Global::Problem::instance()->thermal_dynamic_params());
  //  //const Teuchos::ParameterList& size
  //  //  = Global::Problem::instance()->ProblemSizeParams();

  // add extra parameters (a kind of work-around)
  Teuchos::ParameterList xparams;

  // -------------------------------------------------------------------
  // overrule certain parameters for coupled problems
  // -------------------------------------------------------------------

  // the default time step size
  tdyn->set<double>("TIMESTEP", prbdyn.get<double>("TIMESTEP"));
  // maximum simulation time
  tdyn->set<double>("MAXTIME", prbdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  tdyn->set<int>("NUMSTEP", prbdyn.get<int>("NUMSTEP"));
  // restart
  tdyn->set<int>("RESTARTEVRY", prbdyn.get<int>("RESTARTEVRY"));
  // write results every <RESULTSEVRY> steps
  tdyn->set<int>("RESULTSEVRY", prbdyn.get<int>("RESULTSEVRY"));

  // get the solver number used for thermal solver
  const int linsolvernumber = tdyn->get<int>("LINEAR_SOLVER");
  // check if the THERMAL solver has a valid solver number
  if (linsolvernumber == (-1))
    FOUR_C_THROW(
        "no linear solver defined for thermal solver. Please set LINEAR_SOLVER in THERMAL DYNAMIC "
        "to a valid number!");

  // create a linear solver
  Teuchos::ParameterList solveparams;
  Teuchos::RCP<Core::LinAlg::Solver> solver = Teuchos::make_rcp<Core::LinAlg::Solver>(
      Global::Problem::instance()->solver_params(linsolvernumber), actdis->get_comm(),
      Global::Problem::instance()->solver_params_callback(),
      Teuchos::getIntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::instance()->io_params(), "VERBOSITY"));
  actdis->compute_null_space_if_necessary(solver->params());

  // create marching time integrator
  Teuchos::RCP<Adapter> tmpthr;
  switch (timinttype)
  {
    case Inpar::Thermo::dyna_statics:
    {
      tmpthr =
          Teuchos::make_rcp<Thermo::TimIntStatics>(ioflags, *tdyn, xparams, actdis, solver, output);
      break;
    }
    case Inpar::Thermo::dyna_onesteptheta:
    {
      tmpthr = Teuchos::make_rcp<Thermo::TimIntOneStepTheta>(
          ioflags, *tdyn, xparams, actdis, solver, output);
      break;
    }
    case Inpar::Thermo::dyna_genalpha:
    {
      tmpthr = Teuchos::make_rcp<Thermo::TimIntGenAlpha>(
          ioflags, *tdyn, xparams, actdis, solver, output);
      break;
    }
    case Inpar::Thermo::dyna_expleuler:
    {
      tmpthr = Teuchos::make_rcp<Thermo::TimIntExplEuler>(
          ioflags, *tdyn, xparams, actdis, solver, output);
      break;
    }
    default:
      FOUR_C_THROW("unknown time integration scheme '%s'", timinttype);
      break;
  }

  // link/store thermal field solver
  thermo_ = tmpthr;

  // see you
  return;
}  // setup_tim_int()


/*----------------------------------------------------------------------*
 | integrate                                                bborn 08/09 |
 *----------------------------------------------------------------------*/
void Thermo::Adapter::integrate()
{
  while (not_finished())
  {
    // call the predictor
    prepare_time_step();

    // integrate time step
    Inpar::Thermo::ConvergenceStatus convStatus = solve();

    switch (convStatus)
    {
      case Inpar::Thermo::conv_success:
        update();
        print_step();
        output();
        break;
      case Inpar::Thermo::conv_fail_repeat:
        // do not update and output but try again
        continue;
      default:
        // no other convergence status can be handled at this point, abort
        FOUR_C_THROW("Solver failed.");
    }
  }
  // print monitoring of time consumption
  Teuchos::TimeMonitor::summarize();
}  // Integrate()


/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
