/*----------------------------------------------------------------------*/
/*! \file

\brief Thermo field adapter
\level 1


*/


/*----------------------------------------------------------------------*
 | headers                                                  bborn 08/09 |
 *----------------------------------------------------------------------*/
#include "4C_adapter_thermo.hpp"

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
Adapter::ThermoBaseAlgorithm::ThermoBaseAlgorithm(
    const Teuchos::ParameterList& prbdyn, Teuchos::RCP<Core::FE::Discretization> actdis)
{
  setup_thermo(prbdyn, actdis);
}



/*----------------------------------------------------------------------*
 |                                                          bborn 08/09 |
 *----------------------------------------------------------------------*/
void Adapter::ThermoBaseAlgorithm::setup_thermo(
    const Teuchos::ParameterList& prbdyn, Teuchos::RCP<Core::FE::Discretization> actdis)
{
  const Teuchos::ParameterList& tdyn = Global::Problem::Instance()->thermal_dynamic_params();

  // major switch to different time integrators
  Inpar::THR::DynamicType timinttype =
      Core::UTILS::IntegralValue<Inpar::THR::DynamicType>(tdyn, "DYNAMICTYP");
  switch (timinttype)
  {
    case Inpar::THR::dyna_statics:
    case Inpar::THR::dyna_onesteptheta:
    case Inpar::THR::dyna_genalpha:
    case Inpar::THR::dyna_expleuler:
      setup_tim_int(prbdyn, timinttype, actdis);  // <-- here is the show
      break;
    default:
      FOUR_C_THROW(
          "unknown time integration scheme '%s'", tdyn.get<std::string>("DYNAMICTYP").c_str());
      break;
  }

}  // setup_thermo()


/*----------------------------------------------------------------------*
 | setup of thermal time integration                        bborn 08/09 |
 *----------------------------------------------------------------------*/
void Adapter::ThermoBaseAlgorithm::setup_tim_int(const Teuchos::ParameterList& prbdyn,
    Inpar::THR::DynamicType timinttype, Teuchos::RCP<Core::FE::Discretization> actdis)
{
  // this is not exactly a one hundred meter race, but we need timing
  Teuchos::RCP<Teuchos::Time> t =
      Teuchos::TimeMonitor::getNewTimer("Adapter::ThermoBaseAlgorithm::setup_thermo");
  Teuchos::TimeMonitor monitor(*t);

  // set degrees of freedom in the discretization
  if (not actdis->Filled()) actdis->fill_complete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<Core::IO::DiscretizationWriter> output = actdis->Writer();
  output->WriteMesh(0, 0.0);

  //  // get input parameter lists and copy them, because a few parameters are overwritten
  const Teuchos::RCP<Teuchos::ParameterList> ioflags =
      Teuchos::rcp(new Teuchos::ParameterList(Global::Problem::Instance()->IOParams()));
  const Teuchos::RCP<Teuchos::ParameterList> tdyn = Teuchos::rcp(
      new Teuchos::ParameterList(Global::Problem::Instance()->thermal_dynamic_params()));
  //  //const Teuchos::ParameterList& size
  //  //  = Global::Problem::Instance()->ProblemSizeParams();

  // show default parameters of thermo parameter list
  if ((actdis->Comm()).MyPID() == 0) Input::PrintDefaultParameters(Core::IO::cout, *tdyn);

  // add extra parameters (a kind of work-around)
  Teuchos::RCP<Teuchos::ParameterList> xparams = Teuchos::rcp(new Teuchos::ParameterList());

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
  Teuchos::RCP<Teuchos::ParameterList> solveparams = Teuchos::rcp(new Teuchos::ParameterList());
  Teuchos::RCP<Core::LinAlg::Solver> solver = Teuchos::rcp(
      new Core::LinAlg::Solver(Global::Problem::Instance()->SolverParams(linsolvernumber),
          actdis->Comm(), Global::Problem::Instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::Instance()->IOParams(), "VERBOSITY")));
  actdis->compute_null_space_if_necessary(solver->Params());

  // create marching time integrator
  Teuchos::RCP<Thermo> tmpthr;
  switch (timinttype)
  {
    case Inpar::THR::dyna_statics:
    {
      tmpthr =
          Teuchos::rcp(new THR::TimIntStatics(*ioflags, *tdyn, *xparams, actdis, solver, output));
      break;
    }
    case Inpar::THR::dyna_onesteptheta:
    {
      tmpthr = Teuchos::rcp(
          new THR::TimIntOneStepTheta(*ioflags, *tdyn, *xparams, actdis, solver, output));
      break;
    }
    case Inpar::THR::dyna_genalpha:
    {
      tmpthr =
          Teuchos::rcp(new THR::TimIntGenAlpha(*ioflags, *tdyn, *xparams, actdis, solver, output));
      break;
    }
    case Inpar::THR::dyna_expleuler:
    {
      tmpthr =
          Teuchos::rcp(new THR::TimIntExplEuler(*ioflags, *tdyn, *xparams, actdis, solver, output));
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
void Adapter::Thermo::Integrate()
{
  while (NotFinished())
  {
    // call the predictor
    prepare_time_step();

    // integrate time step
    Inpar::THR::ConvergenceStatus convStatus = Solve();

    switch (convStatus)
    {
      case Inpar::THR::conv_success:
        Update();
        PrintStep();
        Output();
        break;
      case Inpar::THR::conv_fail_repeat:
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
