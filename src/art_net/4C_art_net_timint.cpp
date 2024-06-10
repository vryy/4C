/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for arterial network (time) integration.


\level 3

*----------------------------------------------------------------------*/


#include "4C_art_net_timint.hpp"

#include "4C_art_net_artery_resulttest.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_bio.hpp"
#include "4C_linear_solver_method_linalg.hpp"

FOUR_C_NAMESPACE_OPEN


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                kremheller 03/18|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

Arteries::TimInt::TimInt(Teuchos::RCP<Core::FE::Discretization> actdis, const int linsolvernumber,
    const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& artparams,
    Core::IO::DiscretizationWriter& output)
    : discret_(actdis),
      solver_(Teuchos::null),
      params_(probparams),
      output_(output),
      dtele_(0.0),
      dtsolve_(0.0),
      time_(0.0),
      step_(0),
      uprestart_(params_.get<int>("RESTARTEVRY")),
      upres_(params_.get<int>("RESULTSEVRY")),
      linsolvernumber_(linsolvernumber),
      coupledTo3D_(false)
{
  // -------------------------------------------------------------------
  // get the processor ID from the communicator
  // -------------------------------------------------------------------
  myrank_ = discret_->Comm().MyPID();

  // -------------------------------------------------------------------
  // get the basic parameters first
  // -------------------------------------------------------------------
  // time-step size
  dtp_ = dta_ = params_.get<double>("TIMESTEP");
  // maximum number of timesteps
  stepmax_ = params_.get<int>("NUMSTEP");
  // maximum simulation time
  maxtime_ = dtp_ * double(stepmax_);

  // solve scatra flag
  solvescatra_ = Core::UTILS::IntegralValue<int>(artparams, "SOLVESCATRA");

  if (linsolvernumber_ == -1) FOUR_C_THROW("Set a valid linear solver for arterial network");
}



/*------------------------------------------------------------------------*
 | initialize time integration                            kremheller 03/18 |
 *------------------------------------------------------------------------*/
void Arteries::TimInt::Init(const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& arteryparams, const std::string& scatra_disname)
{
  solver_ = Teuchos::rcp(new Core::LinAlg::Solver(
      Global::Problem::Instance()->SolverParams(linsolvernumber_), discret_->Comm()));

  discret_->compute_null_space_if_necessary(solver_->Params());

  return;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | Start the time integration.                                          |
 |                                                          ismail 06/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void Arteries::TimInt::Integrate(
    bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingParams)
{
  coupledTo3D_ = CoupledTo3D;
  if (CoupledTo3D && CouplingParams.get() == nullptr)
  {
    FOUR_C_THROW(
        "Coupling parameter list is not allowed to be empty, If a 3-D/reduced-D coupling is "
        "defined\n");
  }
  // Prepare the loop
  prepare_time_loop();

  TimeLoop(CoupledTo3D, CouplingParams);

  // print the results of time measurements
  if (!coupledTo3D_)
  {
    Teuchos::TimeMonitor::summarize();
  }

  return;
}  // ArtNetExplicitTimeInt::Integrate

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | prepare the time loop                                kremheller 03/18|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void Arteries::TimInt::prepare_time_loop()
{
  // do nothing
  return;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | contains the time loop                                   ismail 06/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void Arteries::TimInt::TimeLoop(
    bool CoupledTo3D, Teuchos::RCP<Teuchos::ParameterList> CouplingTo3DParams)
{
  coupledTo3D_ = CoupledTo3D;
  // time measurement: time loop
  if (!coupledTo3D_)
  {
    TEUCHOS_FUNC_TIME_MONITOR(" + time loop");
  }

  while (step_ < stepmax_ and time_ < maxtime_)
  {
    prepare_time_step();
    // -------------------------------------------------------------------
    //                       output to screen
    // -------------------------------------------------------------------
    if (myrank_ == 0)
    {
      if (!coupledTo3D_)
      {
        printf("TIME: %11.4E/%11.4E  DT = %11.4E   Solving Artery    STEP = %4d/%4d \n", time_,
            maxtime_, dta_, step_, stepmax_);
      }
      else
      {
        printf(
            "SUBSCALE_TIME: %11.4E/%11.4E  SUBSCALE_DT = %11.4E   Solving Artery    SUBSCALE_STEP "
            "= %4d/%4d \n",
            time_, maxtime_, dta_, step_, stepmax_);
      }
    }

    Solve(CouplingTo3DParams);

    // -------------------------------------------------------------------
    //                         Solve for the scalar transport
    // -------------------------------------------------------------------
    if (solvescatra_)
    {
      SolveScatra();
    }


    // -------------------------------------------------------------------
    //                         update solution
    //        current solution becomes old solution of next timestep
    // -------------------------------------------------------------------
    TimeUpdate();

    // -------------------------------------------------------------------
    //  lift'n'drag forces, statistics time sample and output of solution
    //  and statistics
    // -------------------------------------------------------------------
    if (!CoupledTo3D)
    {
      Output(CoupledTo3D, CouplingTo3DParams);
    }

    // -------------------------------------------------------------------
    //                       update time step sizes
    // -------------------------------------------------------------------
    dtp_ = dta_;

    // -------------------------------------------------------------------
    //                    stop criterium for timeloop
    // -------------------------------------------------------------------
    if (CoupledTo3D)
    {
      break;
    }
  }

}  // Arteries::TimInt::TimeLoop


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 | setup the variables to do a new time step                ismail 06/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void Arteries::TimInt::prepare_time_step()
{
  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  step_ += 1;
  time_ += dta_;
}

FOUR_C_NAMESPACE_CLOSE
