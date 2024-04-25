/*----------------------------------------------------------------------*/
/*! \file
 \brief one-step-theta time integration scheme for porous multiphase flow

   \level 3

 *----------------------------------------------------------------------*/


#include "4C_porofluidmultiphase_timint_ost.hpp"

#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_porofluidmultiphase_ele_action.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Constructor (public)                                   vuong  08/16 |
 *----------------------------------------------------------------------*/
POROFLUIDMULTIPHASE::TimIntOneStepTheta::TimIntOneStepTheta(
    Teuchos::RCP<DRT::Discretization> dis,  //!< discretization
    const int linsolvernumber,              //!< number of linear solver
    const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& poroparams,
    Teuchos::RCP<IO::DiscretizationWriter> output  //!< output writer
    )
    : TimIntImpl(dis, linsolvernumber, probparams, poroparams, output),
      theta_(poroparams.get<double>("THETA"))
{
}



/*----------------------------------------------------------------------*
 |  set parameter for element evaluation                    vuong 06/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntOneStepTheta::SetElementTimeStepParameter() const
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action", POROFLUIDMULTIPHASE::set_timestep_parameter);

  // the total time definitely changes
  eleparams.set<double>("total time", time_);
  // we set the time step and related, just in case we want adaptive time stepping
  eleparams.set<double>("time-step length", dt_);
  eleparams.set<double>("time factor", theta_ * dt_);

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
}


/*-----------------------------------------------------------------------------*
 | set time for evaluation of POINT -Neumann boundary conditions   vuong 08/16 |
 *----------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntOneStepTheta::SetTimeForNeumannEvaluation(
    Teuchos::ParameterList& params)
{
  params.set("total time", time_);
}


/*----------------------------------------------------------------------*
| Print information about current time step to screen      vuong 08/16  |
*-----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntOneStepTheta::PrintTimeStepInfo()
{
  if (myrank_ == 0)
  {
    std::cout << "| TIME: " << std::setw(11) << std::setprecision(4) << std::scientific << time_
              << "/" << std::setw(11) << std::setprecision(4) << std::scientific << maxtime_
              << "  DT = " << std::setw(11) << std::setprecision(4) << std::scientific << dt_
              << "  "
              << "One-Step-Theta (theta = " << std::setw(3) << std::setprecision(2) << theta_
              << ") STEP = " << std::setw(4) << step_ << "/" << std::setw(4) << stepmax_
              << "            |" << '\n';
  }
}


/*----------------------------------------------------------------------*
 | set part of the residual vector belonging to the old timestep        |
 |                                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntOneStepTheta::SetOldPartOfRighthandside()
{
  // hist_ = phin_ + dt*(1-Theta)*phidtn_
  hist_->Update(1.0, *phin_, dt_ * (1.0 - theta_), *phidtn_, 0.0);
}


/*----------------------------------------------------------------------*
 | perform an explicit predictor step                       vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntOneStepTheta::ExplicitPredictor()
{
  phinp_->Update(dt_, *phidtn_, 1.0);
}


/*----------------------------------------------------------------------*
 | add actual Neumann loads                                             |
 | scaled with a factor resulting from time discretization  vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntOneStepTheta::AddNeumannToResidual()
{
  residual_->Update(theta_ * dt_, *neumann_loads_, 1.0);
}


/*----------------------------------------------------------------------------*
 | add global state vectors specific for time-integration scheme  vuong 08/16 |
 *---------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntOneStepTheta::AddTimeIntegrationSpecificVectors()
{
  discret_->SetState("hist", hist_);
  discret_->SetState("phinp_fluid", phinp_);
  discret_->SetState("phin_fluid", phin_);
  discret_->SetState("phidtnp", phidtnp_);
}


/*----------------------------------------------------------------------*
 | compute time derivative                                  vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntOneStepTheta::ComputeTimeDerivative()
{
  // time derivative of phi:
  // phidt(n+1) = (phi(n+1)-phi(n)) / (theta*dt) + (1-(1/theta))*phidt(n)
  const double fact1 = 1.0 / (theta_ * dt_);
  const double fact2 = 1.0 - (1.0 / theta_);
  phidtnp_->Update(fact2, *phidtn_, 0.0);
  phidtnp_->Update(fact1, *phinp_, -fact1, *phin_, 1.0);
}


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntOneStepTheta::Update()
{
  // call base class
  POROFLUIDMULTIPHASE::TimIntImpl::Update();

  // compute time derivative at time n+1
  ComputeTimeDerivative();

  // solution of this step becomes most recent solution of the last step
  phin_->Update(1.0, *phinp_, 0.0);

  // time deriv. of this step becomes most recent time derivative of
  // last step
  phidtn_->Update(1.0, *phidtnp_, 0.0);
}


/*----------------------------------------------------------------------*
 | write additional data required for restart               vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntOneStepTheta::OutputRestart()
{
  // additional state vectors that are needed for One-Step-Theta restart
  output_->WriteVector("phidtn_fluid", phidtn_);
  output_->WriteVector("phin_fluid", phin_);
}


/*----------------------------------------------------------------------*
 |  read restart data                                       vuong 08/16 |
 -----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntOneStepTheta::ReadRestart(const int step)
{
  // call base class
  POROFLUIDMULTIPHASE::TimIntImpl::ReadRestart(step);

  Teuchos::RCP<IO::DiscretizationReader> reader(Teuchos::null);
  reader = Teuchos::rcp(new IO::DiscretizationReader(
      discret_, GLOBAL::Problem::Instance()->InputControlFile(), step));

  time_ = reader->ReadDouble("time");
  step_ = reader->ReadInt("step");

  if (myrank_ == 0)
    std::cout << "Reading POROFLUIDMULTIPHASE restart data (time=" << time_ << " ; step=" << step_
              << ")" << '\n';

  // read state vectors that are needed for One-Step-Theta restart
  reader->ReadVector(phinp_, "phinp_fluid");
  reader->ReadVector(phin_, "phin_fluid");
  reader->ReadVector(phidtn_, "phidtn_fluid");
}

/*--------------------------------------------------------------------*
 | calculate init time derivatives of state variables kremheller 03/17 |
 *--------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntOneStepTheta::CalcInitialTimeDerivative()
{
  // standard general element parameter without stabilization
  SetElementGeneralParameters();

  // we also have to modify the time-parameter list (incremental solve)
  // actually we do not need a time integration scheme for calculating the initial time derivatives,
  // but the rhs of the standard element routine is used as starting point for this special system
  // of equations. Therefore, the rhs vector has to be scaled correctly.
  SetElementTimeStepParameter();

  // call core algorithm
  TimIntImpl::CalcInitialTimeDerivative();

  // and finally undo our temporary settings
  SetElementGeneralParameters();
  SetElementTimeStepParameter();
}

FOUR_C_NAMESPACE_CLOSE
