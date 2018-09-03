/*!------------------------------------------------------------------------------------------------*
\file ad_opt_fluid_adjoint_impl.cpp

\brief adapter for element routines of fluid adjoint equations in topology optimization

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "topopt_fluidAdjoint_timeint.H"
#include "../drt_lib/drt_discret.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
TOPOPT::ADJOINT::FluidAdjointTimeInt::FluidAdjointTimeInt(Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<IO::DiscretizationWriter> output)
    : discret_(discret),
      solver_(solver),
      params_(params),
      output_(output),
      myrank_(discret_->Comm().MyPID()),
      uprestart_(params_->get<int>("write restart every", -1)),
      upres_(params_->get<int>("write solution every", -1)),
      numdim_(params_->get<int>("number of velocity degrees of freedom")),
      dt_(params_->get<double>("time step size")),
      stepmax_(params_->get<int>("max number timesteps")),
      maxtime_(params_->get<double>("total time"))
{
  timealgo_ = DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params_, "time int algo");
  adjointtype_ = params_->get<INPAR::TOPOPT::AdjointType>("adjoint type");

  // check setting of time parameter
  if (fabs(maxtime_ - dt_ * stepmax_) > 1.0e-14)
  {
    dserror(
        "Fix total simulation (T = %f, dt = %f, n = %i\n"
        "so that: T = dt * n",
        maxtime_, dt_, stepmax_);
  }

  // set initial time = endtime so that it fits to the fluid parameter setting
  // potentially we do one step less here than in fluid
  // fluid criteria are:
  // o endtime <= numstep * dt
  // o endtime < maxtime + dt
  //
  // additionally evaluate the correct number of time steps since it is
  // required for evaluation of the fluid velocity at the correct step
  if (timealgo_ == INPAR::FLUID::timeint_stationary)
  {
    time_ = maxtime_ + dt_;
    step_ = 0;
  }
  else
  {
    // for discrete adjoints start with step 0 and end with step stepmax-1
    if (adjointtype_ == INPAR::TOPOPT::discrete_adjoint)
    {
      time_ = maxtime_ + dt_;
      step_ = -1;
    }
    else
    {
      time_ = maxtime_;
      step_ = 0;
    }
  }
}



bool TOPOPT::ADJOINT::FluidAdjointTimeInt::TimeLoopFinished() const
{
  // for discrete adjoints we start with step 0 and end with step stepmax-1
  // since the adjoint solution of step stepmax at time 0 does not influence
  // the optimization since the primal solution u_0 is independent of the
  // optimization variable
  if ((adjointtype_ == INPAR::TOPOPT::discrete_adjoint) and (step_ == stepmax_ - 1))
    return true;
  else if (step_ == stepmax_)
    return true;

  return false;
}
