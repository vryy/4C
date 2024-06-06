/*-----------------------------------------------------------*/
/*! \file

\brief Setting of time-integration parameters in fluid element evaluation

This class provides the fluid element parameter for the time integration,
which are unique and equal for every fluid in the problem. Time integration
with different parameters in more than one fluid field is not yet supported.


\level 1

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_ele_parameter_timint.hpp"

#include "4C_io_pstream.hpp"

#include <iostream>
#include <string>

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::FluidEleParameterTimInt* Discret::ELEMENTS::FluidEleParameterTimInt::Instance(
    Core::UTILS::SingletonAction action)
{
  static auto singleton_owner = Core::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<Discret::ELEMENTS::FluidEleParameterTimInt>(
            new Discret::ELEMENTS::FluidEleParameterTimInt());
      });

  return singleton_owner.Instance(action);
}

//----------------------------------------------------------------------*/
// private constructor of FluidEleParameterTimInt
//----------------------------------------------------------------------*/
Discret::ELEMENTS::FluidEleParameterTimInt::FluidEleParameterTimInt()
    : set_general_fluid_timeparameter_(false),
      is_genalpha_(false),
      is_genalpha_np_(false),
      is_stationary_(false),
      is_one_step_theta_(false),
      is_cont_impl_press_impl_(false),
      is_cont_impl_press_normal_(false),
      ostnew_(false),
      time_(-1.0),
      dt_(0.0),
      timefac_(0.0),
      theta_(0.0),
      omtheta_(0.0),
      gamma_(0.0),
      alpha_f_(0.0),
      alpha_m_(0.0),
      afgdt_(1.0),
      timefacrhs_(1.0),
      timefacpre_(1.0)
{
}

//----------------------------------------------------------------------*/
// set time parameters which are equal for every fluid  rasthofer 11/13 |
//----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidEleParameterTimInt::set_element_time_parameter(
    Teuchos::ParameterList& params)
{
  // second check: timealgo
  // work around to use SetTimeParameter in GenAlpha (Neumann BC)
  if (set_general_fluid_timeparameter_ == false)
  {
    set_general_fluid_timeparameter_ = true;
  }

  // set flag, time integration scheme
  timealgo_ =
      Core::UTILS::GetAsEnum<Inpar::FLUID::TimeIntegrationScheme>(params, "TimeIntegrationScheme");

  // set time integration scheme-specific element parameters
  if (timealgo_ == Inpar::FLUID::timeint_stationary)
  {
    is_genalpha_ = false;
    is_stationary_ = true;
    is_genalpha_np_ = false;
    is_one_step_theta_ = false;
  }
  else if (timealgo_ == Inpar::FLUID::timeint_afgenalpha)
  {
    is_genalpha_ = true;
    is_stationary_ = false;
    is_genalpha_np_ = false;
    is_one_step_theta_ = false;
  }
  else if (timealgo_ == Inpar::FLUID::timeint_npgenalpha)
  {
    is_genalpha_ = true;
    is_stationary_ = false;
    is_genalpha_np_ = true;
    is_one_step_theta_ = false;
  }
  else
  {
    is_genalpha_ = false;
    is_stationary_ = false;
    is_genalpha_np_ = false;
    is_one_step_theta_ = true;
  }


  //----------------------------------------------------------------------
  // get control parameters for time integration
  //----------------------------------------------------------------------

  // get current time: n+alpha_F for generalized-alpha scheme, n+1 otherwise
  time_ = params.get<double>("total time", -1.0);

  // set global variable timefac to zero
  timefac_ = 0.0;

  if (not is_stationary_)
  {
    // get time-step length and time-integration parameters
    dt_ = params.get<double>("dt", -1.0);
    theta_ = params.get<double>("theta", -1.0);
    omtheta_ = params.get<double>("omtheta", -1.0);

    // compute timefactor for left-hand side:
    // one-step-Theta:    timefac = theta*dt
    // BDF2:              timefac = 2/3 * dt
    // generalized-alpha: timefac = (alpha_F/alpha_M) * gamma * dt
    // (For BDF2 and generalized-alpha, theta was already computed
    //  accordingly in FLD::FluidImplicitTimeInt::prepare_time_step().)

    //-----------------------------------------------------------------------
    //       |          timefac         |  timefacpre     |    timefacrhs   |
    // ----------------------------------------------------------------------
    // OST   |                        dt*theta                              |
    //-----------------------------------------------------------------------
    // BDF2  |                        2/3 * dt                              |
    //-----------------------------------------------------------------------
    // Af GA |          alphaF*gamma*dt/alphaM            | gamma*dt/alphaM |
    //----------------------------------------------------------------------
    // NP GA | alphaF*gamma*dt/alphaM   | gamma*dt/alphaM | gamma*dt/alphaM |
    //-----------------------------------------------------------------------

    timefac_ = theta_ * dt_;

    // compute generalized-alpha-related values and set them appropriately
    // otherwise
    if (is_genalpha_)
    {
      gamma_ = params.get<double>("gamma", -1.0);
      alpha_f_ = params.get<double>("alphaF", -1.0);
      alpha_m_ = params.get<double>("alphaM", -1.0);
    }
    else
    {
      gamma_ = theta_;
      alpha_f_ = 1.0;
      alpha_m_ = 1.0;
    }

    // if not generalized-alpha: afgdt = theta * dt_ = timefac_
    // Peter's generalized alpha: timefacmat_u_ for velocity terms
    afgdt_ = alpha_f_ * gamma_ * dt_;

    // timeint_gen_alpha = p(n+1) (Peter's genalpha)
    if (timealgo_ == Inpar::FLUID::timeint_npgenalpha)
    {
      // if not generalized-alpha: timefacrhs_=theta * dt_ = timefac_
      timefacpre_ = gamma_ / alpha_m_ * dt_;
      timefacrhs_ = gamma_ / alpha_m_ * dt_;
    }
    else if (timealgo_ == Inpar::FLUID::timeint_afgenalpha)
    {
      timefacpre_ = gamma_ * alpha_f_ / alpha_m_ * dt_;
      timefacrhs_ = gamma_ / alpha_m_ * dt_;
    }
    else
    {
      // if not generalized-alpha: timefacmat_p_=theta * dt_ = timefac_
      timefacpre_ = gamma_ * alpha_f_ / alpha_m_ * dt_;
      // if not generalized-alpha: timefacrhs_=theta * dt_ = timefac_
      timefacrhs_ = gamma_ * alpha_f_ / alpha_m_ * dt_;

      // set flag, time integration scheme
      ostalgo_ =
          Core::UTILS::GetAsEnum<Inpar::FLUID::OstContAndPress>(params, "ost cont and press");
      ostnew_ = params.get<bool>("ost new", false);

      if (ostnew_)
      {
        // set time integration scheme-specific element parameters
        if (ostalgo_ == Inpar::FLUID::Cont_impl_Press_impl)
        {
          is_cont_impl_press_impl_ = true;
          is_cont_impl_press_normal_ = false;
        }
        else if (ostalgo_ == Inpar::FLUID::Cont_impl_Press_normal)
        {
          is_cont_impl_press_impl_ = false;
          is_cont_impl_press_normal_ = true;
        }
        else
        {
          is_cont_impl_press_impl_ = false;
          is_cont_impl_press_normal_ = false;
        }
      }
    }
  }
  else  // is_stationary == true
  {
    // set timefactor for stationary case to 1.0
    timefac_ = 1.0;
    timefacrhs_ = 1.0;
  }

  if (dt_ < 0.0 or theta_ < 0.0 or time_ < 0.0 or omtheta_ < 0.0 or gamma_ < 0.0 or
      alpha_f_ < 0.0 or alpha_m_ < 0.0)
  {
    std::cout << "dt_: " << dt_ << std::endl;
    std::cout << "theta_ " << theta_ << std::endl;
    std::cout << "time_ " << time_ << std::endl;
    std::cout << "omtheta_ " << omtheta_ << std::endl;
    std::cout << "gamma_ " << gamma_ << std::endl;
    std::cout << "alphaF_ " << alpha_f_ << std::endl;
    std::cout << "alphaM_ " << alpha_m_ << std::endl;
    FOUR_C_THROW("Negative (or no) time-integration parameter or time-step length supplied");
  }
}

//----------------------------------------------------------------------*/
// print fluid time parameter to screen                 rasthofer 11/13 |
//----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidEleParameterTimInt::print_fluid_time_parameter()
{
  std::cout << std::endl
            << "|-----------------------------------------------------------" << std::endl;
  std::cout << "|  Fluid Element Time parameter: " << std::endl;
  std::cout << "|-------------------------------------------------------------------" << std::endl;
  //! time parameters set?
  std::cout << "|    time algorithm:         " << timealgo_ << std::endl;
  //! is stationary?
  std::cout << "|    is stationary?:         " << is_stationary_ << std::endl;
  //! time algorithm
  std::cout << "|    time parameters set:    " << set_general_fluid_timeparameter_ << std::endl;
  //! actual time to evaluate the body BC
  std::cout << "|    time:                   " << time_ << std::endl;
  //! time-step length
  std::cout << "|    time step:              " << dt_ << std::endl;
  //! timefac = dt_ * ("pseudo"-)theta_
  std::cout << "|    time factor:            " << timefac_ << std::endl;
  //! factor for left-hand side due to one-step-theta time-integration scheme
  std::cout << "|    theta:                  " << theta_ << std::endl;
  //! factor for right-hand side due to one-step-theta time-integration scheme
  std::cout << "|    (1-theta):              " << omtheta_ << std::endl;
  //! generalised-alpha parameter (connecting velocity and acceleration)
  std::cout << "|    gamma:                  " << gamma_ << std::endl;
  //! generalised-alpha parameter (velocity)
  std::cout << "|    alpha_F:                " << alpha_f_ << std::endl;
  //! generalised-alpha parameter (acceleration)
  std::cout << "|    alpha_M:                " << alpha_m_ << std::endl;
  //! generalised-alpha parameter, alphaF_*gamma_*dt_
  std::cout << "|    time factor mat_u:      " << afgdt_ << std::endl;
  //! time integration factor for the right hand side (boundary elements)
  std::cout << "|    time factor rhs:        " << timefacrhs_ << std::endl;
  //! time integration factor for the left hand side (pressure)
  std::cout << "|    time factor mat_p:      " << timefacpre_ << std::endl;
  std::cout << std::endl
            << "|-----------------------------------------------------------" << std::endl;
}

FOUR_C_NAMESPACE_CLOSE
