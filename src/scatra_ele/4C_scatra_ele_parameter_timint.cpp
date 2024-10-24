// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_scatra_ele_parameter_timint.hpp"

#include "4C_utils_exceptions.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::Elements::ScaTraEleParameterTimInt* Discret::Elements::ScaTraEleParameterTimInt::instance(
    const std::string& disname)
{
  static auto singleton_map =
      Core::Utils::make_singleton_map<std::string>([](const std::string& disname)
          { return std::unique_ptr<ScaTraEleParameterTimInt>(new ScaTraEleParameterTimInt); });

  return singleton_map[disname].instance(Core::Utils::SingletonAction::create, disname);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::Elements::ScaTraEleParameterTimInt::ScaTraEleParameterTimInt()
    : is_genalpha_(false),
      is_stationary_(false),
      is_incremental_(false),
      time_(-1.0),
      timederivativefac_(-1.0),
      dt_(0.0),
      timefac_(0.0),
      timefacrhs_(0.0),
      timefacrhstau_(0.0),
      alpha_f_(0.0)
{
}

//----------------------------------------------------------------------*/
//----------------------------------------------------------------------*/
void Discret::Elements::ScaTraEleParameterTimInt::set_parameters(Teuchos::ParameterList& parameters)
{
  // get control parameters
  is_stationary_ = parameters.get<bool>("using stationary formulation");
  is_genalpha_ = parameters.get<bool>("using generalized-alpha time integration");
  is_incremental_ = parameters.get<bool>("incremental solver");

  // get current time and time-step length
  time_ = parameters.get<double>("total time");
  dt_ = parameters.get<double>("time-step length");

  // get time factor and alpha_F if required
  // one-step-Theta:    timefac = theta*dt
  // BDF2:              timefac = 2/3 * dt
  // generalized-alpha: timefac = alphaF * (gamma/alpha_M) * dt

  //-----------------------------------------------------
  //       |          timefac         |    timefacrhs   |
  // ----------------------------------------------------
  // OST   |                  dt*theta                  |
  //-----------------------------------------------------
  // BDF2  |               2/3 * dt                     |
  //-----------------------------------------------------
  // Af GA | alphaF*gamma*dt/alphaM   | gamma*dt/alphaM |
  //-----------------------------------------------------

  timefac_ = 1.0;
  alpha_f_ = 1.0;
  timefacrhs_ = 1.0;
  timefacrhstau_ = 1.0;

  if (not is_stationary_)
  {
    timefac_ = parameters.get<double>("time factor");

    if (is_genalpha_)
    {
      alpha_f_ = parameters.get<double>("alpha_F");
      timefac_ *= alpha_f_;
    }
    if (timefac_ < 0.0) FOUR_C_THROW("time factor is negative.");
  }

  if (not is_stationary_)
  {
    if (is_genalpha_)
    {
      timefacrhs_ = timefac_ / alpha_f_;
      timefacrhstau_ = timefacrhs_;
      if (not is_incremental_) timefacrhs_ *= (1.0 - alpha_f_);
    }
    else
    {
      if (not is_incremental_)
        timefacrhs_ = 0.0;
      else
        timefacrhs_ = timefac_;
    }
  }
  else
  {
    if (not is_incremental_) timefacrhs_ = 0.0;
  }

  timederivativefac_ = parameters.get<double>("time derivative factor", -1.0);
}

FOUR_C_NAMESPACE_CLOSE
