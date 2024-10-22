// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_particle_algorithm_gravity.hpp"

#include "4C_global_data.hpp"
#include "4C_utils_function_of_time.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::GravityHandler::GravityHandler(const Teuchos::ParameterList& params)
    : params_(params), gravityrampfctnumber_(params.get<int>("GRAVITY_RAMP_FUNCT"))
{
  // empty constructor
}

void PARTICLEALGORITHM::GravityHandler::init(const std::vector<double>& gravity)
{
  // set gravity acceleration vector
  gravity_ = gravity;

  // safety check
  if (static_cast<int>(gravity_.size()) != 3)
    FOUR_C_THROW("dimension (dim = %d) of gravity acceleration vector is wrong!",
        static_cast<int>(gravity_.size()));
}

void PARTICLEALGORITHM::GravityHandler::setup()
{
  // nothing to do
}

void PARTICLEALGORITHM::GravityHandler::get_gravity_acceleration(
    const double time, std::vector<double>& scaled_gravity)
{
  scaled_gravity = gravity_;

  // evaluate gravity ramp function
  if (gravityrampfctnumber_ > 0)
  {
    const double fac = Global::Problem::instance()
                           ->function_by_id<Core::Utils::FunctionOfTime>(gravityrampfctnumber_ - 1)
                           .evaluate(time);

    for (int dim = 0; dim < 3; ++dim) scaled_gravity[dim] *= fac;
  }
}

FOUR_C_NAMESPACE_CLOSE
