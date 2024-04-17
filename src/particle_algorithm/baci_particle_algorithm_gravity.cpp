/*---------------------------------------------------------------------------*/
/*! \file
\brief gravity handler for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_algorithm_gravity.hpp"

#include "baci_global_data.hpp"
#include "baci_utils_function_of_time.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::GravityHandler::GravityHandler(const Teuchos::ParameterList& params)
    : params_(params), gravityrampfctnumber_(params.get<int>("GRAVITY_RAMP_FUNCT"))
{
  // empty constructor
}

void PARTICLEALGORITHM::GravityHandler::Init(const std::vector<double>& gravity)
{
  // set gravity acceleration vector
  gravity_ = gravity;

  // safety check
  if (static_cast<int>(gravity_.size()) != 3)
    dserror("dimension (dim = %d) of gravity acceleration vector is wrong!",
        static_cast<int>(gravity_.size()));
}

void PARTICLEALGORITHM::GravityHandler::Setup()
{
  // nothing to do
}

void PARTICLEALGORITHM::GravityHandler::GetGravityAcceleration(
    const double time, std::vector<double>& scaled_gravity)
{
  scaled_gravity = gravity_;

  // evaluate gravity ramp function
  if (gravityrampfctnumber_ > 0)
  {
    const double fac = GLOBAL::Problem::Instance()
                           ->FunctionById<CORE::UTILS::FunctionOfTime>(gravityrampfctnumber_ - 1)
                           .Evaluate(time);

    for (int dim = 0; dim < 3; ++dim) scaled_gravity[dim] *= fac;
  }
}

FOUR_C_NAMESPACE_CLOSE
