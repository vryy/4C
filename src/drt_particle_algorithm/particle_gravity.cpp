/*---------------------------------------------------------------------------*/
/*!
\file particle_gravity.cpp

\brief gravity handler for particle simulations

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_gravity.H"

#include "../drt_lib/drt_globalproblem.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::GravityHandler::GravityHandler(const Teuchos::ParameterList& params)
    : params_(params), gravityrampfctnumber_(params.get<int>("GRAVITY_RAMP_FUNCT"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init gravity handler                                       sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::GravityHandler::Init(const std::vector<double>& gravity)
{
  // set gravity acceleration vector
  gravity_ = gravity;

  // safety check
  if ((int)gravity_.size() != 3)
    dserror("dimension (dim = %d) of gravity acceleration vector is wrong!", (int)gravity_.size());
}

/*---------------------------------------------------------------------------*
 | setup gravity handler                                      sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::GravityHandler::Setup()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | write restart of gravity handler                           sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::GravityHandler::WriteRestart(const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of gravity handler                            sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::GravityHandler::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | get gravity acceleration                                   sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::GravityHandler::GetGravityAcceleration(
    const double time, std::vector<double>& scaled_gravity)
{
  scaled_gravity = gravity_;

  // evaluate gravity ramp function
  if (gravityrampfctnumber_ > 0)
  {
    const double fac =
        DRT::Problem::Instance()->Funct(gravityrampfctnumber_ - 1).EvaluateTime(time);

    for (int dim = 0; dim < 3; ++dim) scaled_gravity[dim] *= fac;
  }
}
