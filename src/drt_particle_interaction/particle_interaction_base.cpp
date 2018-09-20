/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_base.cpp

\brief base particle interaction handler

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
#include "particle_interaction_base.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::ParticleInteractionBase::ParticleInteractionBase(
    const Epetra_Comm& comm, const Teuchos::ParameterList& params)
    : comm_(comm), myrank_(comm.MyPID()), params_(params), time_(0.0), dt_(0.0)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init particle interaction handler                          sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionBase::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup particle interaction handler                         sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionBase::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // init vector
  gravity_.resize(3, 0.0);
}

/*---------------------------------------------------------------------------*
 | write restart of particle interaction handler              sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionBase::WriteRestart(
    const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of particle interaction handler               sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | set current time                                           sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionBase::SetCurrentTime(const double currenttime)
{
  time_ = currenttime;
}

/*---------------------------------------------------------------------------*
 | set current step size                                      sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionBase::SetCurrentStepSize(const double currentstepsize)
{
  dt_ = currentstepsize;
}

/*---------------------------------------------------------------------------*
 | set gravity                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionBase::SetGravity(std::vector<double>& gravity)
{
  gravity_ = gravity;
}

/*---------------------------------------------------------------------------*
 | maximum particle radius (on this processor)                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::ParticleInteractionBase::MaxParticleRadius() const
{
  // init value of maximum radius
  double maxrad = 0.0;

  // iterate over particle types
  for (auto& typeIt : particlecontainerbundle_->GetRefToAllContainersMap())
  {
    // get container of owned particles of current particle type
    auto statusIt = (typeIt.second).find(PARTICLEENGINE::Owned);
    if (statusIt == (typeIt.second).end())
      dserror("particle status '%s' not found!",
          PARTICLEENGINE::EnumToStatusName(PARTICLEENGINE::Owned).c_str());
    PARTICLEENGINE::ParticleContainerShrdPtr container = statusIt->second;

    // get number of particles stored in container
    int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get pointer to particle radius
    double* rad = container->GetPtrToParticleState(PARTICLEENGINE::Radius, 0);

    // iterate over owned particles of current type
    for (int i = 0; i < particlestored; ++i)
      if (rad[i] > maxrad) maxrad = rad[i];
  }

  return maxrad;
}
