/*---------------------------------------------------------------------------*/
/*!

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

#include "particle_interaction_material_handler.H"

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
  // init particle material handler
  InitParticleMaterialHandler();
}

/*---------------------------------------------------------------------------*
 | setup particle interaction handler                         sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionBase::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEALGORITHM::WallHandlerInterface> particlewallinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // set interface to particle wall hander
  particlewallinterface_ = particlewallinterface;

  // setup particle material handler
  particlematerial_->Setup();

  // init vector
  gravity_.resize(3, 0.0);
}

/*---------------------------------------------------------------------------*
 | write restart of particle interaction handler              sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionBase::WriteRestart(
    const int step, const double time) const
{
  // write restart of particle material handler
  particlematerial_->WriteRestart(step, time);
}

/*---------------------------------------------------------------------------*
 | read restart of particle interaction handler               sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // read restart of particle material handler
  particlematerial_->ReadRestart(reader);
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
 | init particle material handler                             sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::ParticleInteractionBase::InitParticleMaterialHandler()
{
  // create particle material handler
  particlematerial_ = std::make_shared<PARTICLEINTERACTION::MaterialHandler>(params_);

  // init particle material handler
  particlematerial_->Init();
}

/*---------------------------------------------------------------------------*
 | maximum particle radius (on this processor)                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::ParticleInteractionBase::MaxParticleRadius() const
{
  // init value of maximum radius
  double maxrad = 0.0;

  // iterate over particle types
  for (const auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    // get maximum stored value of state
    double currmaxrad = container->GetMaxValueOfState(PARTICLEENGINE::Radius);

    // compare to current maximum
    maxrad = std::max(maxrad, currmaxrad);
  }

  return maxrad;
}
