/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_dem_contact.cpp

\brief contact handler for discrete element method (DEM) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_dem_contact.H"

#include "particle_interaction_material_handler.H"
#include "particle_interaction_utils.H"

#include "particle_interaction_dem_neighbor_pairs.H"
#include "particle_interaction_dem_contact_normal.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMContact::DEMContact(const Teuchos::ParameterList& params)
    : params_dem_(params), dt_(0.0)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | destructor                                                 sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMContact::~DEMContact()
{
  // note: destructor declaration here since at compile-time a complete type
  // of class T as used in class member std::unique_ptr<T> ptr_T_ is required
}

/*---------------------------------------------------------------------------*
 | init contact handler                                       sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::Init()
{
  // init normal contact handler
  InitNormalContactHandler();
}

/*---------------------------------------------------------------------------*
 | setup contact handler                                      sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
    const std::shared_ptr<PARTICLEINTERACTION::DEMNeighborPairs> neighborpairs)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // set particle material handler
  particlematerial_ = particlematerial;

  // set neighbor pair handler
  neighborpairs_ = neighborpairs;

  //! get maximum density of all materials
  const double maxdensity = GetMaxDensityOfAllMaterials();

  // setup normal contact handler
  contactnormal_->Setup(maxdensity);
}

/*---------------------------------------------------------------------------*
 | write restart of contact handler                           sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::WriteRestart(const int step, const double time) const
{
  // write restart of normal contact handler
  contactnormal_->WriteRestart(step, time);
}

/*---------------------------------------------------------------------------*
 | read restart of contact handler                            sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // read restart of normal contact handler
  contactnormal_->ReadRestart(reader);
}

/*---------------------------------------------------------------------------*
 | set current step size                                      sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::SetCurrentStepSize(const double currentstepsize)
{
  dt_ = currentstepsize;
}

/*---------------------------------------------------------------------------*
 | check critical time step (on this processor)               sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::CheckCriticalTimeStep() const
{
  // init value of minimum mass
  double minmass = 0.0;

  // iterate over particle types
  for (const auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    // get minimum stored value of state
    double currminmass = container->GetMinValueOfState(PARTICLEENGINE::Mass);

    if ((not(minmass > 0.0)) or (currminmass < minmass)) minmass = currminmass;
  }

  // currently no particles in simulation domain
  if (minmass == 0.0) return;

  // get time critical contact stiffness
  const double k_tcrit = contactnormal_->GetTimeCriticalStiffness();

  // critical time step size based on particle-particle contact
  const double safety = 0.75;
  const double factor = 0.34;  // todo set to 0.22 if also tangential contact is considered!
  const double dt_crit = safety * factor * std::sqrt(minmass / k_tcrit);

  // checks time step
  if (dt_ > dt_crit) dserror("time step %f larger than critical time step %f!", dt_, dt_crit);
}

/*---------------------------------------------------------------------------*
 | add contact contribution to force field                    sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::AddForceContribution() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::DEMContact::AddForceContribution");

  // iterate over neighbor pairs
  for (const auto& neighborpair : neighborpairs_->GetRefToNeighborPairData())
  {
    // access values of local index tuples of particle i and j
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = neighborpair.tuple_i_;

    PARTICLEENGINE::TypeEnum type_j;
    PARTICLEENGINE::StatusEnum status_j;
    int particle_j;
    std::tie(type_j, status_j, particle_j) = neighborpair.tuple_j_;

    // get corresponding particle containers
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    PARTICLEENGINE::ParticleContainer* container_j =
        particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

    // declare pointer variables for particle i and j
    const double *vel_i, *rad_i;
    double* force_i;

    const double *vel_j, *rad_j;
    double* force_j;

    // get pointer to particle states
    vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);
    rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
    force_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Force, particle_i);

    vel_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_j);
    rad_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_j);
    force_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Force, particle_j);

    // relative velocity in contact point c between particle i and j
    double vel_rel[3];
    UTILS::vec_set(vel_rel, vel_i);
    UTILS::vec_sub(vel_rel, vel_j);

    // magnitude of relative velocity in normal direction
    const double vel_rel_normal = UTILS::vec_dot(vel_rel, neighborpair.e_ji_);

    // calculate normal contact force
    double normalcontactforce(0.0);
    contactnormal_->NormalContactForce(
        neighborpair.gap_, rad_i, rad_j, vel_rel_normal, neighborpair.m_eff_, normalcontactforce);

    // add normal contact force contribution
    UTILS::vec_addscale(force_i, normalcontactforce, neighborpair.e_ji_);
    if (status_j == PARTICLEENGINE::Owned)
      UTILS::vec_addscale(force_j, -normalcontactforce, neighborpair.e_ji_);
  }
}

/*---------------------------------------------------------------------------*
 | init normal contact handler                                sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::InitNormalContactHandler()
{
  // get type of discrete element method contact normal
  INPAR::PARTICLE::NormalContact contactnormaltype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::NormalContact>(params_dem_, "NORMALCONTACTLAW");

  // create normal contact handler
  switch (contactnormaltype)
  {
    case INPAR::PARTICLE::NormalLinSpring:
    {
      contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLinearSpring>(
          new PARTICLEINTERACTION::DEMContactNormalLinearSpring(params_dem_));
      break;
    }
    case INPAR::PARTICLE::NormalLinSpringDamp:
    {
      contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp>(
          new PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp(params_dem_));
      break;
    }
    case INPAR::PARTICLE::NormalHertz:
    {
      contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalHertz>(
          new PARTICLEINTERACTION::DEMContactNormalHertz(params_dem_));
      break;
    }
    case INPAR::PARTICLE::NormalLeeHerrmann:
    {
      contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLeeHerrmann>(
          new PARTICLEINTERACTION::DEMContactNormalLeeHerrmann(params_dem_));
      break;
    }
    case INPAR::PARTICLE::NormalKuwabaraKono:
    {
      contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalKuwabaraKono>(
          new PARTICLEINTERACTION::DEMContactNormalKuwabaraKono(params_dem_));
      break;
    }
    case INPAR::PARTICLE::NormalTsuji:
    {
      contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalTsuji>(
          new PARTICLEINTERACTION::DEMContactNormalTsuji(params_dem_));
      break;
    }
    default:
    {
      dserror("unknown normal contact law type!");
      break;
    }
  }

  // init normal contact handler
  contactnormal_->Init();
}

/*---------------------------------------------------------------------------*
 | get maximum density of all materials                       sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::DEMContact::GetMaxDensityOfAllMaterials() const
{
  // init value of maximum density
  double maxdensity(0.0);

  // iterate over particle types
  for (const auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // get material for current particle type
    const MAT::PAR::ParticleMaterialBase* material =
        particlematerial_->GetPtrToParticleMatParameter(typeEnum);

    // compare to current maximum
    maxdensity = std::max(maxdensity, material->initDensity_);
  }

  return maxdensity;
}
