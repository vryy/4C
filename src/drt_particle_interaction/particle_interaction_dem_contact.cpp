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
#include "particle_interaction_dem_neighbor_pairs.H"
#include "particle_interaction_dem_contact_normal.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_dserror.H"

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
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    // get minimum stored value of state
    double currminmass = container->GetMinValueOfState(PARTICLEENGINE::Mass);

    if ((not(minmass > 0.0)) or (currminmass < minmass)) minmass = currminmass;
  }

  // currently no particles in simulation domain
  if (minmass == 0) return;

  // get time critical contact stiffness
  const double k_tcrit = contactnormal_->GetTimeCriticalStiffness();

  // critical time step size based on particle-particle contact
  const double safety = 0.75;
  const double factor = 0.34;  // todo set to 0.22 if also tangential contact is considered!
  double dt_crit = safety * factor * std::sqrt(minmass / k_tcrit);

  // checks time step
  if (dt_ > dt_crit) dserror("time step %f larger than critical time step %f!", dt_, dt_crit);
}

/*---------------------------------------------------------------------------*
 | add contact contribution to force field                    sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::AddForceContribution() const
{
  // get reference to neighbor pair data
  const DEMNeighborPairData& neighborpairdata = neighborpairs_->GetRefToNeighborPairData();

  // iterate over particle types
  for (auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
    {
      // get reference to vector of neighbor pairs of current particle
      const std::vector<DEMNeighborPair>& currentNeighborPairs =
          (neighborpairdata[type_i])[particle_i];

      // check for neighbor pairs of current particle
      if (currentNeighborPairs.empty()) continue;

      // declare pointer variables for particle i
      const double *vel_i, *rad_i;
      double* force_i;

      // get pointer to particle states
      vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);
      rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
      force_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Force, particle_i);

      // iterate over neighbor pairs
      for (auto& neighborIt : currentNeighborPairs)
      {
        // access values of local index tuple of neighboring particle
        PARTICLEENGINE::TypeEnum type_j;
        PARTICLEENGINE::StatusEnum status_j;
        int particle_j;
        std::tie(type_j, status_j, particle_j) = neighborIt.first;

        // get reference to current particle pair
        const DEMParticlePair& particlepair = neighborIt.second;

        // get container of neighboring particles of current particle type and state
        PARTICLEENGINE::ParticleContainerShrdPtr container_j =
            particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

        // declare pointer variables for neighbor particle j
        const double *vel_j, *rad_j;

        // get pointer to particle states
        vel_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_j);
        rad_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_j);

        // evaluate relative velocity in normal direction
        double vel_rel_normal(0.0), vel_rel[3];

        // iterate over spatial dimension
        for (int dim = 0; dim < 3; ++dim)
        {
          // relative velocity between particles
          vel_rel[dim] = vel_i[dim] - vel_j[dim];

          // relative velocity between particles in normal direction
          vel_rel_normal += vel_rel[dim] * particlepair.e_ji_[dim];
        }

        // calculate normal contact force
        double normalcontactforce(0.0);
        contactnormal_->NormalContactForce(particlepair.gap_, rad_i, rad_j, vel_rel_normal,
            particlepair.m_eff_, normalcontactforce);

        // add contact contribution
        for (int i = 0; i < 3; ++i) force_i[i] += normalcontactforce * particlepair.e_ji_[i];
      }
    }
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
    case INPAR::PARTICLE::LinSpring:
    {
      contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLinearSpring>(
          new PARTICLEINTERACTION::DEMContactNormalLinearSpring(params_dem_));
      break;
    }
    case INPAR::PARTICLE::LinSpringDamp:
    {
      contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp>(
          new PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp(params_dem_));
      break;
    }
    case INPAR::PARTICLE::Hertz:
    {
      contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalHertz>(
          new PARTICLEINTERACTION::DEMContactNormalHertz(params_dem_));
      break;
    }
    case INPAR::PARTICLE::LeeHerrmann:
    {
      contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLeeHerrmann>(
          new PARTICLEINTERACTION::DEMContactNormalLeeHerrmann(params_dem_));
      break;
    }
    case INPAR::PARTICLE::KuwabaraKono:
    {
      contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalKuwabaraKono>(
          new PARTICLEINTERACTION::DEMContactNormalKuwabaraKono(params_dem_));
      break;
    }
    case INPAR::PARTICLE::Tsuji:
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
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // get material for current particle type
    const MAT::PAR::ParticleMaterialBase* material =
        particlematerial_->GetPtrToParticleMatParameter(typeEnum);

    if (material->initDensity_ > maxdensity) maxdensity = material->initDensity_;
  }

  return maxdensity;
}
