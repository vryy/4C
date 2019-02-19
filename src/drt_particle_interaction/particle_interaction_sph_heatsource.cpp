/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_heatsource.cpp

\brief heat source handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_heatsource.H"

#include "particle_interaction_material_handler.H"
#include "particle_interaction_sph_neighbor_pairs.H"

#include "particle_interaction_utils.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_dserror.H"

#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHHeatSourceBase::SPHHeatSourceBase(const Teuchos::ParameterList& params)
    : params_sph_(params), heatsourcefctnumber_(params.get<int>("HEATSOURCE_FUNCT"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init heat source handler                                   sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHHeatSourceBase::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup heat source handler                                  sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHHeatSourceBase::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
    const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // set particle material handler
  particlematerial_ = particlematerial;

  // set neighbor pair handler
  neighborpairs_ = neighborpairs;

  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  // allocate memory to hold particle types
  thermomaterial_.resize(typevectorsize);

  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    thermomaterial_[type_i] = dynamic_cast<const MAT::PAR::ParticleMaterialThermo*>(
        particlematerial_->GetPtrToParticleMatParameter(type_i));

    // determine absorbing particle types
    if (thermomaterial_[type_i]->thermalAbsorptivity_ > 0.0) absorbingtypes_.insert(type_i);
  }

  // safety check
  if (absorbingtypes_.count(PARTICLEENGINE::BoundaryPhase))
    dserror("no heat source evaluation for boundary particles!");
}

/*---------------------------------------------------------------------------*
 | write restart of heat source handler                       sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHHeatSourceBase::WriteRestart(const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of heat source handler                        sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHHeatSourceBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHHeatSourceVolume::SPHHeatSourceVolume(const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::SPHHeatSourceBase(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | evaluate heat source                                       sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHHeatSourceVolume::EvaluateHeatSource(const double& evaltime) const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHHeatSourceVolume::EvaluateHeatSource");

  // init vector containing evaluated function
  std::vector<double> funct(1);

  // get reference to function
  DRT::UTILS::Function& function = DRT::Problem::Instance()->Funct(heatsourcefctnumber_ - 1);

  // safety check
  if (function.NumberComponents() != 1)
    dserror("dimension of function defining heat source is not one!");

  // iterate over absorbing particle types
  for (const auto& type_i : absorbingtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get material for current particle type
    const MAT::PAR::ParticleMaterialThermo* thermomaterial_i = thermomaterial_[type_i];

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
    {
      // declare pointer variables for particle i
      const double* pos_i;
      double* tempdot_i;

      // get pointer to particle states
      pos_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i);
      tempdot_i = container_i->GetPtrToParticleState(PARTICLEENGINE::TemperatureDot, particle_i);

      // evaluate function defining heat source
      funct = function.EvaluateTimeDerivative(0, &(pos_i[0]), evaltime, 0);

      // add contribution of heat source
      tempdot_i[0] +=
          thermomaterial_i->thermalAbsorptivity_ * funct[0] * thermomaterial_i->invThermalCapacity_;
    }
  }
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHHeatSourceSurface::SPHHeatSourceSurface(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::SPHHeatSourceBase(params), eval_direction_(false)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init heat source handler                                   sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHHeatSourceSurface::Init()
{
  // call base class init
  SPHHeatSourceBase::Init();

  // init heat source direction vector
  double value;
  std::istringstream directionstream(
      Teuchos::getNumericStringParameter(params_sph_, "HEATSOURCE_DIRECTION"));

  while (directionstream >> value) direction_.push_back(value);

  // safety check
  if ((int)direction_.size() != 3)
    dserror(
        "dimension (dim = %d) of heat source direction vector is wrong!", (int)direction_.size());

  // normalize heat source direction vector
  const double direction_norm = UTILS::vec_norm2(&direction_[0]);
  if (direction_norm > 0.0)
  {
    eval_direction_ = true;
    UTILS::vec_setscale(&direction_[0], 1.0 / direction_norm, &direction_[0]);
  }
}

/*---------------------------------------------------------------------------*
 | evaluate heat source                                       sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHHeatSourceSurface::EvaluateHeatSource(const double& evaltime) const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHHeatSourceSurface::EvaluateHeatSource");

  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--absorbingtypes_.end()) + 1;

  // colorfield gradient of absorbing interface particles
  std::vector<std::vector<std::vector<double>>> cfg_i(typevectorsize);

  // iterate over absorbing particle types
  for (const auto& type_i : absorbingtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container_i->ParticlesStored();

    // allocate memory
    cfg_i[type_i].assign(particlestored, std::vector<double>(3, 0.0));
  }

  // get reference to index of neighbor pairs for each type
  const SPHIndexOfNeighborPairs& indexofneighborpairs =
      neighborpairs_->GetRefToIndexOfNeighborPairs();

  // determine relevant neighbor pair indices
  std::vector<int> relindices;

  // iterate over absorbing particle types
  for (const auto& type_i : absorbingtypes_)
    relindices.insert(
        relindices.end(), indexofneighborpairs[type_i].begin(), indexofneighborpairs[type_i].end());

  // sort and erase duplicate indices of relevant neighbor pairs
  if (absorbingtypes_.size() > 1)
  {
    std::sort(relindices.begin(), relindices.end());
    relindices.erase(std::unique(relindices.begin(), relindices.end()), relindices.end());
  }

  // iterate over relevant neighbor pairs
  for (const int neighborpairindex : relindices)
  {
    const SPHNeighborPair& neighborpair =
        neighborpairs_->GetRefToNeighborPairData()[neighborpairindex];

    // access values of local index tuples of particle i and j
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = neighborpair.tuple_i_;

    PARTICLEENGINE::TypeEnum type_j;
    PARTICLEENGINE::StatusEnum status_j;
    int particle_j;
    std::tie(type_j, status_j, particle_j) = neighborpair.tuple_j_;

    // no evaluation for boundary particles
    if (type_i == PARTICLEENGINE::BoundaryPhase or type_j == PARTICLEENGINE::BoundaryPhase)
      continue;

    // check for absorbing particles
    bool isabsorbing_i = absorbingtypes_.count(type_i);
    bool isabsorbing_j = absorbingtypes_.count(type_j);

    // no evaluation for both absorbing particles
    if (isabsorbing_i and isabsorbing_j) continue;

    // get corresponding particle containers
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    PARTICLEENGINE::ParticleContainer* container_j =
        particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

    // get material for particle types
    const MAT::PAR::ParticleMaterialBase* material_i =
        particlematerial_->GetPtrToParticleMatParameter(type_i);

    const MAT::PAR::ParticleMaterialBase* material_j =
        particlematerial_->GetPtrToParticleMatParameter(type_j);

    // declare pointer variables for particle i and j
    const double *mass_i, *dens_i;
    const double *mass_j, *dens_j;

    // get pointer to particle states
    mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);

    dens_i = (type_i == PARTICLEENGINE::RigidPhase)
                 ? &(material_i->initDensity_)
                 : container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);

    mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);

    dens_j = (type_j == PARTICLEENGINE::RigidPhase)
                 ? &(material_j->initDensity_)
                 : container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j);

    const double fac =
        (UTILS::pow<2>(mass_i[0] / dens_i[0]) + UTILS::pow<2>(mass_j[0] / dens_j[0])) /
        (dens_i[0] + dens_j[0]);

    // evaluate contribution of neighboring particle j
    if (isabsorbing_i)
    {
      // sum contribution of neighboring particle j
      UTILS::vec_addscale(&cfg_i[type_i][particle_i][0], fac * dens_i[0] * neighborpair.dWdrij_,
          neighborpair.e_ij_);
    }

    // evaluate contribution of neighboring particle i
    if (isabsorbing_j and status_j == PARTICLEENGINE::Owned)
    {
      // sum contribution of neighboring particle i
      UTILS::vec_addscale(&cfg_i[type_j][particle_j][0], -fac * dens_j[0] * neighborpair.dWdrji_,
          neighborpair.e_ij_);
    }
  }

  // init vector containing evaluated function
  std::vector<double> funct(1);

  // get reference to function
  DRT::UTILS::Function& function = DRT::Problem::Instance()->Funct(heatsourcefctnumber_ - 1);

  // safety check
  if (function.NumberComponents() != 1)
    dserror("dimension of function defining heat source is not one!");

  // iterate over absorbing particle types
  for (const auto& type_i : absorbingtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get material for current particle type
    const MAT::PAR::ParticleMaterialThermo* thermomaterial_i = thermomaterial_[type_i];

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
    {
      // norm of colorfield gradient of absorbing interface particles
      const double f_i = UTILS::vec_norm2(&cfg_i[type_i][particle_i][0]);

      // no heat source contribution to current particle
      if (not(f_i > 0.0)) continue;

      // projection of colorfield gradient with heat source direction
      const double f_i_proj =
          eval_direction_ ? UTILS::vec_dot(&direction_[0], &cfg_i[type_i][particle_i][0]) : f_i;

      // declare pointer variables for particle i
      const double* pos_i;
      double* tempdot_i;

      // get pointer to particle states
      pos_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i);
      tempdot_i = container_i->GetPtrToParticleState(PARTICLEENGINE::TemperatureDot, particle_i);

      // evaluate function defining heat source
      funct = function.EvaluateTimeDerivative(0, &(pos_i[0]), evaltime, 0);

      // add contribution of heat source
      tempdot_i[0] += f_i_proj * thermomaterial_i->thermalAbsorptivity_ * funct[0] *
                      thermomaterial_i->invThermalCapacity_;
    }
  }
}
