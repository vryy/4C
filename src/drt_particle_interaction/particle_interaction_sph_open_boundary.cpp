/*---------------------------------------------------------------------------*/
/*! \file
\brief open boundary handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_open_boundary.H"

#include "particle_interaction_sph_kernel.H"
#include "particle_interaction_material_handler.H"
#include "particle_interaction_sph_equationofstate.H"
#include "particle_interaction_sph_equationofstate_bundle.H"
#include "particle_interaction_sph_neighbor_pairs.H"

#include "particle_interaction_utils.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"
#include "../drt_particle_engine/particle_object.H"

#include "../drt_lib/drt_dserror.H"

#include "../drt_lib/drt_globalproblem.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHOpenBoundaryBase::SPHOpenBoundaryBase(const Teuchos::ParameterList& params)
    : params_sph_(params),
      prescribedstatefunctid_(-1),
      fluidphase_(PARTICLEENGINE::Phase1),
      openboundaryphase_(PARTICLEENGINE::DirichletPhase)
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHOpenBoundaryBase::Init()
{
  // nothing to do
}

void PARTICLEINTERACTION::SPHOpenBoundaryBase::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
    const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
    const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle,
    const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // set kernel handler
  kernel_ = kernel;

  // set particle material handler
  particlematerial_ = particlematerial;

  // set equation of state handler
  equationofstatebundle_ = equationofstatebundle;

  // set neighbor pair handler
  neighborpairs_ = neighborpairs;

  // safety check
  for (const auto& type_i : {fluidphase_, openboundaryphase_})
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i))
      dserror("no particle container for particle type '%s' found!",
          PARTICLEENGINE::EnumToTypeName(type_i).c_str());
}

void PARTICLEINTERACTION::SPHOpenBoundaryBase::CheckOpenBoundaryPhaseChange(
    const double maxinteractiondistance)
{
  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  std::vector<std::set<int>> particlestoremove(typevectorsize);
  std::vector<std::vector<std::pair<int, PARTICLEENGINE::ParticleObjShrdPtr>>> particlestoinsert(
      typevectorsize);

  // global ids to be freed
  std::vector<int> globalidstofree;

  // generated particle objects
  std::vector<PARTICLEENGINE::ParticleObjShrdPtr> particlesgenerated;

  // get initial particle spacing
  const double initialparticlespacing = params_sph_.get<double>("INITIALPARTICLESPACING");

  // number of particles per direction
  const int numparticleperdir = std::round(maxinteractiondistance / initialparticlespacing);

  // tolerance for phase change of open boundary particle to fluid particle
  const double toleranceopenboundarytofluid = 0.1 * initialparticlespacing;

  // get container of owned particles of open boundary phase
  PARTICLEENGINE::ParticleContainer* container_i =
      particlecontainerbundle_->GetSpecificContainer(openboundaryphase_, PARTICLEENGINE::Owned);

  // iterate over particles in container
  for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
  {
    // get pointer to particle states
    double* pos_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i);

    // compute distance of open boundary particle from plane
    std::vector<double> temp(3);
    UTILS::vec_set(&temp[0], pos_i);
    UTILS::vec_sub(&temp[0], &planepoint_[0]);

    const double distancefromplane = UTILS::vec_dot(&temp[0], &outwardnormal_[0]);

    // open boundary particle traveled over plane
    if (distancefromplane < -toleranceopenboundarytofluid)
    {
      // duplicate open boundary particle and change phase
      int globalid(0);
      PARTICLEENGINE::ParticleStates particlestates;
      container_i->GetParticle(particle_i, globalid, particlestates);

      // construct and store generated particle object
      particlesgenerated.emplace_back(
          std::make_shared<PARTICLEENGINE::ParticleObject>(fluidphase_, -1, particlestates));

      // shift open boundary particle back
      UTILS::vec_addscale(pos_i, numparticleperdir * initialparticlespacing, &outwardnormal_[0]);
    }
    // open boundary particle more than maximum interaction distance away from plane
    else if (distancefromplane > maxinteractiondistance)
    {
      // get global id of particle i
      const int* globalid_i = container_i->GetPtrToParticleGlobalID(particle_i);

      // store global id of open boundary particle to be freed
      globalidstofree.push_back(globalid_i[0]);

      // store index of open boundary particle to be removed from container
      particlestoremove[openboundaryphase_].insert(particle_i);
    }
  }

  // get container of owned particles of fluid phase
  PARTICLEENGINE::ParticleContainer* container_j =
      particlecontainerbundle_->GetSpecificContainer(fluidphase_, PARTICLEENGINE::Owned);

  // iterate over particles in container
  for (int particle_j = 0; particle_j < container_j->ParticlesStored(); ++particle_j)
  {
    // get pointer to particle states
    double* pos_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Position, particle_j);

    // compute distance of fluid particle from plane
    std::vector<double> temp(3);
    UTILS::vec_set(&temp[0], pos_j);
    UTILS::vec_sub(&temp[0], &planepoint_[0]);

    const double distancefromplane = UTILS::vec_dot(&temp[0], &outwardnormal_[0]);

    // fluid particle traveled over plane
    if (distancefromplane > 0.0)
    {
      int globalid(0);
      PARTICLEENGINE::ParticleStates particlestates;
      container_j->GetParticle(particle_j, globalid, particlestates);

      PARTICLEENGINE::ParticleObjShrdPtr particleobject =
          std::make_shared<PARTICLEENGINE::ParticleObject>(
              openboundaryphase_, globalid, particlestates);

      // append particle to be insert
      particlestoinsert[openboundaryphase_].push_back(std::make_pair(-1, particleobject));

      // store index of fluid particle to be removed from container
      particlestoremove[fluidphase_].insert(particle_j);
    }
  }

  // get unique global ids for all particles
  particleengineinterface_->GetUniqueGlobalIdsForAllParticles(particlesgenerated);

  // append particle to be insert
  for (auto& particleobject : particlesgenerated)
    particlestoinsert[particleobject->ReturnParticleType()].push_back(
        std::make_pair(-1, particleobject));

  // free unique global ids
  particleengineinterface_->FreeUniqueGlobalIds(globalidstofree);

  // hand over particles to be removed
  particleengineinterface_->HandOverParticlesToBeRemoved(particlestoremove);

  // hand over particles to be inserted
  particleengineinterface_->HandOverParticlesToBeInserted(particlestoinsert);
}

PARTICLEINTERACTION::SPHOpenBoundaryDirichlet::SPHOpenBoundaryDirichlet(
    const Teuchos::ParameterList& params)
    : SPHOpenBoundaryBase::SPHOpenBoundaryBase(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHOpenBoundaryDirichlet::Init()
{
  // call base class init
  SPHOpenBoundaryBase::Init();

  // init function id of prescribed state
  prescribedstatefunctid_ = params_sph_.get<int>("DIRICHLET_FUNCT");

  // safety check
  if (not(prescribedstatefunctid_ > 0)) dserror("no function id of prescribed state set!");

  // init outward normal
  {
    double value;
    std::istringstream stream(
        Teuchos::getNumericStringParameter(params_sph_, "DIRICHLET_OUTWARD_NORMAL"));
    while (stream >> value) outwardnormal_.push_back(value);

    // safety check
    if (static_cast<int>(outwardnormal_.size()) != 3)
      dserror("dimension (dim = %d) of outward normal is wrong!",
          static_cast<int>(outwardnormal_.size()));

    // normalize outward normal
    const double norm = UTILS::vec_norm2(&outwardnormal_[0]);
    if (not(norm > 0.0)) dserror("no outward normal set!");
    UTILS::vec_setscale(&outwardnormal_[0], 1.0 / norm, &outwardnormal_[0]);
  }

  // init plain point
  {
    double value;
    std::istringstream stream(
        Teuchos::getNumericStringParameter(params_sph_, "DIRICHLET_PLANE_POINT"));
    while (stream >> value) planepoint_.push_back(value);

    // safety check
    if (static_cast<int>(planepoint_.size()) != 3)
      dserror(
          "dimension (dim = %d) of plane point is wrong!", static_cast<int>(planepoint_.size()));
  }

  // init fluid phase and open boundary phase
  fluidphase_ = PARTICLEENGINE::Phase1;
  openboundaryphase_ = PARTICLEENGINE::DirichletPhase;
}

void PARTICLEINTERACTION::SPHOpenBoundaryDirichlet::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
    const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
    const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle,
    const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs)
{
  // call base class setup
  SPHOpenBoundaryBase::Setup(
      particleengineinterface, kernel, particlematerial, equationofstatebundle, neighborpairs);

  // setup states of ghosted particles to refresh
  {
    std::vector<PARTICLEENGINE::StateEnum> states{
        PARTICLEENGINE::Density, PARTICLEENGINE::Pressure};

    statestorefresh_.push_back(std::make_pair(openboundaryphase_, states));
  }
}

void PARTICLEINTERACTION::SPHOpenBoundaryDirichlet::PrescribeOpenBoundaryStates(
    const double& evaltime)
{
  // get container of owned particles of open boundary phase
  PARTICLEENGINE::ParticleContainer* container_i =
      particlecontainerbundle_->GetSpecificContainer(openboundaryphase_, PARTICLEENGINE::Owned);

  // get number of particles stored in container
  const int particlestored = container_i->ParticlesStored();

  // no owned particles of current particle type
  if (particlestored <= 0) return;

  // get reference to function
  DRT::UTILS::Function& function = DRT::Problem::Instance()->Funct(prescribedstatefunctid_ - 1);

  // safety check
  if (function.NumberComponents() != 1)
    dserror("dimension of function governing velocity condition is not one!");

  // iterate over particles in container
  for (int particle_i = 0; particle_i < particlestored; ++particle_i)
  {
    // get pointer to particle states
    const double* pos_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i);
    double* vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);

    // evaluate function to set velocity
    UTILS::vec_setscale(vel_i, -function.Evaluate(0, pos_i, evaltime), &outwardnormal_[0]);
  }
}

void PARTICLEINTERACTION::SPHOpenBoundaryDirichlet::InterpolateOpenBoundaryStates()
{
  // get container of owned particles of open boundary phase
  PARTICLEENGINE::ParticleContainer* container_k =
      particlecontainerbundle_->GetSpecificContainer(openboundaryphase_, PARTICLEENGINE::Owned);

  // get material for current particle type
  const MAT::PAR::ParticleMaterialBase* material_k =
      particlematerial_->GetPtrToParticleMatParameter(openboundaryphase_);

  // get equation of state for current particle type
  const PARTICLEINTERACTION::SPHEquationOfStateBase* equationofstate_k =
      equationofstatebundle_->GetPtrToSpecificEquationOfState(openboundaryphase_);

  // get number of particles stored in container
  const int particlestored = container_k->ParticlesStored();

  std::vector<double> sumj_Vj_Wij(particlestored, 0.0);
  std::vector<double> sumj_Vj_Wij_pj(particlestored, 0.0);

  // get relevant particle pair indices
  std::vector<int> relindices;
  neighborpairs_->GetRelevantParticlePairIndicesForDisjointCombination(
      {openboundaryphase_}, {fluidphase_}, relindices);

  // iterate over relevant particle pairs
  for (const int particlepairindex : relindices)
  {
    const SPHParticlePair& particlepair =
        neighborpairs_->GetRefToParticlePairData()[particlepairindex];

    // access values of local index tuples of particle i and j
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = particlepair.tuple_i_;

    PARTICLEENGINE::TypeEnum type_j;
    PARTICLEENGINE::StatusEnum status_j;
    int particle_j;
    std::tie(type_j, status_j, particle_j) = particlepair.tuple_j_;

    // get corresponding particle containers
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    PARTICLEENGINE::ParticleContainer* container_j =
        particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

    // get pointer to particle states
    const double* mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
    const double* dens_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);
    const double* press_i =
        container_i->GetPtrToParticleState(PARTICLEENGINE::Pressure, particle_i);

    // get pointer to particle states
    const double* mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);
    const double* dens_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j);
    const double* press_j =
        container_j->GetPtrToParticleState(PARTICLEENGINE::Pressure, particle_j);

    // evaluate contribution of neighboring particle j
    if (type_i == openboundaryphase_)
    {
      const double fac = mass_j[0] / dens_j[0] * particlepair.Wij_;
      sumj_Vj_Wij[particle_i] += fac;
      sumj_Vj_Wij_pj[particle_i] += fac * press_j[0];
    }

    // evaluate contribution of neighboring particle i
    if (type_j == openboundaryphase_ and status_j == PARTICLEENGINE::Owned)
    {
      const double fac = mass_i[0] / dens_i[0] * particlepair.Wji_;
      sumj_Vj_Wij[particle_j] += fac;
      sumj_Vj_Wij_pj[particle_j] += fac * press_i[0];
    }
  }

  // iterate over particles in container
  for (int particle_k = 0; particle_k < particlestored; ++particle_k)
  {
    // get pointer to particle states
    double* dens_k = container_k->GetPtrToParticleState(PARTICLEENGINE::Density, particle_k);
    double* press_k = container_k->GetPtrToParticleState(PARTICLEENGINE::Pressure, particle_k);

    // interpolate pressure
    press_k[0] = (sumj_Vj_Wij[particle_k] > 0.0)
                     ? sumj_Vj_Wij_pj[particle_k] / sumj_Vj_Wij[particle_k]
                     : 0.0;

    // compute density
    dens_k[0] = equationofstate_k->PressureToDensity(press_k[0], material_k->initDensity_);
  }

  // refresh states of ghosted particles
  particleengineinterface_->RefreshParticlesOfSpecificStatesAndTypes(statestorefresh_);
}

PARTICLEINTERACTION::SPHOpenBoundaryNeumann::SPHOpenBoundaryNeumann(
    const Teuchos::ParameterList& params)
    : SPHOpenBoundaryBase::SPHOpenBoundaryBase(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHOpenBoundaryNeumann::Init()
{
  // call base class init
  SPHOpenBoundaryBase::Init();

  // init function id of prescribed state
  prescribedstatefunctid_ = params_sph_.get<int>("NEUMANN_FUNCT");

  // init outward normal
  {
    double value;
    std::istringstream stream(
        Teuchos::getNumericStringParameter(params_sph_, "NEUMANN_OUTWARD_NORMAL"));
    while (stream >> value) outwardnormal_.push_back(value);

    // safety check
    if (static_cast<int>(outwardnormal_.size()) != 3)
      dserror("dimension (dim = %d) of outward normal is wrong!",
          static_cast<int>(outwardnormal_.size()));

    // normalize outward normal
    const double direction_norm = UTILS::vec_norm2(&outwardnormal_[0]);
    if (not(direction_norm > 0.0)) dserror("no outward normal set!");
    UTILS::vec_setscale(&outwardnormal_[0], 1.0 / direction_norm, &outwardnormal_[0]);
  }

  // init plain point
  {
    double value;
    std::istringstream stream(
        Teuchos::getNumericStringParameter(params_sph_, "NEUMANN_PLANE_POINT"));
    while (stream >> value) planepoint_.push_back(value);

    // safety check
    if (static_cast<int>(planepoint_.size()) != 3)
      dserror(
          "dimension (dim = %d) of plane point is wrong!", static_cast<int>(planepoint_.size()));
  }

  // init fluid phase and open boundary phase
  fluidphase_ = PARTICLEENGINE::Phase1;
  openboundaryphase_ = PARTICLEENGINE::NeumannPhase;
}

void PARTICLEINTERACTION::SPHOpenBoundaryNeumann::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
    const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
    const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle,
    const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs)
{
  // call base class setup
  SPHOpenBoundaryBase::Setup(
      particleengineinterface, kernel, particlematerial, equationofstatebundle, neighborpairs);

  // setup states of ghosted particles to refresh
  {
    std::vector<PARTICLEENGINE::StateEnum> states{PARTICLEENGINE::Velocity};

    statestorefresh_.push_back(std::make_pair(openboundaryphase_, states));
  }
}

void PARTICLEINTERACTION::SPHOpenBoundaryNeumann::PrescribeOpenBoundaryStates(
    const double& evaltime)
{
  // get container of owned particles of open boundary phase
  PARTICLEENGINE::ParticleContainer* container_i =
      particlecontainerbundle_->GetSpecificContainer(openboundaryphase_, PARTICLEENGINE::Owned);

  // get material for current particle type
  const MAT::PAR::ParticleMaterialBase* material_i =
      particlematerial_->GetPtrToParticleMatParameter(openboundaryphase_);

  // get equation of state for current particle type
  const PARTICLEINTERACTION::SPHEquationOfStateBase* equationofstate_i =
      equationofstatebundle_->GetPtrToSpecificEquationOfState(openboundaryphase_);

  // get number of particles stored in container
  const int particlestored = container_i->ParticlesStored();

  // no owned particles of current particle type
  if (particlestored <= 0) return;

  if (prescribedstatefunctid_ > 0)
  {
    // get reference to function
    DRT::UTILS::Function& function = DRT::Problem::Instance()->Funct(prescribedstatefunctid_ - 1);

    // safety check
    if (function.NumberComponents() != 1)
      dserror("dimension of function governing pressure condition is not one!");

    // iterate over particles in container
    for (int particle_i = 0; particle_i < particlestored; ++particle_i)
    {
      // get pointer to particle states
      const double* pos_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i);
      double* press_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Pressure, particle_i);

      // evaluate function to set pressure
      press_i[0] = function.Evaluate(0, pos_i, evaltime);
    }
  }
  else
  {
    // clear pressure state
    container_i->ClearState(PARTICLEENGINE::Pressure);
  }

  // iterate over particles in container
  for (int particle_i = 0; particle_i < particlestored; ++particle_i)
  {
    // get pointer to particle states
    const double* press_i =
        container_i->GetPtrToParticleState(PARTICLEENGINE::Pressure, particle_i);
    double* dens_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);

    // compute density
    dens_i[0] = equationofstate_i->PressureToDensity(press_i[0], material_i->initDensity_);
  }
}

void PARTICLEINTERACTION::SPHOpenBoundaryNeumann::InterpolateOpenBoundaryStates()
{
  // nothing to do
}
