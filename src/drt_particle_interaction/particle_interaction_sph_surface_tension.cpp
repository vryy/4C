/*---------------------------------------------------------------------------*/
/*! \file
\brief surface tension handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_surface_tension.H"

#include "particle_interaction_sph_kernel.H"
#include "particle_interaction_material_handler.H"
#include "particle_interaction_sph_equationofstate.H"
#include "particle_interaction_sph_equationofstate_bundle.H"
#include "particle_interaction_sph_neighbor_pairs.H"

#include "particle_interaction_utils.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_dserror.H"

#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHSurfaceTensionBase::SPHSurfaceTensionBase(
    const Teuchos::ParameterList& params)
    : params_sph_(params),
      time_(0.0),
      surfacetensionrampfctnumber_(params.get<int>("SURFACETENSION_RAMP_FUNCT"))
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHSurfaceTensionBase::Init()
{
  // init with potential fluid particle types
  fluidtypes_ = {PARTICLEENGINE::Phase1, PARTICLEENGINE::Phase2};

  // init with potential boundary particle types
  boundarytypes_ = {PARTICLEENGINE::BoundaryPhase, PARTICLEENGINE::RigidPhase};
}

void PARTICLEINTERACTION::SPHSurfaceTensionBase::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
    const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
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

  // set neighbor pair handler
  neighborpairs_ = neighborpairs;

  // update with actual fluid particle types
  for (const auto& type_i : fluidtypes_)
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i)) fluidtypes_.erase(type_i);

  // update with actual boundary particle types
  for (const auto& type_i : boundarytypes_)
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i))
      boundarytypes_.erase(type_i);
}

void PARTICLEINTERACTION::SPHSurfaceTensionBase::WriteRestart(
    const int step, const double time) const
{
  // nothing to do
}

void PARTICLEINTERACTION::SPHSurfaceTensionBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

void PARTICLEINTERACTION::SPHSurfaceTensionBase::SetCurrentTime(const double currenttime)
{
  time_ = currenttime;
}

PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::SPHSurfaceTensionContinuumSurfaceForce(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::SPHSurfaceTensionBase(params),
      alpha0_(params_sph_.get<double>("SURFACETENSIONCOEFFICIENT")),
      staticcontactangle_(params_sph_.get<double>("STATICCONTACTANGLE")),
      alphaT_(params_sph_.get<double>("SURFACETENSIONTEMPFAC")),
      reftemp_(params_sph_.get<double>("SURFACETENSIONREFTEMP"))
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::Init()
{
  // call base class init
  SPHSurfaceTensionBase::Init();

  // safety check
  if (not(alpha0_ > 0.0)) dserror("constant factor of surface tension coefficient not positive!");

  if (alphaT_ != 0.0)
  {
    if (DRT::INPUT::IntegralValue<INPAR::PARTICLE::TemperatureEvaluationScheme>(
            params_sph_, "TEMPERATUREEVALUATION") == INPAR::PARTICLE::NoTemperatureEvaluation)
      dserror("temperature evaluation needed for temperature dependent surface tension!");

    if (DRT::INPUT::IntegralValue<int>(params_sph_, "TEMPERATUREGRADIENT") == false)
      dserror("temperature gradient evaluation needed for temperature dependent surface tension!");
  }
}

void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
    const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
    const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs)
{
  // call base class setup
  SPHSurfaceTensionBase::Setup(particleengineinterface, kernel, particlematerial, neighborpairs);

  // setup colorfield gradient and wall distance of ghosted particles to refresh
  {
    std::vector<PARTICLEENGINE::StateEnum> states{
        PARTICLEENGINE::ColorfieldGradient, PARTICLEENGINE::WallDistance};

    // iterate over fluid particle types
    for (const auto& type_i : fluidtypes_)
      cfgwalldisttorefresh_.push_back(std::make_pair(type_i, states));
  }

  // setup colorfield gradient and interface normal of ghosted particles to refresh
  {
    std::vector<PARTICLEENGINE::StateEnum> states{
        PARTICLEENGINE::ColorfieldGradient, PARTICLEENGINE::InterfaceNormal};

    // iterate over fluid particle types
    for (const auto& type_i : fluidtypes_)
      cfgintnormtorefresh_.push_back(std::make_pair(type_i, states));
  }

  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  // allocate memory to hold particle types
  f_i_.resize(typevectorsize, std::vector<std::vector<double>>(2));
}

void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::
    InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes)
{
  bool haveboundarytypes = false;

  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type = typeIt.first;

    if (boundarytypes_.count(type)) haveboundarytypes = true;
  }

  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type = typeIt.first;

    // set of particle states for current particle type
    std::set<PARTICLEENGINE::StateEnum>& particlestates = typeIt.second;

    // current particle type is not a fluid particle type
    if (not fluidtypes_.count(type)) continue;

    // states for surface tension evaluation scheme
    particlestates.insert({PARTICLEENGINE::ColorfieldGradient, PARTICLEENGINE::InterfaceNormal});

    if (haveboundarytypes)
      particlestates.insert({PARTICLEENGINE::UnitWallNormal, PARTICLEENGINE::WallDistance});
  }
}

void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::AddAccelerationContribution()
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::AddAccelerationContribution");

  // compute colorfield gradient
  ComputeColorfieldGradient();

  if (not boundarytypes_.empty())
  {
    // compute unit wall normal
    ComputeUnitWallNormal();

    // compute wall distance
    ComputeWallDistance();

    // refresh colorfield gradient and wall distance
    particleengineinterface_->RefreshParticlesOfSpecificStatesAndTypes(cfgwalldisttorefresh_);

    // compute wall correction factor
    ComputeWallCorrectionFactor();

    // extrapolate colorfield gradient at triple point
    ExtrapolateColorfieldGradientAtTriplePoint();
  }

  // compute interface normal
  ComputeInterfaceNormal();

  if (not boundarytypes_.empty())
  {
    // correct normal vector of particles close to triple point
    CorrectTriplePointNormal();
  }

  // refresh colorfield gradient and interface normal
  particleengineinterface_->RefreshParticlesOfSpecificStatesAndTypes(cfgintnormtorefresh_);

  // compute surface tension contribution
  ComputeSurfaceTensionContribution();

  // compute temperature gradient driven contribution
  if (alphaT_ != 0.0) ComputeTempGradDrivenContribution();
}

void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::ComputeColorfieldGradient() const
{
  // iterate over fluid particle types
  for (const auto& type_i : fluidtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear colorfield gradient state
    container_i->ClearState(PARTICLEENGINE::ColorfieldGradient);
  }

  // get relevant particle pair indices
  std::vector<int> relindices;
  neighborpairs_->GetRelevantParticlePairIndicesForEqualCombination(fluidtypes_, relindices);

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

    // no evaluation for particles of same type
    if (type_i == type_j) continue;

    // get corresponding particle containers
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    PARTICLEENGINE::ParticleContainer* container_j =
        particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

    // get pointer to particle states
    const double* mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
    const double* dens_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);
    double* colorfieldgrad_i =
        container_i->GetPtrToParticleState(PARTICLEENGINE::ColorfieldGradient, particle_i);

    const double* mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);
    const double* dens_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j);
    double* colorfieldgrad_j =
        container_j->GetPtrToParticleState(PARTICLEENGINE::ColorfieldGradient, particle_j);

    const double fac =
        (UTILS::pow<2>(mass_i[0] / dens_i[0]) + UTILS::pow<2>(mass_j[0] / dens_j[0])) /
        (dens_i[0] + dens_j[0]);

    // sum contribution of neighboring particle j
    UTILS::vec_addscale(colorfieldgrad_i,
        fac * UTILS::pow<2>(dens_i[0]) / mass_i[0] * particlepair.dWdrij_, particlepair.e_ij_);

    // sum contribution of neighboring particle i
    if (status_j == PARTICLEENGINE::Owned)
      UTILS::vec_addscale(colorfieldgrad_j,
          -fac * UTILS::pow<2>(dens_j[0]) / mass_j[0] * particlepair.dWdrji_, particlepair.e_ij_);
  }
}

void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::ComputeUnitWallNormal() const
{
  // iterate over fluid particle types
  for (const auto& type_i : fluidtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear unit wall normal state
    container_i->ClearState(PARTICLEENGINE::UnitWallNormal);
  }

  // get relevant particle pair indices
  std::vector<int> relindices;
  neighborpairs_->GetRelevantParticlePairIndicesForDisjointCombination(
      boundarytypes_, fluidtypes_, relindices);

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

    // evaluate contribution of neighboring boundary particle j
    if (fluidtypes_.count(type_i))
    {
      // get container of owned particles of current particle type
      PARTICLEENGINE::ParticleContainer* container_i =
          particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

      // get material for current particle type
      const MAT::PAR::ParticleMaterialBase* material_i =
          particlematerial_->GetPtrToParticleMatParameter(type_i);

      // get pointer to particle states
      const double* mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
      const double* dens_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);
      double* wallnormal_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::UnitWallNormal, particle_i);

      // (current) inverse volume of particle i
      const double inv_V_i = dens_i[0] / mass_i[0];

      // (initial) volume of boundary particle j
      const double V_j = mass_i[0] / material_i->initDensity_;

      // sum contribution of neighboring boundary particle j
      UTILS::vec_addscale(
          wallnormal_i, -UTILS::pow<2>(V_j) * inv_V_i * particlepair.dWdrij_, particlepair.e_ij_);
    }

    // evaluate contribution of neighboring boundary particle i
    if (fluidtypes_.count(type_j) and status_j == PARTICLEENGINE::Owned)
    {
      // get container of owned particles of current particle type
      PARTICLEENGINE::ParticleContainer* container_j =
          particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

      // get material for current particle type
      const MAT::PAR::ParticleMaterialBase* material_j =
          particlematerial_->GetPtrToParticleMatParameter(type_j);

      // get pointer to particle states
      const double* mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);
      const double* dens_j =
          container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j);
      double* wallnormal_j =
          container_j->GetPtrToParticleState(PARTICLEENGINE::UnitWallNormal, particle_j);

      // (current) inverse volume of particle j
      const double inv_V_j = dens_j[0] / mass_j[0];

      // (initial) volume of boundary particle i
      const double V_i = mass_j[0] / material_j->initDensity_;

      // sum contribution of neighboring boundary particle i
      UTILS::vec_addscale(
          wallnormal_j, UTILS::pow<2>(V_i) * inv_V_j * particlepair.dWdrji_, particlepair.e_ij_);
    }
  }

  // iterate over fluid particle types
  for (const auto& type_i : fluidtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
    {
      // get pointer to particle state
      double* wallnormal_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::UnitWallNormal, particle_i);

      // norm of wall normal
      const double wallnormal_i_norm = UTILS::vec_norm2(wallnormal_i);

      // no interacting boundary particle
      if (not(wallnormal_i_norm > 0.0)) continue;

      // scale unit wall normal
      UTILS::vec_setscale(wallnormal_i, 1.0 / wallnormal_i_norm, wallnormal_i);
    }
  }
}

void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::ComputeWallDistance() const
{
  // iterate over fluid particle types
  for (const auto& type_i : fluidtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // set support radius of particles as initial wall distance
    container_i->UpdateState(0.0, PARTICLEENGINE::WallDistance, 1.0, PARTICLEENGINE::Radius);
  }

  // get relevant particle pair indices
  std::vector<int> relindices;
  neighborpairs_->GetRelevantParticlePairIndicesForDisjointCombination(
      boundarytypes_, fluidtypes_, relindices);

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

    // evaluate distance to neighboring boundary particle j
    if (fluidtypes_.count(type_i))
    {
      // get container of owned particles of current particle type
      PARTICLEENGINE::ParticleContainer* container_i =
          particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

      // get pointer to particle states
      const double* wallnormal_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::UnitWallNormal, particle_i);
      double* walldistance_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::WallDistance, particle_i);

      // no interacting boundary particle
      if (not(UTILS::vec_norm2(wallnormal_i) > 0.0)) continue;

      // distance of particle i to neighboring boundary particle j
      const double currentwalldistance =
          particlepair.absdist_ * UTILS::vec_dot(wallnormal_i, particlepair.e_ij_);

      // update wall distance of particle i
      walldistance_i[0] = std::min(walldistance_i[0], currentwalldistance);
    }

    // evaluate distance to neighboring boundary particle i
    if (fluidtypes_.count(type_j) and status_j == PARTICLEENGINE::Owned)
    {
      // get container of owned particles of current particle type
      PARTICLEENGINE::ParticleContainer* container_j =
          particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

      // get pointer to particle states
      const double* wallnormal_j =
          container_j->GetPtrToParticleState(PARTICLEENGINE::UnitWallNormal, particle_j);
      double* walldistance_j =
          container_j->GetPtrToParticleState(PARTICLEENGINE::WallDistance, particle_j);

      // no interacting boundary particle
      if (not(UTILS::vec_norm2(wallnormal_j) > 0.0)) continue;

      // distance of particle j to neighboring boundary particle i
      const double currentwalldistance =
          -particlepair.absdist_ * UTILS::vec_dot(wallnormal_j, particlepair.e_ij_);

      // update wall distance of particle j
      walldistance_j[0] = std::min(walldistance_j[0], currentwalldistance);
    }
  }
}

void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::ComputeWallCorrectionFactor()
{
  // get initial particle spacing
  const double initialparticlespacing = params_sph_.get<double>("INITIALPARTICLESPACING");

  // iterate over fluid particle types
  for (const auto& type_i : fluidtypes_)
  {
    // iterate over particle statuses
    for (auto& status_i : {PARTICLEENGINE::Owned, PARTICLEENGINE::Ghosted})
    {
      // get container of particles of current particle type and status
      PARTICLEENGINE::ParticleContainer* container_i =
          particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

      // get number of particles stored in container
      const int particlestored = container_i->ParticlesStored();

      // allocate memory
      f_i_[type_i][status_i].assign(particlestored, 0.0);

      // iterate over particles in container
      for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
      {
        // get pointer to particle states
        const double* rad_i =
            container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
        const double* walldistance_i =
            container_i->GetPtrToParticleState(PARTICLEENGINE::WallDistance, particle_i);

        // corrected wall distance and maximum correction distance
        const double dw_i = walldistance_i[0] - initialparticlespacing;
        const double dmax_i = kernel_->SmoothingLength(rad_i[0]);

        // determine correction factor
        double& f_i = f_i_[type_i][status_i][particle_i];
        f_i = 1.0;
        if (dw_i < 0.0)
          f_i = 0.0;
        else if (dw_i < dmax_i)
          f_i = dw_i / dmax_i;
      }
    }
  }
}

void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::
    ExtrapolateColorfieldGradientAtTriplePoint() const
{
  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  std::vector<std::vector<double>> sumj_fj_Vj_Wij(typevectorsize);
  std::vector<std::vector<std::vector<double>>> sumj_fj_Vj_Wij_CFGj(typevectorsize);

  // iterate over fluid particle types
  for (const auto& type_i : fluidtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container_i->ParticlesStored();

    // allocate memory
    sumj_fj_Vj_Wij[type_i].assign(particlestored, 0.0);
    sumj_fj_Vj_Wij_CFGj[type_i].assign(particlestored, std::vector<double>(3, 0.0));
  }

  // get relevant particle pair indices
  std::vector<int> relindices;
  neighborpairs_->GetRelevantParticlePairIndicesForEqualCombination(fluidtypes_, relindices);

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
    const double* rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
    const double* mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
    const double* dens_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);
    const double* colorfieldgrad_i =
        container_i->GetPtrToParticleState(PARTICLEENGINE::ColorfieldGradient, particle_i);
    const double* walldistance_i =
        container_i->GetPtrToParticleState(PARTICLEENGINE::WallDistance, particle_i);

    // get pointer to particle states
    const double* rad_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_j);
    const double* mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);
    const double* dens_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j);
    const double* colorfieldgrad_j =
        container_j->GetPtrToParticleState(PARTICLEENGINE::ColorfieldGradient, particle_j);
    const double* walldistance_j =
        container_j->GetPtrToParticleState(PARTICLEENGINE::WallDistance, particle_j);

    // change sign of colorfield gradient for different particle types
    double signfac = (type_i == type_j) ? 1.0 : -1.0;

    // get reference to correction factor
    const double& f_i = f_i_[type_i][status_i][particle_i];
    const double& f_j = f_i_[type_j][status_j][particle_j];

    // evaluate contribution of neighboring particle j
    if (walldistance_i[0] < rad_i[0] and f_i < 1.0 and f_j > 0.0 and
        UTILS::vec_norm2(colorfieldgrad_i) > 0.0)
    {
      const double fac = f_j * mass_j[0] / dens_j[0] * particlepair.Wij_;

      // initial estimate
      UTILS::vec_addscale(
          &sumj_fj_Vj_Wij_CFGj[type_i][particle_i][0], signfac * fac, colorfieldgrad_j);

      // correction factor
      sumj_fj_Vj_Wij[type_i][particle_i] += fac;
    }

    // evaluate contribution of neighboring particle i
    if (walldistance_j[0] < rad_j[0] and f_j < 1.0 and f_i > 0.0 and
        status_j == PARTICLEENGINE::Owned and UTILS::vec_norm2(colorfieldgrad_j) > 0.0)
    {
      const double fac = f_i * mass_i[0] / dens_i[0] * particlepair.Wji_;

      // initial estimate
      UTILS::vec_addscale(
          &sumj_fj_Vj_Wij_CFGj[type_j][particle_j][0], signfac * fac, colorfieldgrad_i);

      // correction factor
      sumj_fj_Vj_Wij[type_j][particle_j] += fac;
    }
  }

  // iterate over fluid particle types
  for (const auto& type_i : fluidtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
    {
      // evaluation only for particles with contributions from neighboring particles
      if (not(sumj_fj_Vj_Wij[type_i][particle_i] > 0.0) or
          not(UTILS::vec_norm2(&sumj_fj_Vj_Wij_CFGj[type_i][particle_i][0]) > 0.0))
        continue;

      // get pointer to particle states
      double* colorfieldgrad_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::ColorfieldGradient, particle_i);

      // norm of colorfield gradient
      const double colorfieldgrad_i_norm = UTILS::vec_norm2(colorfieldgrad_i);

      // get reference to correction factor
      const double& f_i = f_i_[type_i][PARTICLEENGINE::Owned][particle_i];

      // determine modified colorfield gradient
      UTILS::vec_setscale(colorfieldgrad_i, f_i, colorfieldgrad_i);
      const double fac = (1.0 - f_i) / sumj_fj_Vj_Wij[type_i][particle_i];
      UTILS::vec_addscale(colorfieldgrad_i, fac, &sumj_fj_Vj_Wij_CFGj[type_i][particle_i][0]);

      // scale modified colorfield gradient to magnitude of original colorfield gradient
      UTILS::vec_setscale(colorfieldgrad_i,
          colorfieldgrad_i_norm / UTILS::vec_norm2(colorfieldgrad_i), colorfieldgrad_i);
    }
  }
}

void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::ComputeInterfaceNormal() const
{
  // iterate over fluid particle types
  for (const auto& type_i : fluidtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear interface normal state
    container_i->ClearState(PARTICLEENGINE::InterfaceNormal);

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
    {
      // get pointer to particle states
      const double* rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
      const double* colorfieldgrad_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::ColorfieldGradient, particle_i);
      double* interfacenormal_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::InterfaceNormal, particle_i);

      // norm of colorfield gradient
      const double colorfieldgrad_i_norm = UTILS::vec_norm2(colorfieldgrad_i);

      // scale colorfield gradient
      if (colorfieldgrad_i_norm > (1.0e-10 * rad_i[0]))
        UTILS::vec_setscale(interfacenormal_i, 1.0 / colorfieldgrad_i_norm, colorfieldgrad_i);
    }
  }
}

void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::CorrectTriplePointNormal() const
{
  // iterate over fluid particle types
  for (const auto& type_i : fluidtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
    {
      // get pointer to particle states
      const double* rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
      const double* wallnormal_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::UnitWallNormal, particle_i);
      const double* walldistance_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::WallDistance, particle_i);
      double* interfacenormal_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::InterfaceNormal, particle_i);

      // evaluation only for particles close to boundary
      if (not(walldistance_i[0] < rad_i[0])) continue;

      // evaluation only for non-zero interface normal
      if (not(UTILS::vec_norm2(interfacenormal_i) > 0.0)) continue;

      // get reference to correction factor
      const double& f_i = f_i_[type_i][PARTICLEENGINE::Owned][particle_i];

      // no correction for current particle i
      if (f_i == 1.0) continue;

      // determine wall tangential
      double walltangential_i[3];
      UTILS::vec_set(walltangential_i, interfacenormal_i);
      UTILS::vec_addscale(
          walltangential_i, -UTILS::vec_dot(interfacenormal_i, wallnormal_i), wallnormal_i);

      // scale unit wall tangential
      UTILS::vec_setscale(
          walltangential_i, 1.0 / UTILS::vec_norm2(walltangential_i), walltangential_i);

      // convert static contact angle in radians
      const double theta_0 = (type_i == PARTICLEENGINE::Phase1)
                                 ? staticcontactangle_ * M_PI / 180.0
                                 : (180 - staticcontactangle_) * M_PI / 180.0;

      // determine triple point normal
      double triplepointnormal_i[3];
      UTILS::vec_setscale(triplepointnormal_i, std::sin(theta_0), walltangential_i);
      UTILS::vec_addscale(triplepointnormal_i, std::cos(theta_0), wallnormal_i);

      // determine corrected normal
      double correctednormal_i[3];
      UTILS::vec_setscale(correctednormal_i, f_i, interfacenormal_i);
      UTILS::vec_addscale(correctednormal_i, (1.0 - f_i), triplepointnormal_i);

      // overwrite interface normal with scaled corrected normal
      UTILS::vec_setscale(
          interfacenormal_i, 1.0 / UTILS::vec_norm2(correctednormal_i), correctednormal_i);
    }
  }
}

void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::
    ComputeSurfaceTensionContribution() const
{
  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  std::vector<std::vector<double>> interfacenormal_i_norm(typevectorsize);
  std::vector<std::vector<double>> sumj_nij_Vj_eij_dWij(typevectorsize);
  std::vector<std::vector<double>> sumj_Vj_Wij(typevectorsize);

  // iterate over fluid particle types
  for (const auto& type_i : fluidtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container_i->ParticlesStored();

    // allocate memory
    interfacenormal_i_norm[type_i].assign(particlestored, 0.0);
    sumj_nij_Vj_eij_dWij[type_i].assign(particlestored, 0.0);
    sumj_Vj_Wij[type_i].assign(particlestored, 0.0);

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
    {
      // get pointer to particle states
      const double* rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
      const double* mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
      const double* dens_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);
      const double* interfacenormal_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::InterfaceNormal, particle_i);

      // compute norm of interface normal
      interfacenormal_i_norm[type_i][particle_i] = UTILS::vec_norm2(interfacenormal_i);

      // evaluation only for non-zero interface normal
      if (not(interfacenormal_i_norm[type_i][particle_i] > 0.0)) continue;

      // evaluate kernel
      const double Wii = kernel_->W0(rad_i[0]);

      // add self-interaction
      sumj_Vj_Wij[type_i][particle_i] += Wii * mass_i[0] / dens_i[0];
    }
  }

  // get relevant particle pair indices
  std::vector<int> relindices;
  neighborpairs_->GetRelevantParticlePairIndicesForEqualCombination(fluidtypes_, relindices);

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
    const double* interfacenormal_i =
        container_i->GetPtrToParticleState(PARTICLEENGINE::InterfaceNormal, particle_i);

    const double* mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);
    const double* dens_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j);
    const double* interfacenormal_j =
        container_j->GetPtrToParticleState(PARTICLEENGINE::InterfaceNormal, particle_j);

    // evaluation only for non-zero interface normal
    if (not(interfacenormal_i_norm[type_i][particle_i] > 0.0)) continue;

    if (not((status_j == PARTICLEENGINE::Owned and
                interfacenormal_i_norm[type_j][particle_j] > 0.0) or
            UTILS::vec_norm2(interfacenormal_j) > 0.0))
      continue;

    // change sign of interface normal for different particle types
    double signfac = (type_i == type_j) ? 1.0 : -1.0;

    double n_ij[3];
    UTILS::vec_set(n_ij, interfacenormal_i);
    UTILS::vec_addscale(n_ij, -signfac, interfacenormal_j);

    const double fac = UTILS::vec_dot(n_ij, particlepair.e_ij_);

    // initial curvature estimate and correction factor
    const double V_j = mass_j[0] / dens_j[0];
    sumj_nij_Vj_eij_dWij[type_i][particle_i] += fac * V_j * particlepair.dWdrij_;
    sumj_Vj_Wij[type_i][particle_i] += V_j * particlepair.Wij_;

    if (status_j == PARTICLEENGINE::Owned)
    {
      const double V_i = mass_i[0] / dens_i[0];
      sumj_nij_Vj_eij_dWij[type_j][particle_j] += signfac * fac * V_i * particlepair.dWdrji_;
      sumj_Vj_Wij[type_j][particle_j] += V_i * particlepair.Wji_;
    }
  }

  // evaluate surface tension ramp function
  double timefac = 1.0;
  if (surfacetensionrampfctnumber_ > 0)
    timefac = DRT::Problem::Instance()->Funct(surfacetensionrampfctnumber_ - 1).EvaluateTime(time_);

  // iterate over fluid particle types
  for (const auto& type_i : fluidtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
    {
      // get pointer to particle states
      const double* rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
      const double* dens_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);
      const double* colorfieldgrad_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::ColorfieldGradient, particle_i);
      double* acc_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Acceleration, particle_i);

      const double* temp_i = (alphaT_ != 0.0) ? container_i->GetPtrToParticleState(
                                                    PARTICLEENGINE::Temperature, particle_i)
                                              : nullptr;

      // only add meaningful contributions
      if (std::abs(sumj_Vj_Wij[type_i][particle_i]) > (1.0e-10 * rad_i[0]))
      {
        // compute curvature
        const double curvature_i =
            -sumj_nij_Vj_eij_dWij[type_i][particle_i] / sumj_Vj_Wij[type_i][particle_i];

        // evaluate surface tension coefficient
        double alpha = alpha0_;
        if (alphaT_ != 0.0) alpha += alphaT_ * (temp_i[0] - reftemp_);

        // add contribution to acceleration
        UTILS::vec_addscale(acc_i, -timefac * alpha * curvature_i / dens_i[0], colorfieldgrad_i);
      }
    }
  }
}

void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::
    ComputeTempGradDrivenContribution() const
{
  // evaluate surface tension ramp function
  double timefac = 1.0;
  if (surfacetensionrampfctnumber_ > 0)
    timefac = DRT::Problem::Instance()->Funct(surfacetensionrampfctnumber_ - 1).EvaluateTime(time_);

  // iterate over fluid particle types
  for (const auto& type_i : fluidtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
    {
      // get pointer to particle states
      const double* dens_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);
      const double* colorfieldgrad_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::ColorfieldGradient, particle_i);
      const double* interfacenormal_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::InterfaceNormal, particle_i);
      const double* tempgrad_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::TemperatureGradient, particle_i);
      double* acc_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Acceleration, particle_i);

      // evaluation only for non-zero interface normal
      if (not(UTILS::vec_norm2(interfacenormal_i) > 0.0)) continue;

      // projection of temperature gradient onto tangential plane defined by interface normal
      double tempgrad_i_proj[3];
      UTILS::vec_set(tempgrad_i_proj, tempgrad_i);
      UTILS::vec_addscale(
          tempgrad_i_proj, -UTILS::vec_dot(tempgrad_i, interfacenormal_i), interfacenormal_i);

      // add contribution to acceleration
      UTILS::vec_addscale(acc_i, timefac * alphaT_ * UTILS::vec_norm2(colorfieldgrad_i) / dens_i[0],
          tempgrad_i_proj);
    }
  }
}
