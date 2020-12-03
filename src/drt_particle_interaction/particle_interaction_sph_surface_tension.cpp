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
PARTICLEINTERACTION::SPHSurfaceTension::SPHSurfaceTension(const Teuchos::ParameterList& params)
    : params_sph_(params),
      time_(0.0),
      surfacetensionrampfctnumber_(params.get<int>("SURFACETENSION_RAMP_FUNCT")),
      alpha0_(params_sph_.get<double>("SURFACETENSIONCOEFFICIENT")),
      alphamin_(params_sph_.get<double>("SURFACETENSIONMINIMUM")),
      staticcontactangle_(params_sph_.get<double>("STATICCONTACTANGLE")),
      alphaT_(params_sph_.get<double>("SURFACETENSIONTEMPFAC")),
      reftemp_(params_sph_.get<double>("SURFACETENSIONREFTEMP"))
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHSurfaceTension::Init()
{
  // init with potential fluid particle types
  fluidtypes_ = {PARTICLEENGINE::Phase1, PARTICLEENGINE::Phase2};

  // init with potential boundary particle types
  boundarytypes_ = {PARTICLEENGINE::BoundaryPhase, PARTICLEENGINE::RigidPhase};

  // safety check
  if (not(alpha0_ > 0.0)) dserror("constant factor of surface tension coefficient not positive!");

  if (not(alpha0_ > alphamin_))
    dserror("constant part smaller than minimum surface tension coefficient!");

  if (alphaT_ != 0.0)
  {
    if (DRT::INPUT::IntegralValue<INPAR::PARTICLE::TemperatureEvaluationScheme>(
            params_sph_, "TEMPERATUREEVALUATION") == INPAR::PARTICLE::NoTemperatureEvaluation)
      dserror("temperature evaluation needed for temperature dependent surface tension!");

    if (DRT::INPUT::IntegralValue<int>(params_sph_, "TEMPERATUREGRADIENT") == false)
      dserror("temperature gradient evaluation needed for temperature dependent surface tension!");
  }
}

void PARTICLEINTERACTION::SPHSurfaceTension::Setup(
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
  const auto fluidtypes = fluidtypes_;
  for (const auto& type_i : fluidtypes)
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i)) fluidtypes_.erase(type_i);

  // update with actual boundary particle types
  const auto boundarytypes = boundarytypes_;
  for (const auto& type_i : boundarytypes)
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i))
      boundarytypes_.erase(type_i);

  // setup interface normal of ghosted particles to refresh
  {
    std::vector<PARTICLEENGINE::StateEnum> states{PARTICLEENGINE::InterfaceNormal};

    // iterate over fluid particle types
    for (const auto& type_i : fluidtypes_)
      intnormtorefresh_.push_back(std::make_pair(type_i, states));
  }
}

void PARTICLEINTERACTION::SPHSurfaceTension::SetCurrentTime(const double currenttime)
{
  time_ = currenttime;
}

void PARTICLEINTERACTION::SPHSurfaceTension::InsertParticleStatesOfParticleTypes(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
    const
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
    particlestates.insert({PARTICLEENGINE::ColorfieldGradient, PARTICLEENGINE::InterfaceNormal,
        PARTICLEENGINE::Curvature});

    if (haveboundarytypes)
      particlestates.insert({PARTICLEENGINE::UnitWallNormal, PARTICLEENGINE::WallDistance});
  }
}

void PARTICLEINTERACTION::SPHSurfaceTension::ComputeInterfaceQuantities()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHSurfaceTension::ComputeInterfaceQuantities");

  // compute colorfield gradient
  ComputeColorfieldGradient();

  // compute interface normal
  ComputeInterfaceNormal();

  if (not boundarytypes_.empty())
  {
    // compute unit wall normal
    ComputeUnitWallNormal();

    // compute wall distance
    ComputeWallDistance();

    // correct normal vector of particles close to triple point
    CorrectTriplePointNormal();
  }

  // refresh interface normal
  particleengineinterface_->RefreshParticlesOfSpecificStatesAndTypes(intnormtorefresh_);
}

void PARTICLEINTERACTION::SPHSurfaceTension::AddAccelerationContribution()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHSurfaceTension::AddAccelerationContribution");

  // compute curvature
  ComputeCurvature();

  // compute surface tension contribution
  ComputeSurfaceTensionContribution();

  // compute temperature gradient driven contribution
  if (alphaT_ != 0.0) ComputeTempGradDrivenContribution();
}

void PARTICLEINTERACTION::SPHSurfaceTension::ComputeColorfieldGradient() const
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

    // (current) volume of particle i and j
    const double V_i = mass_i[0] / dens_i[0];
    const double V_j = mass_j[0] / dens_j[0];

    const double fac = (UTILS::pow<2>(V_i) + UTILS::pow<2>(V_j)) / (dens_i[0] + dens_j[0]);

    // sum contribution of neighboring particle j
    UTILS::vec_addscale(
        colorfieldgrad_i, dens_i[0] / V_i * fac * particlepair.dWdrij_, particlepair.e_ij_);

    // sum contribution of neighboring particle i
    if (status_j == PARTICLEENGINE::Owned)
      UTILS::vec_addscale(
          colorfieldgrad_j, -dens_j[0] / V_j * fac * particlepair.dWdrji_, particlepair.e_ij_);
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
      const double* rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
      double* colorfieldgrad_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::ColorfieldGradient, particle_i);

      // norm of colorfield gradient
      const double colorfieldgrad_i_norm = UTILS::vec_norm2(colorfieldgrad_i);

      // clear colorfield gradient
      if (not(colorfieldgrad_i_norm > (1.0e-10 * rad_i[0]))) UTILS::vec_clear(colorfieldgrad_i);
    }
  }
}

void PARTICLEINTERACTION::SPHSurfaceTension::ComputeInterfaceNormal() const
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

      // set interface normal
      if (colorfieldgrad_i_norm > (1.0e-10 * rad_i[0]))
        UTILS::vec_setscale(interfacenormal_i, 1.0 / colorfieldgrad_i_norm, colorfieldgrad_i);
    }
  }
}

void PARTICLEINTERACTION::SPHSurfaceTension::ComputeUnitWallNormal() const
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

    // get corresponding particle containers
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    PARTICLEENGINE::ParticleContainer* container_j =
        particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

    // evaluate contribution of neighboring boundary particle j
    if (fluidtypes_.count(type_i))
    {
      // get material for current particle type
      const MAT::PAR::ParticleMaterialBase* material_j =
          particlematerial_->GetPtrToParticleMatParameter(type_j);

      // get pointer to particle states
      const double* mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
      const double* dens_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);
      double* wallnormal_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::UnitWallNormal, particle_i);

      const double* mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);

      // (current) volume of particle i
      const double V_i = mass_i[0] / dens_i[0];

      // (initial) volume of boundary particle j
      const double V_j = mass_j[0] / material_j->initDensity_;

      const double fac =
          (UTILS::pow<2>(V_i) + UTILS::pow<2>(V_j)) / (dens_i[0] + material_j->initDensity_);

      // sum contribution of neighboring boundary particle j
      UTILS::vec_addscale(
          wallnormal_i, dens_i[0] / V_i * fac * particlepair.dWdrij_, particlepair.e_ij_);
    }

    // evaluate contribution of neighboring boundary particle i
    if (fluidtypes_.count(type_j) and status_j == PARTICLEENGINE::Owned)
    {
      // get material for current particle type
      const MAT::PAR::ParticleMaterialBase* material_i =
          particlematerial_->GetPtrToParticleMatParameter(type_i);

      // get pointer to particle states
      const double* mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);
      const double* dens_j =
          container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j);
      double* wallnormal_j =
          container_j->GetPtrToParticleState(PARTICLEENGINE::UnitWallNormal, particle_j);

      const double* mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);

      // (initial) volume of boundary particle i
      const double V_i = mass_i[0] / material_i->initDensity_;

      // (current) volume of particle j
      const double V_j = mass_j[0] / dens_j[0];

      const double fac =
          (UTILS::pow<2>(V_i) + UTILS::pow<2>(V_j)) / (material_i->initDensity_ + dens_j[0]);

      // sum contribution of neighboring boundary particle i
      UTILS::vec_addscale(
          wallnormal_j, -dens_j[0] / V_j * fac * particlepair.dWdrij_, particlepair.e_ij_);
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
      const double* rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
      double* wallnormal_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::UnitWallNormal, particle_i);

      // norm of wall normal
      const double wallnormal_i_norm = UTILS::vec_norm2(wallnormal_i);

      // scale or clear unit wall normal
      if (wallnormal_i_norm > (1.0e-10 * rad_i[0]))
        UTILS::vec_setscale(wallnormal_i, 1.0 / wallnormal_i_norm, wallnormal_i);
      else
        UTILS::vec_clear(wallnormal_i);
    }
  }
}

void PARTICLEINTERACTION::SPHSurfaceTension::ComputeWallDistance() const
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
          -particlepair.absdist_ * UTILS::vec_dot(wallnormal_i, particlepair.e_ij_);

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
          particlepair.absdist_ * UTILS::vec_dot(wallnormal_j, particlepair.e_ij_);

      // update wall distance of particle j
      walldistance_j[0] = std::min(walldistance_j[0], currentwalldistance);
    }
  }
}

void PARTICLEINTERACTION::SPHSurfaceTension::CorrectTriplePointNormal() const
{
  // get initial particle spacing
  const double initialparticlespacing = params_sph_.get<double>("INITIALPARTICLESPACING");

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

      // corrected wall distance and maximum correction distance
      const double dw_i = walldistance_i[0] - initialparticlespacing;
      const double dmax_i = kernel_->SmoothingLength(rad_i[0]);

      // determine correction factor
      double f_i = 1.0;
      if (dw_i < 0.0)
        f_i = 0.0;
      else if (dw_i < dmax_i)
        f_i = dw_i / dmax_i;

      // no correction for current particle i
      if (f_i == 1.0) continue;

      // determine wall tangential
      double walltangential_i[3];
      UTILS::vec_set(walltangential_i, interfacenormal_i);
      UTILS::vec_addscale(
          walltangential_i, -UTILS::vec_dot(interfacenormal_i, wallnormal_i), wallnormal_i);

      // norm of wall tangential
      const double walltangential_i_norm = UTILS::vec_norm2(walltangential_i);

      // scale or clear unit wall tangential
      if (walltangential_i_norm > (1.0e-10 * rad_i[0]))
        UTILS::vec_setscale(walltangential_i, 1.0 / walltangential_i_norm, walltangential_i);
      else
        UTILS::vec_clear(walltangential_i);

      // convert static contact angle in radians
      const double theta_0 = (type_i == PARTICLEENGINE::Phase1)
                                 ? staticcontactangle_ * M_PI / 180.0
                                 : (180 - staticcontactangle_) * M_PI / 180.0;

      // determine triple point normal
      double triplepointnormal_i[3];
      UTILS::vec_setscale(triplepointnormal_i, std::sin(theta_0), walltangential_i);
      UTILS::vec_addscale(triplepointnormal_i, -std::cos(theta_0), wallnormal_i);

      // determine corrected normal
      double correctednormal_i[3];
      UTILS::vec_setscale(correctednormal_i, f_i, interfacenormal_i);
      UTILS::vec_addscale(correctednormal_i, (1.0 - f_i), triplepointnormal_i);

      // norm of corrected normal
      const double correctednormal_i_norm = UTILS::vec_norm2(correctednormal_i);

      // scale or clear interface normal
      if (correctednormal_i_norm > (1.0e-10 * rad_i[0]))
        UTILS::vec_setscale(interfacenormal_i, 1.0 / correctednormal_i_norm, correctednormal_i);
      else
        UTILS::vec_clear(interfacenormal_i);

      // overwrite interface normal with scaled corrected normal
      UTILS::vec_setscale(
          interfacenormal_i, 1.0 / UTILS::vec_norm2(correctednormal_i), correctednormal_i);
    }
  }
}

void PARTICLEINTERACTION::SPHSurfaceTension::ComputeCurvature() const
{
  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  std::vector<std::vector<double>> sumj_nij_Vj_eij_dWij(typevectorsize);
  std::vector<std::vector<double>> sumj_Vj_Wij(typevectorsize);

  // iterate over fluid particle types
  for (const auto& type_i : fluidtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear curvature state
    container_i->ClearState(PARTICLEENGINE::Curvature);

    // get number of particles stored in container
    const int particlestored = container_i->ParticlesStored();

    // allocate memory
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

      // evaluation only for non-zero interface normal
      if (not(UTILS::vec_norm2(interfacenormal_i) > 0.0)) continue;

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

    // evaluation only for non-zero interface normals
    if (not(UTILS::vec_norm2(interfacenormal_i) > 0.0) or
        not(UTILS::vec_norm2(interfacenormal_j) > 0.0))
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
      const double* interfacenormal_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::InterfaceNormal, particle_i);
      double* curvature_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::Curvature, particle_i);

      // evaluation only for non-zero interface normal
      if (not(UTILS::vec_norm2(interfacenormal_i) > 0.0)) continue;

      // compute curvature
      curvature_i[0] = -sumj_nij_Vj_eij_dWij[type_i][particle_i] / sumj_Vj_Wij[type_i][particle_i];
    }
  }
}

void PARTICLEINTERACTION::SPHSurfaceTension::ComputeSurfaceTensionContribution() const
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
      const double* curvature_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::Curvature, particle_i);
      const double* colorfieldgrad_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::ColorfieldGradient, particle_i);
      double* acc_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Acceleration, particle_i);

      const double* temp_i = (alphaT_ != 0.0) ? container_i->GetPtrToParticleState(
                                                    PARTICLEENGINE::Temperature, particle_i)
                                              : nullptr;

      // evaluate surface tension coefficient
      double alpha = alpha0_;
      if (alphaT_ != 0.0)
      {
        alpha += alphaT_ * (temp_i[0] - reftemp_);
        alpha = std::max(alpha, alphamin_);
      }

      // add contribution to acceleration
      UTILS::vec_addscale(acc_i, -timefac * alpha * curvature_i[0] / dens_i[0], colorfieldgrad_i);
    }
  }
}

void PARTICLEINTERACTION::SPHSurfaceTension::ComputeTempGradDrivenContribution() const
{
  // evaluate surface tension ramp function
  double timefac = 1.0;
  if (surfacetensionrampfctnumber_ > 0)
    timefac = DRT::Problem::Instance()->Funct(surfacetensionrampfctnumber_ - 1).EvaluateTime(time_);

  // temperature in transition from linear to constant regime of surface tension coefficient
  const double transitiontemp = reftemp_ + (alphamin_ - alpha0_) / alphaT_;

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
      const double* temp_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::Temperature, particle_i);
      const double* tempgrad_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::TemperatureGradient, particle_i);
      double* acc_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Acceleration, particle_i);

      // evaluation only for non-zero interface normal
      if (not(UTILS::vec_norm2(interfacenormal_i) > 0.0)) continue;

      // no evaluation in the regime of constant surface tension coefficient
      if (temp_i[0] > transitiontemp) continue;

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
