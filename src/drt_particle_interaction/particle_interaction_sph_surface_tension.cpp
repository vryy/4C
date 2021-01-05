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
#include "particle_interaction_sph_surface_tension_interface_viscosity.H"
#include "particle_interaction_sph_surface_tension_recoilpressure_evaporation.H"
#include "particle_interaction_sph_surface_tension_barrier_force.H"

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
      liquidtype_(PARTICLEENGINE::Phase1),
      gastype_(PARTICLEENGINE::Phase2),
      time_(0.0),
      timerampfct_(params.get<int>("SURFACETENSION_RAMP_FUNCT")),
      alpha0_(params_sph_.get<double>("SURFACETENSIONCOEFFICIENT")),
      alphamin_(params_sph_.get<double>("SURFACETENSIONMINIMUM")),
      staticcontactangle_(params_sph_.get<double>("STATICCONTACTANGLE")),
      alphaT_(params_sph_.get<double>("SURFACETENSIONTEMPFAC")),
      surf_ref_temp_(params_sph_.get<double>("SURFACETENSIONREFTEMP"))
{
  // empty constructor
}

PARTICLEINTERACTION::SPHSurfaceTension::~SPHSurfaceTension() = default;

void PARTICLEINTERACTION::SPHSurfaceTension::Init()
{
  // init interface viscosity handler
  InitInterfaceViscosityHandler();

  // init evaporation induced recoil pressure handler
  InitRecoilPressureEvaporationHandler();

  // init barrier force handler
  InitBarrierForceHandler();

  // init fluid particle types
  fluidtypes_ = {liquidtype_, gastype_};

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

  // set neighbor pair handler
  neighborpairs_ = neighborpairs;

  // setup interface viscosity handler
  if (interfaceviscosity_)
    interfaceviscosity_->Setup(
        particleengineinterface, kernel, particlematerial, equationofstatebundle, neighborpairs);

  // setup evaporation induced recoil pressure handler
  if (recoilpressureevaporation_) recoilpressureevaporation_->Setup(particleengineinterface);

  // setup barrier force handler
  if (barrierforce_) barrierforce_->Setup(particleengineinterface, neighborpairs);

  // safety check
  for (const auto& type_i : fluidtypes_)
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i))
      dserror("no particle container for particle type '%s' found!",
          PARTICLEENGINE::EnumToTypeName(type_i).c_str());

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
      particlestates.insert({PARTICLEENGINE::WallColorfield, PARTICLEENGINE::WallInterfaceNormal});
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
    // compute wall colorfield and wall interface normal
    ComputeWallColorfieldAndWallInterfaceNormal();

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

  // compute interface viscosity contribution
  if (interfaceviscosity_) interfaceviscosity_->ComputeInterfaceViscosityContribution();

  // compute evaporation induced recoil pressure contribution
  if (recoilpressureevaporation_) recoilpressureevaporation_->ComputeRecoilPressureContribution();

  // compute barrier force contribution
  if (barrierforce_) barrierforce_->ComputeBarrierForceContribution();
}

void PARTICLEINTERACTION::SPHSurfaceTension::InitInterfaceViscosityHandler()
{
  // create interface viscosity handler
  if (DRT::INPUT::IntegralValue<int>(params_sph_, "INTERFACE_VISCOSITY"))
    interfaceviscosity_ = std::unique_ptr<PARTICLEINTERACTION::SPHInterfaceViscosity>(
        new PARTICLEINTERACTION::SPHInterfaceViscosity(params_sph_));

  // init interface viscosity handler
  if (interfaceviscosity_) interfaceviscosity_->Init();
}

void PARTICLEINTERACTION::SPHSurfaceTension::InitRecoilPressureEvaporationHandler()
{
  // create evaporation induced recoil pressure handler
  if (DRT::INPUT::IntegralValue<int>(params_sph_, "VAPOR_RECOIL"))
    recoilpressureevaporation_ = std::unique_ptr<PARTICLEINTERACTION::SPHRecoilPressureEvaporation>(
        new PARTICLEINTERACTION::SPHRecoilPressureEvaporation(params_sph_));

  // init evaporation induced recoil pressure handler
  if (recoilpressureevaporation_) recoilpressureevaporation_->Init();
}

void PARTICLEINTERACTION::SPHSurfaceTension::InitBarrierForceHandler()
{
  // create barrier force handler
  if (DRT::INPUT::IntegralValue<int>(params_sph_, "BARRIER_FORCE"))
    barrierforce_ = std::unique_ptr<PARTICLEINTERACTION::SPHBarrierForce>(
        new PARTICLEINTERACTION::SPHBarrierForce(params_sph_));

  // init barrier force handler
  if (barrierforce_) barrierforce_->Init();
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
    const double* mass_i = container_i->GetPtrToState(PARTICLEENGINE::Mass, particle_i);
    const double* dens_i = container_i->GetPtrToState(PARTICLEENGINE::Density, particle_i);
    double* cfg_i = container_i->GetPtrToState(PARTICLEENGINE::ColorfieldGradient, particle_i);

    const double* mass_j = container_j->GetPtrToState(PARTICLEENGINE::Mass, particle_j);
    const double* dens_j = container_j->GetPtrToState(PARTICLEENGINE::Density, particle_j);
    double* cfg_j = container_j->GetPtrToState(PARTICLEENGINE::ColorfieldGradient, particle_j);

    // (current) volume of particle i and j
    const double V_i = mass_i[0] / dens_i[0];
    const double V_j = mass_j[0] / dens_j[0];

    const double fac = (UTILS::pow<2>(V_i) + UTILS::pow<2>(V_j)) / (dens_i[0] + dens_j[0]);

    // sum contribution of neighboring particle j
    UTILS::vec_addscale(cfg_i, dens_i[0] / V_i * fac * particlepair.dWdrij_, particlepair.e_ij_);

    // sum contribution of neighboring particle i
    if (status_j == PARTICLEENGINE::Owned)
      UTILS::vec_addscale(cfg_j, -dens_j[0] / V_j * fac * particlepair.dWdrji_, particlepair.e_ij_);
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
      const double* rad_i = container_i->GetPtrToState(PARTICLEENGINE::Radius, particle_i);
      double* cfg_i = container_i->GetPtrToState(PARTICLEENGINE::ColorfieldGradient, particle_i);

      // norm of colorfield gradient
      const double cfg_i_norm = UTILS::vec_norm2(cfg_i);

      // clear colorfield gradient
      if (not(cfg_i_norm > (1.0e-10 * rad_i[0]))) UTILS::vec_clear(cfg_i);
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
      const double* rad_i = container_i->GetPtrToState(PARTICLEENGINE::Radius, particle_i);
      const double* cfg_i =
          container_i->GetPtrToState(PARTICLEENGINE::ColorfieldGradient, particle_i);
      double* ifn_i = container_i->GetPtrToState(PARTICLEENGINE::InterfaceNormal, particle_i);

      // norm of colorfield gradient
      const double cfg_i_norm = UTILS::vec_norm2(cfg_i);

      // set interface normal
      if (cfg_i_norm > (1.0e-10 * rad_i[0])) UTILS::vec_setscale(ifn_i, 1.0 / cfg_i_norm, cfg_i);
    }
  }
}

void PARTICLEINTERACTION::SPHSurfaceTension::ComputeWallColorfieldAndWallInterfaceNormal() const
{
  // iterate over fluid particle types
  for (const auto& type_i : fluidtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear wall colorfield state
    container_i->ClearState(PARTICLEENGINE::WallColorfield);

    // clear wall interface normal state
    container_i->ClearState(PARTICLEENGINE::WallInterfaceNormal);
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
      const double* mass_i = container_i->GetPtrToState(PARTICLEENGINE::Mass, particle_i);
      const double* dens_i = container_i->GetPtrToState(PARTICLEENGINE::Density, particle_i);
      double* wallcf_i = container_i->GetPtrToState(PARTICLEENGINE::WallColorfield, particle_i);
      double* wallifn_i =
          container_i->GetPtrToState(PARTICLEENGINE::WallInterfaceNormal, particle_i);

      const double* mass_j = container_j->GetPtrToState(PARTICLEENGINE::Mass, particle_j);

      // (current) volume of particle i
      const double V_i = mass_i[0] / dens_i[0];

      // (initial) volume of boundary particle j
      const double V_j = mass_j[0] / material_j->initDensity_;

      const double fac = (UTILS::pow<2>(V_i) + UTILS::pow<2>(V_j)) * dens_i[0] /
                         (V_i * (dens_i[0] + material_j->initDensity_));

      // sum contribution of neighboring boundary particle j
      wallcf_i[0] += fac * particlepair.Wij_;
      UTILS::vec_addscale(wallifn_i, fac * particlepair.dWdrij_, particlepair.e_ij_);
    }

    // evaluate contribution of neighboring boundary particle i
    if (fluidtypes_.count(type_j) and status_j == PARTICLEENGINE::Owned)
    {
      // get material for current particle type
      const MAT::PAR::ParticleMaterialBase* material_i =
          particlematerial_->GetPtrToParticleMatParameter(type_i);

      // get pointer to particle states
      const double* mass_j = container_j->GetPtrToState(PARTICLEENGINE::Mass, particle_j);
      const double* dens_j = container_j->GetPtrToState(PARTICLEENGINE::Density, particle_j);
      double* wallcf_j = container_j->GetPtrToState(PARTICLEENGINE::WallColorfield, particle_j);
      double* wallifn_j =
          container_j->GetPtrToState(PARTICLEENGINE::WallInterfaceNormal, particle_j);

      const double* mass_i = container_i->GetPtrToState(PARTICLEENGINE::Mass, particle_i);

      // (initial) volume of boundary particle i
      const double V_i = mass_i[0] / material_i->initDensity_;

      // (current) volume of particle j
      const double V_j = mass_j[0] / dens_j[0];

      const double fac = (UTILS::pow<2>(V_i) + UTILS::pow<2>(V_j)) * dens_j[0] /
                         (V_j * (material_i->initDensity_ + dens_j[0]));

      // sum contribution of neighboring boundary particle i
      wallcf_j[0] += fac * particlepair.Wji_;
      UTILS::vec_addscale(wallifn_j, -fac * particlepair.dWdrji_, particlepair.e_ij_);
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
      const double* rad_i = container_i->GetPtrToState(PARTICLEENGINE::Radius, particle_i);
      double* wallifn_i =
          container_i->GetPtrToState(PARTICLEENGINE::WallInterfaceNormal, particle_i);

      // norm of wall interface normal
      const double wallifn_i_norm = UTILS::vec_norm2(wallifn_i);

      // scale or clear wall interface normal
      if (wallifn_i_norm > (1.0e-10 * rad_i[0]))
        UTILS::vec_setscale(wallifn_i, 1.0 / wallifn_i_norm, wallifn_i);
      else
        UTILS::vec_clear(wallifn_i);
    }
  }
}

void PARTICLEINTERACTION::SPHSurfaceTension::CorrectTriplePointNormal() const
{
  // iterate over fluid particle types
  for (const auto& type_i : fluidtypes_)
  {
    // static contact angle with respect to liquid particle type
    const double staticcontactangle =
        (type_i == liquidtype_) ? staticcontactangle_ : (180 - staticcontactangle_);

    // convert static contact angle in radians
    const double theta_0 = staticcontactangle * M_PI / 180.0;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
    {
      // get pointer to particle states
      const double* rad_i = container_i->GetPtrToState(PARTICLEENGINE::Radius, particle_i);
      const double* wallifn_i =
          container_i->GetPtrToState(PARTICLEENGINE::WallInterfaceNormal, particle_i);
      const double* wallcf_i =
          container_i->GetPtrToState(PARTICLEENGINE::WallColorfield, particle_i);
      double* ifn_i = container_i->GetPtrToState(PARTICLEENGINE::InterfaceNormal, particle_i);

      // evaluation only for non-zero wall interface normal
      if (not(UTILS::vec_norm2(wallifn_i) > 0.0)) continue;

      // evaluation only for non-zero interface normal
      if (not(UTILS::vec_norm2(ifn_i) > 0.0)) continue;

      // determine correction factor
      double f_i = UTILS::complintrans(wallcf_i[0], 0.0, 0.2);

      // determine wall interface tangential
      double wallift_i[3];
      UTILS::vec_set(wallift_i, ifn_i);
      UTILS::vec_addscale(wallift_i, -UTILS::vec_dot(ifn_i, wallifn_i), wallifn_i);

      // norm of wall interface tangential
      const double wallift_i_norm = UTILS::vec_norm2(wallift_i);

      // scale or clear wall interface tangential
      if (wallift_i_norm > (1.0e-10 * rad_i[0]))
        UTILS::vec_setscale(wallift_i, 1.0 / wallift_i_norm, wallift_i);
      else
        UTILS::vec_clear(wallift_i);

      // determine triple point normal
      double tpn_i[3];
      UTILS::vec_setscale(tpn_i, std::sin(theta_0), wallift_i);
      UTILS::vec_addscale(tpn_i, -std::cos(theta_0), wallifn_i);

      // determine corrected interface normal
      double corifn_i[3];
      UTILS::vec_setscale(corifn_i, f_i, ifn_i);
      UTILS::vec_addscale(corifn_i, (1.0 - f_i), tpn_i);

      // norm of corrected interface normal
      const double corifn_i_norm = UTILS::vec_norm2(corifn_i);

      // scale or clear interface normal
      if (corifn_i_norm > (1.0e-10 * rad_i[0]))
        UTILS::vec_setscale(ifn_i, 1.0 / corifn_i_norm, corifn_i);
      else
        UTILS::vec_clear(ifn_i);
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
      const double* rad_i = container_i->GetPtrToState(PARTICLEENGINE::Radius, particle_i);
      const double* mass_i = container_i->GetPtrToState(PARTICLEENGINE::Mass, particle_i);
      const double* dens_i = container_i->GetPtrToState(PARTICLEENGINE::Density, particle_i);
      const double* ifn_i = container_i->GetPtrToState(PARTICLEENGINE::InterfaceNormal, particle_i);

      // evaluation only for non-zero interface normal
      if (not(UTILS::vec_norm2(ifn_i) > 0.0)) continue;

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
    const double* mass_i = container_i->GetPtrToState(PARTICLEENGINE::Mass, particle_i);
    const double* dens_i = container_i->GetPtrToState(PARTICLEENGINE::Density, particle_i);
    const double* ifn_i = container_i->GetPtrToState(PARTICLEENGINE::InterfaceNormal, particle_i);

    const double* mass_j = container_j->GetPtrToState(PARTICLEENGINE::Mass, particle_j);
    const double* dens_j = container_j->GetPtrToState(PARTICLEENGINE::Density, particle_j);
    const double* ifn_j = container_j->GetPtrToState(PARTICLEENGINE::InterfaceNormal, particle_j);

    // evaluation only for non-zero interface normals
    if (not(UTILS::vec_norm2(ifn_i) > 0.0) or not(UTILS::vec_norm2(ifn_j) > 0.0)) continue;

    // change sign of interface normal for different particle types
    double signfac = (type_i == type_j) ? 1.0 : -1.0;

    double n_ij[3];
    UTILS::vec_set(n_ij, ifn_i);
    UTILS::vec_addscale(n_ij, -signfac, ifn_j);

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
      const double* ifn_i = container_i->GetPtrToState(PARTICLEENGINE::InterfaceNormal, particle_i);
      double* curv_i = container_i->GetPtrToState(PARTICLEENGINE::Curvature, particle_i);

      // evaluation only for non-zero interface normal
      if (not(UTILS::vec_norm2(ifn_i) > 0.0)) continue;

      // compute curvature
      curv_i[0] = -sumj_nij_Vj_eij_dWij[type_i][particle_i] / sumj_Vj_Wij[type_i][particle_i];
    }
  }
}

void PARTICLEINTERACTION::SPHSurfaceTension::ComputeSurfaceTensionContribution() const
{
  // evaluate surface tension time ramp function
  double timefac = 1.0;
  if (timerampfct_ > 0)
    timefac = DRT::Problem::Instance()->Funct(timerampfct_ - 1).EvaluateTime(time_);

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
      const double* dens_i = container_i->GetPtrToState(PARTICLEENGINE::Density, particle_i);
      const double* curv_i = container_i->GetPtrToState(PARTICLEENGINE::Curvature, particle_i);
      const double* cfg_i =
          container_i->GetPtrToState(PARTICLEENGINE::ColorfieldGradient, particle_i);
      const double* ifn_i = container_i->GetPtrToState(PARTICLEENGINE::InterfaceNormal, particle_i);
      double* acc_i = container_i->GetPtrToState(PARTICLEENGINE::Acceleration, particle_i);

      const double* temp_i =
          (alphaT_ != 0.0) ? container_i->GetPtrToState(PARTICLEENGINE::Temperature, particle_i)
                           : nullptr;

      // evaluation only for non-zero interface normal
      if (not(UTILS::vec_norm2(ifn_i) > 0.0)) continue;

      // evaluate surface tension coefficient
      double alpha = alpha0_;
      if (alphaT_ != 0.0)
      {
        alpha += alphaT_ * (temp_i[0] - surf_ref_temp_);
        alpha = std::max(alpha, alphamin_);
      }

      // add contribution to acceleration
      UTILS::vec_addscale(acc_i, -timefac * alpha * curv_i[0] / dens_i[0], cfg_i);
    }
  }
}

void PARTICLEINTERACTION::SPHSurfaceTension::ComputeTempGradDrivenContribution() const
{
  // evaluate surface tension time ramp function
  double timefac = 1.0;
  if (timerampfct_ > 0)
    timefac = DRT::Problem::Instance()->Funct(timerampfct_ - 1).EvaluateTime(time_);

  // temperature in transition from linear to constant regime of surface tension coefficient
  const double transitiontemp = surf_ref_temp_ + (alphamin_ - alpha0_) / alphaT_;

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
      const double* dens_i = container_i->GetPtrToState(PARTICLEENGINE::Density, particle_i);
      const double* cfg_i =
          container_i->GetPtrToState(PARTICLEENGINE::ColorfieldGradient, particle_i);
      const double* ifn_i = container_i->GetPtrToState(PARTICLEENGINE::InterfaceNormal, particle_i);
      const double* temp_i = container_i->GetPtrToState(PARTICLEENGINE::Temperature, particle_i);
      const double* tempgrad_i =
          container_i->GetPtrToState(PARTICLEENGINE::TemperatureGradient, particle_i);
      double* acc_i = container_i->GetPtrToState(PARTICLEENGINE::Acceleration, particle_i);

      // evaluation only for non-zero interface normal
      if (not(UTILS::vec_norm2(ifn_i) > 0.0)) continue;

      // no evaluation in the regime of constant surface tension coefficient
      if (temp_i[0] > transitiontemp) continue;

      // projection of temperature gradient onto tangential plane defined by interface normal
      double tempgrad_i_proj[3];
      UTILS::vec_set(tempgrad_i_proj, tempgrad_i);
      UTILS::vec_addscale(tempgrad_i_proj, -UTILS::vec_dot(tempgrad_i, ifn_i), ifn_i);

      // add contribution to acceleration
      UTILS::vec_addscale(
          acc_i, timefac * alphaT_ * UTILS::vec_norm2(cfg_i) / dens_i[0], tempgrad_i_proj);
    }
  }
}
