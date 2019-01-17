/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_density.cpp

\brief density handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_density.H"

#include "particle_interaction_sph_kernel.H"
#include "particle_interaction_material_handler.H"
#include "particle_interaction_sph_equationofstate.H"
#include "particle_interaction_sph_equationofstate_bundle.H"
#include "particle_interaction_sph_neighbor_pairs.H"
#include "particle_interaction_sph_density_correction.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHDensityBase::SPHDensityBase(const Teuchos::ParameterList& params)
    : params_sph_(params), dt_(0.0), computecolorfield_(false), applytransportvelocity_(false)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init density handler                                       sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityBase::Init()
{
  // check if transport velocity formulation is applied
  if (DRT::INPUT::IntegralValue<INPAR::PARTICLE::TransportVelocityFormulation>(
          params_sph_, "TRANSPORTVELOCITYFORMULATION") !=
      INPAR::PARTICLE::TransportVelocityFormulation::NoTransportVelocity)
    applytransportvelocity_ = true;
}

/*---------------------------------------------------------------------------*
 | setup density handler                                      sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityBase::Setup(
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

  // setup density of ghosted particles to refresh
  {
    std::vector<PARTICLEENGINE::StateEnum> states{PARTICLEENGINE::Density};

    for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
    {
      // no refreshing of density states for boundary or rigid particles
      if (typeEnum == PARTICLEENGINE::BoundaryPhase or typeEnum == PARTICLEENGINE::RigidPhase)
        continue;

      densitytorefresh_.push_back(std::make_pair(typeEnum, states));
    }
  }
}

/*---------------------------------------------------------------------------*
 | write restart of density handler                           sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityBase::WriteRestart(const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of density handler                            sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | set current step size                                      sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityBase::SetCurrentStepSize(const double currentstepsize)
{
  dt_ = currentstepsize;
}

/*---------------------------------------------------------------------------*
 | evaluate sum of weighted mass and colorfield               sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityBase::SumWeightedMassAndColorfield() const
{
  // get reference to neighbor pair data
  const SPHNeighborPairData& neighborpairdata = neighborpairs_->GetRefToNeighborPairData();

  // iterate over particle types
  for (auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // no density summation for boundary or rigid particles
    if (type_i == PARTICLEENGINE::BoundaryPhase or type_i == PARTICLEENGINE::RigidPhase) continue;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get material for current particle type
    const MAT::PAR::ParticleMaterialBase* material_i =
        particlematerial_->GetPtrToParticleMatParameter(type_i);

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
    {
      // get reference to vector of neighbor pairs of current particle
      const std::vector<SPHNeighborPair>& currentNeighborPairs =
          (neighborpairdata[type_i])[particle_i];

      // check for neighbor pairs of current particle
      if (currentNeighborPairs.empty()) continue;

      // declare pointer variables for particle i
      const double *rad_i, *mass_i, *dens_i;
      double *denssum_i, *colorfield_i;

      // get pointer to particle states
      rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
      mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
      denssum_i = container_i->GetPtrToParticleState(PARTICLEENGINE::DensitySum, particle_i);

      if (computecolorfield_)
      {
        dens_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);
        colorfield_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Colorfield, particle_i);
      }

      // initialize sum of evaluated kernel values for particle i due to neighbor particles j
      double sumj_Wij(0.0);
      double sumj_Wij_dens_j_mass_j(0.0);

      // iterate over neighbor pairs
      for (auto& neighborIt : currentNeighborPairs)
      {
        // access values of local index tuple of neighboring particle
        PARTICLEENGINE::TypeEnum type_j;
        PARTICLEENGINE::StatusEnum status_j;
        int particle_j;
        std::tie(type_j, status_j, particle_j) = neighborIt.first;

        // get reference to current particle pair
        const SPHParticlePair& particlepair = neighborIt.second;

        // get container of neighboring particles of current particle type and state
        PARTICLEENGINE::ParticleContainerShrdPtr container_j =
            particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

        // declare pointer variables for neighbor particle j
        const double *mass_j, *dens_j;

        // get pointer to particle states
        if (computecolorfield_)
        {
          mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);

          if (type_j == PARTICLEENGINE::BoundaryPhase or type_j == PARTICLEENGINE::RigidPhase)
            dens_j = &(material_i->initDensity_);
          else
            dens_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j);
        }

        // sum contribution of neighbor particle j
        sumj_Wij += particlepair.Wij_;

        if (computecolorfield_)
          sumj_Wij_dens_j_mass_j += (particlepair.Wij_ / dens_j[0]) * mass_j[0];
      }

      // evaluate kernel
      const double Wii = kernel_->W(0.0, rad_i[0]);

      // add self-interaction and contributions of neighbor particles
      denssum_i[0] = (Wii + sumj_Wij) * mass_i[0];

      if (computecolorfield_)
        colorfield_i[0] = (Wii / dens_i[0]) * mass_i[0] + sumj_Wij_dens_j_mass_j;
    }
  }
}

/*---------------------------------------------------------------------------*
 | evaluate continuity equation                               sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityBase::ContinuityEquation() const
{
  // get reference to neighbor pair data
  const SPHNeighborPairData& neighborpairdata = neighborpairs_->GetRefToNeighborPairData();

  // iterate over particle types
  for (auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // no density integration for boundary or rigid particles
    if (type_i == PARTICLEENGINE::BoundaryPhase or type_i == PARTICLEENGINE::RigidPhase) continue;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get material for current particle type
    const MAT::PAR::ParticleMaterialBase* material_i =
        particlematerial_->GetPtrToParticleMatParameter(type_i);

    // iterate over particles in container
    for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
    {
      // get reference to vector of neighbor pairs of current particle
      const std::vector<SPHNeighborPair>& currentNeighborPairs =
          (neighborpairdata[type_i])[particle_i];

      // check for neighbor pairs of current particle
      if (currentNeighborPairs.empty()) continue;

      // declare pointer variables for particle i
      const double *vel_i, *dens_i;
      double* densdot_i;

      // get pointer to particle states
      dens_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);
      densdot_i = container_i->GetPtrToParticleState(PARTICLEENGINE::DensityDot, particle_i);

      if (applytransportvelocity_)
        vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::ModifiedVelocity, particle_i);
      else
        vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);

      // initialize sum of evaluated kernel values for particle i due to neighbor particles j
      double sumj_densdot_ij(0.0);

      // iterate over neighbor pairs
      for (auto& neighborIt : currentNeighborPairs)
      {
        // access values of local index tuple of neighboring particle
        PARTICLEENGINE::TypeEnum type_j;
        PARTICLEENGINE::StatusEnum status_j;
        int particle_j;
        std::tie(type_j, status_j, particle_j) = neighborIt.first;

        // get reference to current particle pair
        const SPHParticlePair& particlepair = neighborIt.second;

        // get container of neighboring particles of current particle type and state
        PARTICLEENGINE::ParticleContainerShrdPtr container_j =
            particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

        // declare pointer variables for neighbor particle j
        const double *vel_j, *mass_j, *dens_j;

        // get pointer to particle states
        mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);

        if (type_j == PARTICLEENGINE::BoundaryPhase or type_j == PARTICLEENGINE::RigidPhase)
        {
          vel_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_j);

          dens_j = &(material_i->initDensity_);
        }
        else
        {
          if (applytransportvelocity_)
            vel_j =
                container_j->GetPtrToParticleState(PARTICLEENGINE::ModifiedVelocity, particle_j);
          else
            vel_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_j);

          dens_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j);
        }

        // relative velocity (use modified velocities in case of transport velocity formulation)
        double vel_ij[3];
        for (int i = 0; i < 3; ++i) vel_ij[i] = vel_i[i] - vel_j[i];

        // sum contribution of neighbor particle j
        sumj_densdot_ij += dens_i[0] * (mass_j[0] / dens_j[0]) * particlepair.dWdrij_ *
                           (particlepair.e_ij_[0] * vel_ij[0] + particlepair.e_ij_[1] * vel_ij[1] +
                               particlepair.e_ij_[2] * vel_ij[2]);
      }

      // add contributions of neighbor particles
      densdot_i[0] = sumj_densdot_ij;
    }
  }
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHDensitySummation::SPHDensitySummation(const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::SPHDensityBase(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | insert density evaluation dependent states                 sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensitySummation::InsertParticleStatesOfParticleTypes(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
    const
{
  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type = typeIt.first;

    // set of particle states for current particle type
    std::set<PARTICLEENGINE::StateEnum>& particlestates = typeIt.second;

    // no states for boundary or rigid particles
    if (type == PARTICLEENGINE::BoundaryPhase or type == PARTICLEENGINE::RigidPhase) continue;

    // states for density evaluation scheme
    particlestates.insert(PARTICLEENGINE::DensitySum);
  }
}

/*---------------------------------------------------------------------------*
 | compute density field                                      sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensitySummation::ComputeDensity() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHDensitySummation::ComputeDensity");

  // evaluate sum of weighted mass and colorfield
  SumWeightedMassAndColorfield();

  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // no density integration for boundary or rigid particles
    if (typeEnum == PARTICLEENGINE::BoundaryPhase or typeEnum == PARTICLEENGINE::RigidPhase)
      continue;

    // update density of all particles
    particlecontainerbundle_->UpdateStateSpecificContainer(
        0.0, PARTICLEENGINE::Density, 1.0, PARTICLEENGINE::DensitySum, typeEnum);
  }

  // refresh density of ghosted particles
  particleengineinterface_->RefreshParticlesOfSpecificStatesAndTypes(densitytorefresh_);
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHDensityIntegration::SPHDensityIntegration(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::SPHDensityBase(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | insert density evaluation dependent states                 sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityIntegration::InsertParticleStatesOfParticleTypes(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
    const
{
  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type = typeIt.first;

    // set of particle states for current particle type
    std::set<PARTICLEENGINE::StateEnum>& particlestates = typeIt.second;

    // no states for boundary or rigid particles
    if (type == PARTICLEENGINE::BoundaryPhase or type == PARTICLEENGINE::RigidPhase) continue;

    // states for density evaluation scheme
    particlestates.insert(PARTICLEENGINE::DensityDot);
  }
}

/*---------------------------------------------------------------------------*
 | compute density field                                      sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityIntegration::ComputeDensity() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHDensityIntegration::ComputeDensity");

  // evaluate continuity equation
  ContinuityEquation();

  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // no density integration for boundary or rigid particles
    if (typeEnum == PARTICLEENGINE::BoundaryPhase or typeEnum == PARTICLEENGINE::RigidPhase)
      continue;

    // update density of all particles
    particlecontainerbundle_->UpdateStateSpecificContainer(
        1.0, PARTICLEENGINE::Density, dt_, PARTICLEENGINE::DensityDot, typeEnum);
  }

  // refresh density of ghosted particles
  particleengineinterface_->RefreshParticlesOfSpecificStatesAndTypes(densitytorefresh_);
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHDensityPredictCorrect::SPHDensityPredictCorrect(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::SPHDensityBase(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | destructor                                                 sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHDensityPredictCorrect::~SPHDensityPredictCorrect()
{
  // note: destructor declaration here since at compile-time a complete type
  // of class T as used in class member std::unique_ptr<T> ptr_T_ is required
}

/*---------------------------------------------------------------------------*
 | init density handler                                       sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityPredictCorrect::Init()
{
  // call base class init
  SPHDensityBase::Init();

  // compute colorfield set to true
  computecolorfield_ = true;

  // init density correction handler
  InitDensityCorrectionHandler();
}

/*---------------------------------------------------------------------------*
 | setup density handler                                      sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityPredictCorrect::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
    const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
    const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle,
    const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs)
{
  // call base class setup
  SPHDensityBase::Setup(
      particleengineinterface, kernel, particlematerial, equationofstatebundle, neighborpairs);

  // setup density correction handler
  densitycorrection_->Setup();
}

/*---------------------------------------------------------------------------*
 | write restart of density handler                           sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityPredictCorrect::WriteRestart(
    const int step, const double time) const
{
  // call base class function
  SPHDensityBase::WriteRestart(step, time);

  // write restart of density correction handler
  densitycorrection_->WriteRestart(step, time);
}

/*---------------------------------------------------------------------------*
 | read restart of density handler                            sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityPredictCorrect::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // call base class function
  SPHDensityBase::ReadRestart(reader);

  // read restart of density correction handler
  densitycorrection_->ReadRestart(reader);
}

/*---------------------------------------------------------------------------*
 | insert density evaluation dependent states                 sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityPredictCorrect::InsertParticleStatesOfParticleTypes(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
    const
{
  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type = typeIt.first;

    // set of particle states for current particle type
    std::set<PARTICLEENGINE::StateEnum>& particlestates = typeIt.second;

    // no states for boundary or rigid particles
    if (type == PARTICLEENGINE::BoundaryPhase or type == PARTICLEENGINE::RigidPhase) continue;

    // states for density evaluation scheme
    particlestates.insert(
        {PARTICLEENGINE::DensityDot, PARTICLEENGINE::DensitySum, PARTICLEENGINE::Colorfield});
  }
}

/*---------------------------------------------------------------------------*
 | compute density field                                      sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityPredictCorrect::ComputeDensity() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHDensityPredictCorrect::ComputeDensity");

  // evaluate continuity equation
  ContinuityEquation();

  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // no density integration for boundary or rigid particles
    if (typeEnum == PARTICLEENGINE::BoundaryPhase or typeEnum == PARTICLEENGINE::RigidPhase)
      continue;

    // update density of all particles
    particlecontainerbundle_->UpdateStateSpecificContainer(
        1.0, PARTICLEENGINE::Density, dt_, PARTICLEENGINE::DensityDot, typeEnum);
  }

  // refresh density of ghosted particles
  particleengineinterface_->RefreshParticlesOfSpecificStatesAndTypes(densitytorefresh_);

  // evaluate sum of weighted mass and colorfield
  SumWeightedMassAndColorfield();

  // correct density of interior/surface particles
  CorrectDensity();

  // refresh density of ghosted particles
  particleengineinterface_->RefreshParticlesOfSpecificStatesAndTypes(densitytorefresh_);
}

/*---------------------------------------------------------------------------*
 | init density correction handler                            sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityPredictCorrect::InitDensityCorrectionHandler()
{
  // get type of density correction scheme
  INPAR::PARTICLE::DensityCorrectionScheme densitycorrectionscheme =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::DensityCorrectionScheme>(
          params_sph_, "DENSITYCORRECTION");

  // create density correction handler
  switch (densitycorrectionscheme)
  {
    case INPAR::PARTICLE::InteriorCorrection:
    {
      densitycorrection_ = std::unique_ptr<PARTICLEINTERACTION::SPHDensityCorrectionInterior>(
          new PARTICLEINTERACTION::SPHDensityCorrectionInterior());
      break;
    }
    case INPAR::PARTICLE::NormalizedCorrection:
    {
      densitycorrection_ = std::unique_ptr<PARTICLEINTERACTION::SPHDensityCorrectionNormalized>(
          new PARTICLEINTERACTION::SPHDensityCorrectionNormalized());
      break;
    }
    case INPAR::PARTICLE::RandlesCorrection:
    {
      densitycorrection_ = std::unique_ptr<PARTICLEINTERACTION::SPHDensityCorrectionRandles>(
          new PARTICLEINTERACTION::SPHDensityCorrectionRandles());
      break;
    }
    default:
    {
      dserror("no density correction scheme set via parameter 'DENSITYCORRECTION'!");
      break;
    }
  }

  // init density correction handler
  densitycorrection_->Init();
}

/*---------------------------------------------------------------------------*
 | correct density of interior/surface particles              sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityPredictCorrect::CorrectDensity() const
{
  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // no pressure computation for boundary or rigid particles
    if (typeEnum == PARTICLEENGINE::BoundaryPhase or typeEnum == PARTICLEENGINE::RigidPhase)
      continue;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // declare pointer variables for particle
    const double *denssum, *colorfield;
    double* dens;

    // get pointer to particle state
    dens = container->GetPtrToParticleState(PARTICLEENGINE::Density, 0);
    denssum = container->GetPtrToParticleState(PARTICLEENGINE::DensitySum, 0);
    colorfield = container->GetPtrToParticleState(PARTICLEENGINE::Colorfield, 0);

    // get material for current particle type
    const MAT::PAR::ParticleMaterialBase* material =
        particlematerial_->GetPtrToParticleMatParameter(typeEnum);

    // get equation of state for current particle type
    std::shared_ptr<SPHEquationOfStateBase> equationofstate =
        equationofstatebundle_->GetSpecificEquationOfState(typeEnum);

    // iterate over owned particles of current type
    for (int i = 0; i < particlestored; ++i)
    {
      if (colorfield[i] >= 1.0)
      {
        // set corrected density of interior particles
        densitycorrection_->CorrectedDensityInterior(&denssum[i], &dens[i]);
      }
      else
      {
        double dens_bc = 0.0;
        if (densitycorrection_->ComputeDensityBC())
        {
          double press_bc = 0.0;
          dens_bc = equationofstate->PressureToDensity(press_bc, material->initDensity_);
        }

        // set corrected density of free surface particles
        densitycorrection_->CorrectedDensityFreeSurface(
            &denssum[i], &colorfield[i], &dens_bc, &dens[i]);
      }
    }
  }
}
