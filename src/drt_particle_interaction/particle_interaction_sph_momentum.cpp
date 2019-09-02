/*---------------------------------------------------------------------------*/
/*! \file
\brief momentum handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_momentum.H"

#include "particle_interaction_sph_kernel.H"
#include "particle_interaction_material_handler.H"
#include "particle_interaction_sph_equationofstate.H"
#include "particle_interaction_sph_equationofstate_bundle.H"
#include "particle_interaction_sph_neighbor_pairs.H"
#include "particle_interaction_sph_momentum_formulation.H"
#include "particle_interaction_sph_artificialviscosity.H"

#include "particle_interaction_utils.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHMomentum::SPHMomentum(const Teuchos::ParameterList& params)
    : params_sph_(params),
      boundaryparticleinteraction_(
          DRT::INPUT::IntegralValue<INPAR::PARTICLE::BoundaryParticleInteraction>(
              params_sph_, "BOUNDARYPARTICLEINTERACTION")),
      transportvelocityformulation_(
          DRT::INPUT::IntegralValue<INPAR::PARTICLE::TransportVelocityFormulation>(
              params_sph_, "TRANSPORTVELOCITYFORMULATION")),
      applytransportvelocity_(false),
      norelativevelocitycontribution_(DRT::INPUT::IntegralValue<int>(params_sph_, "NO_RELVEL_TERM"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | destructor                                                 sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHMomentum::~SPHMomentum()
{
  // note: destructor declaration here since at compile-time a complete type
  // of class T as used in class member std::unique_ptr<T> ptr_T_ is required
}

/*---------------------------------------------------------------------------*
 | init momentum handler                                      sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentum::Init()
{
  // init momentum formulation handler
  InitMomentumFormulationHandler();

  // init artificial viscosity handler
  InitArtificialViscosityHandler();

  // check if transport velocity formulation is applied
  if (transportvelocityformulation_ !=
      INPAR::PARTICLE::TransportVelocityFormulation::NoTransportVelocity)
    applytransportvelocity_ = true;

  // safety checks
  if ((not applytransportvelocity_) and norelativevelocitycontribution_)
    dserror(
        "the parameter 'NO_RELVEL_TERM' only makes sense if transport velocity formulation "
        "'TRANSPORT_VELOCITY' is applied!");
}

/*---------------------------------------------------------------------------*
 | setup momentum handler                                     sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentum::Setup(
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

  // setup momentum formulation handler
  momentumformulation_->Setup();

  // setup artificial viscosity handler
  artificialviscosity_->Setup();

  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  // allocate memory to hold particle types
  fluidmaterial_.resize(typevectorsize);

  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // no fluid material for boundary or rigid particles
    if (type_i == PARTICLEENGINE::BoundaryPhase or type_i == PARTICLEENGINE::RigidPhase) continue;

    fluidmaterial_[type_i] = dynamic_cast<const MAT::PAR::ParticleMaterialSPHFluid*>(
        particlematerial_->GetPtrToParticleMatParameter(type_i));
  }
}

/*---------------------------------------------------------------------------*
 | write restart of momentum handler                          sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentum::WriteRestart(const int step, const double time) const
{
  // write restart of momentum formulation handler
  momentumformulation_->WriteRestart(step, time);

  // write restart of artificial viscosity handler
  artificialviscosity_->WriteRestart(step, time);
}

/*---------------------------------------------------------------------------*
 | read restart of momentum handler                           sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentum::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // read restart of momentum formulation handler
  momentumformulation_->ReadRestart(reader);

  // read restart of artificial viscosity handler
  artificialviscosity_->ReadRestart(reader);
}

/*---------------------------------------------------------------------------*
 | insert momentum evaluation dependent states                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentum::InsertParticleStatesOfParticleTypes(
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

    // additional states for transport velocity formulation
    if (applytransportvelocity_)
      particlestates.insert(
          {PARTICLEENGINE::ModifiedVelocity, PARTICLEENGINE::ModifiedAcceleration});
  }
}

/*---------------------------------------------------------------------------*
 | add momentum contribution to acceleration field            sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentum::AddAccelerationContribution() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHMomentum::AddAccelerationContribution");

  // momentum equation (particle contribution)
  MomentumEquationParticleContribution();
}

/*---------------------------------------------------------------------------*
 | init momentum formulation handler                          sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentum::InitMomentumFormulationHandler()
{
  // get type of smoothed particle hydrodynamics momentum formulation
  INPAR::PARTICLE::MomentumFormulationType momentumformulationtype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::MomentumFormulationType>(
          params_sph_, "MOMENTUMFORMULATION");

  // create momentum formulation handler
  switch (momentumformulationtype)
  {
    case INPAR::PARTICLE::AdamiMomentumFormulation:
    {
      momentumformulation_ = std::unique_ptr<PARTICLEINTERACTION::SPHMomentumFormulationAdami>(
          new PARTICLEINTERACTION::SPHMomentumFormulationAdami());
      break;
    }
    case INPAR::PARTICLE::MonaghanMomentumFormulation:
    {
      momentumformulation_ = std::unique_ptr<PARTICLEINTERACTION::SPHMomentumFormulationMonaghan>(
          new PARTICLEINTERACTION::SPHMomentumFormulationMonaghan());
      break;
    }
    default:
    {
      dserror("unknown acceleration formulation type!");
      break;
    }
  }

  // init momentum formulation handler
  momentumformulation_->Init();
}

/*---------------------------------------------------------------------------*
 | init artificial viscosity handler                          sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentum::InitArtificialViscosityHandler()
{
  // create artificial viscosity handler
  artificialviscosity_ = std::unique_ptr<PARTICLEINTERACTION::SPHArtificialViscosity>(
      new PARTICLEINTERACTION::SPHArtificialViscosity());

  // init artificial viscosity handler
  artificialviscosity_->Init();
}

/*---------------------------------------------------------------------------*
 | momentum equation (particle contribution)                  sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentum::MomentumEquationParticleContribution() const
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::SPHMomentum::MomentumEquationParticleContribution");

  // iterate over particle pairs
  for (auto& particlepair : neighborpairs_->GetRefToParticlePairData())
  {
    // access values of local index tuples of particle i and j
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = particlepair.tuple_i_;

    PARTICLEENGINE::TypeEnum type_j;
    PARTICLEENGINE::StatusEnum status_j;
    int particle_j;
    std::tie(type_j, status_j, particle_j) = particlepair.tuple_j_;

    // check for boundary or rigid particles
    bool isboundaryrigid_i =
        (type_i == PARTICLEENGINE::BoundaryPhase or type_i == PARTICLEENGINE::RigidPhase);
    bool isboundaryrigid_j =
        (type_j == PARTICLEENGINE::BoundaryPhase or type_j == PARTICLEENGINE::RigidPhase);

    // no momentum evaluation for both boundary or rigid particles
    if (isboundaryrigid_i and (isboundaryrigid_j or status_j == PARTICLEENGINE::Ghosted)) continue;

    // get corresponding particle containers
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    PARTICLEENGINE::ParticleContainer* container_j =
        particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

    // get material for particle types
    const MAT::PAR::ParticleMaterialSPHFluid* material_i =
        (isboundaryrigid_i) ? fluidmaterial_[type_j] : fluidmaterial_[type_i];
    const MAT::PAR::ParticleMaterialSPHFluid* material_j =
        (isboundaryrigid_j) ? fluidmaterial_[type_i] : fluidmaterial_[type_j];

    // get equation of state for particle types
    const PARTICLEINTERACTION::SPHEquationOfStateBase* equationofstate_i;
    if (not isboundaryrigid_i)
      equationofstate_i = equationofstatebundle_->GetPtrToSpecificEquationOfState(type_i);

    const PARTICLEINTERACTION::SPHEquationOfStateBase* equationofstate_j;
    if (not isboundaryrigid_j)
      equationofstate_j = equationofstatebundle_->GetPtrToSpecificEquationOfState(type_j);

    // declare pointer variables for particle i and j
    const double *vel_i, *rad_i, *mass_i, *dens_i, *press_i;
    const double* mod_vel_i = nullptr;
    double *acc_i = nullptr, *mod_acc_i = nullptr;

    const double *vel_j, *rad_j, *mass_j, *dens_j, *press_j;
    const double* mod_vel_j = nullptr;
    double *acc_j = nullptr, *mod_acc_j = nullptr;

    double temp_dens(0.0);

    // get pointer to particle states
    rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);

    if (isboundaryrigid_i)
    {
      mass_i = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);
      press_i = container_i->GetPtrToParticleState(PARTICLEENGINE::BoundaryPressure, particle_i);
      vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::BoundaryVelocity, particle_i);

      temp_dens = equationofstate_j->PressureToDensity(press_i[0], material_j->initDensity_);
      dens_i = &temp_dens;
    }
    else
    {
      mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
      dens_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);
      press_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Pressure, particle_i);
      vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);
      acc_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Acceleration, particle_i);

      if (applytransportvelocity_)
      {
        mod_vel_i =
            container_i->GetPtrToParticleState(PARTICLEENGINE::ModifiedVelocity, particle_i);
        mod_acc_i =
            container_i->GetPtrToParticleState(PARTICLEENGINE::ModifiedAcceleration, particle_i);
      }
    }

    // get pointer to particle states
    rad_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_j);

    if (isboundaryrigid_j)
    {
      mass_j = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
      press_j = container_j->GetPtrToParticleState(PARTICLEENGINE::BoundaryPressure, particle_j);
      vel_j = container_j->GetPtrToParticleState(PARTICLEENGINE::BoundaryVelocity, particle_j);

      temp_dens = equationofstate_i->PressureToDensity(press_j[0], material_i->initDensity_);
      dens_j = &temp_dens;
    }
    else
    {
      mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);
      dens_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j);
      press_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Pressure, particle_j);
      vel_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_j);

      if (status_j == PARTICLEENGINE::Owned)
        acc_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Acceleration, particle_j);

      if (applytransportvelocity_)
      {
        mod_vel_j =
            container_j->GetPtrToParticleState(PARTICLEENGINE::ModifiedVelocity, particle_j);
        mod_acc_j =
            container_j->GetPtrToParticleState(PARTICLEENGINE::ModifiedAcceleration, particle_j);
      }
    }

    // determine weather viscous contribution are evaluated
    bool evaluateviscouscontributions = true;
    if ((isboundaryrigid_i or isboundaryrigid_j) and
        boundaryparticleinteraction_ == INPAR::PARTICLE::FreeSlipBoundaryParticle)
      evaluateviscouscontributions = false;

    // evaluate specific coefficient
    double speccoeff_ij(0.0);
    double speccoeff_ji(0.0);
    momentumformulation_->SpecificCoefficient(dens_i, dens_j, mass_i, mass_j, particlepair.dWdrij_,
        particlepair.dWdrji_, &speccoeff_ij, &speccoeff_ji);

    // evaluate pressure gradient
    momentumformulation_->PressureGradient(dens_i, dens_j, press_i, press_j, speccoeff_ij,
        speccoeff_ji, particlepair.e_ij_, acc_i, acc_j);

    // evaluate shear forces
    if (evaluateviscouscontributions)
    {
      // get factor from kernel space dimension
      int kernelfac = 0;
      kernel_->KernelSpaceDimension(kernelfac);
      kernelfac += 2;

      // evaluate shear forces
      momentumformulation_->ShearForces(dens_i, dens_j, vel_i, vel_j, kernelfac,
          material_i->dynamicViscosity_, material_j->dynamicViscosity_, material_i->bulkViscosity_,
          material_j->bulkViscosity_, particlepair.absdist_, speccoeff_ij, speccoeff_ji,
          particlepair.e_ij_, acc_i, acc_j);
    }

    // apply transport velocity formulation
    if (applytransportvelocity_)
    {
      if (transportvelocityformulation_ ==
          INPAR::PARTICLE::TransportVelocityFormulation::StandardTransportVelocity)
      {
        // evaluate background pressure (standard formulation)
        momentumformulation_->StandardBackgroundPressure(dens_i, dens_j,
            material_i->backgroundPressure_, material_j->backgroundPressure_, speccoeff_ij,
            speccoeff_ji, particlepair.e_ij_, mod_acc_i, mod_acc_j);
      }
      else if (transportvelocityformulation_ ==
               INPAR::PARTICLE::TransportVelocityFormulation::GeneralizedTransportVelocity)
      {
        // modified support radius and first derivative of kernel
        const double mod_rad_i = (mod_acc_i) ? kernel_->SmoothingLength(rad_i[0]) : 0.0;
        const double mod_dWdrij =
            (mod_acc_i) ? kernel_->dWdrij(particlepair.absdist_, mod_rad_i) : 0.0;

        double mod_rad_j = 0.0;
        double mod_dWdrji = 0.0;
        if (mod_acc_j and mod_acc_i and rad_i[0] == rad_j[0])
        {
          mod_rad_j = mod_rad_i;
          mod_dWdrji = mod_dWdrij;
        }
        else if (mod_acc_j)
        {
          mod_rad_j = kernel_->SmoothingLength(rad_j[0]);
          mod_dWdrji = kernel_->dWdrij(particlepair.absdist_, mod_rad_j);
        }

        // modified background pressure
        const double mod_bg_press_i =
            (mod_acc_i) ? std::min(std::abs(10.0 * press_i[0]), material_i->backgroundPressure_)
                        : 0.0;
        const double mod_bg_press_j =
            (mod_acc_j) ? std::min(std::abs(10.0 * press_j[0]), material_j->backgroundPressure_)
                        : 0.0;

        // evaluate background pressure (generalized formulation)
        momentumformulation_->GeneralizedBackgroundPressure(dens_i, dens_j, mass_i, mass_j,
            mod_bg_press_i, mod_bg_press_j, mod_dWdrij, mod_dWdrji, particlepair.e_ij_, mod_acc_i,
            mod_acc_j);
      }

      // evaluate convection of momentum with relative velocity
      if (not norelativevelocitycontribution_)
      {
        // evaluate modified velocity contribution
        momentumformulation_->ModifiedVelocityContribution(dens_i, dens_j, vel_i, vel_j, mod_vel_i,
            mod_vel_j, speccoeff_ij, speccoeff_ji, particlepair.e_ij_, acc_i, acc_j);
      }
    }

    // evaluate artificial viscosity
    if (evaluateviscouscontributions and
        (material_i->artificialViscosity_ > 0.0 or material_j->artificialViscosity_ > 0.0))
    {
      // get smoothing length
      const double h_i = kernel_->SmoothingLength(rad_i[0]);
      const double h_j = (rad_i[0] == rad_j[0]) ? h_i : kernel_->SmoothingLength(rad_j[0]);

      // particle averaged smoothing length
      const double h_ij = 0.5 * (h_i + h_j);

      // get speed of sound
      const double c_i = material_i->SpeedOfSound();
      const double c_j = (type_i == type_j) ? c_i : material_j->SpeedOfSound();

      // particle averaged speed of sound
      const double c_ij = 0.5 * (c_i + c_j);

      // particle averaged density
      const double dens_ij = 0.5 * (dens_i[0] + dens_j[0]);

      // evaluate artificial viscosity
      artificialviscosity_->ArtificialViscosity(vel_i, vel_j, mass_i, mass_j,
          material_i->artificialViscosity_, material_j->artificialViscosity_, particlepair.dWdrij_,
          particlepair.dWdrji_, dens_ij, h_ij, c_ij, particlepair.absdist_, particlepair.e_ij_,
          acc_i, acc_j);
    }
  }
}
