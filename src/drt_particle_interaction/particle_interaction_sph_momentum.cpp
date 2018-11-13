/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_momentum.cpp

\brief momentum handler for smoothed particle hydrodynamics (SPH) interactions

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
#include "particle_interaction_sph_momentum.H"

#include "particle_interaction_sph_kernel.H"
#include "particle_interaction_material_handler.H"
#include "particle_interaction_sph_equationofstate.H"
#include "particle_interaction_sph_equationofstate_bundle.H"
#include "particle_interaction_sph_neighbor_pairs.H"
#include "particle_interaction_sph_momentum_formulation.H"
#include "particle_interaction_sph_artificialviscosity.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_dserror.H"

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
      norelativevelocitycontribution_(
          DRT::INPUT::IntegralValue<int>(params_sph_, "NO_RELVEL_TERM")),
      dampingfactor_(params_sph_.get<double>("VISCOUS_DAMPING"))
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
  // declare temp variables
  double temp(0.0);

  // iterate over particle types
  for (auto& typeIt : neighborpairs_->GetRefToNeighborPairsMap())
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type_i = typeIt.first;

    // no acceleration evaluation for boundary or rigid particles
    if (type_i == PARTICLEENGINE::BoundaryPhase or type_i == PARTICLEENGINE::RigidPhase) continue;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get material for current particle type
    const MAT::PAR::ParticleMaterialSPHFluid* material_i =
        dynamic_cast<const MAT::PAR::ParticleMaterialSPHFluid*>(
            particlematerial_->GetPtrToParticleMatParameter(type_i));

    // get equation of state for current particle type
    std::shared_ptr<SPHEquationOfStateBase> equationofstate_i =
        equationofstatebundle_->GetSpecificEquationOfState(type_i);

    // particles of current type with neighbors
    const auto& currparticles = typeIt.second;

    // iterate over particles of current type
    for (auto& particleIt : currparticles)
    {
      // get local index of particle i
      const int particle_i = particleIt.first;

      // declare pointer variables for particle i
      const double *vel_i, *rad_i, *mass_i, *dens_i, *press_i, *mod_vel_i;
      double *acc_i, *mod_acc_i;

      // get pointer to particle states
      vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);
      acc_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Acceleration, particle_i);
      rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
      mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
      dens_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);
      press_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Pressure, particle_i);

      if (applytransportvelocity_)
      {
        mod_vel_i =
            container_i->GetPtrToParticleState(PARTICLEENGINE::ModifiedVelocity, particle_i);
        mod_acc_i =
            container_i->GetPtrToParticleState(PARTICLEENGINE::ModifiedAcceleration, particle_i);
      }

      // add contributions due to artificial viscous damping
      if (dampingfactor_ > 0.0)
        for (int i = 0; i < 3; ++i) acc_i[i] -= dampingfactor_ * vel_i[i];

      // iterate over particle types of neighboring particles
      for (auto& neighborTypeIt : particleIt.second)
      {
        // get type of neighboring particles
        PARTICLEENGINE::TypeEnum type_j = neighborTypeIt.first;

        // get material for current particle type
        const MAT::PAR::ParticleMaterialSPHFluid* material_j = NULL;
        if (type_i == type_j or type_j == PARTICLEENGINE::BoundaryPhase or
            type_j == PARTICLEENGINE::RigidPhase)
          material_j = material_i;
        else
          material_j = dynamic_cast<const MAT::PAR::ParticleMaterialSPHFluid*>(
              particlematerial_->GetPtrToParticleMatParameter(type_j));

        // determine weather viscous contribution are evaluated
        bool evaluateviscouscontributions = true;
        if ((type_j == PARTICLEENGINE::BoundaryPhase or type_j == PARTICLEENGINE::RigidPhase) and
            boundaryparticleinteraction_ == INPAR::PARTICLE::FreeSlipBoundaryParticle)
          evaluateviscouscontributions = false;

        // iterate over particle status of neighboring particles
        for (auto& neighborStatusIt : neighborTypeIt.second)
        {
          // get status of neighboring particles of current type
          PARTICLEENGINE::StatusEnum status_j = neighborStatusIt.first;

          // get container of neighboring particles of current particle type and state
          PARTICLEENGINE::ParticleContainerShrdPtr container_j =
              particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

          // iterate over neighboring particles of current type and status
          for (auto& neighborParticleIt : neighborStatusIt.second)
          {
            // get local index of neighbor particle j
            const int particle_j = neighborParticleIt.first;

            // get reference to particle pair
            const ParticlePairSPH& particlepair = neighborParticleIt.second;

            // declare pointer variables for neighbor particle j
            const double *vel_j, *rad_j, *mass_j, *dens_j, *press_j, *mod_vel_j;

            // get pointer to particle states
            if (type_j == PARTICLEENGINE::BoundaryPhase or type_j == PARTICLEENGINE::RigidPhase)
            {
              mass_j = mass_i;
              dens_j = &temp;
              press_j =
                  container_j->GetPtrToParticleState(PARTICLEENGINE::BoundaryPressure, particle_j);
              vel_j =
                  container_j->GetPtrToParticleState(PARTICLEENGINE::BoundaryVelocity, particle_j);

              temp = equationofstate_i->PressureToDensity(press_j[0], material_i->initDensity_);
            }
            else
            {
              mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);
              dens_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j);
              press_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Pressure, particle_j);
              vel_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_j);

              if (applytransportvelocity_)
                mod_vel_j = container_j->GetPtrToParticleState(
                    PARTICLEENGINE::ModifiedVelocity, particle_j);
            }

            // evaluate specific coefficient
            double specificcoefficient_ij(0.0);
            momentumformulation_->SpecificCoefficient(
                dens_i, dens_j, mass_i, mass_j, particlepair.dWdrij_, specificcoefficient_ij);

            // evaluate pressure gradient
            momentumformulation_->PressureGradient(dens_i, dens_j, press_i, press_j,
                specificcoefficient_ij, particlepair.e_ij_, acc_i);

            // evaluate shear forces
            if (evaluateviscouscontributions)
            {
              // get dynamic shear viscosity
              const double visc_i = material_i->dynamicViscosity_;
              const double visc_j = material_j->dynamicViscosity_;

              // get bulk viscosity
              const double bulk_visc_i = material_i->bulkViscosity_;
              const double bulk_visc_j = material_j->bulkViscosity_;

              // get factor from kernel space dimension
              int kernelfac = 0;
              kernel_->KernelSpaceDimension(kernelfac);
              kernelfac += 2;

              momentumformulation_->ShearForces(dens_i, dens_j, vel_i, vel_j, kernelfac, visc_i,
                  visc_j, bulk_visc_i, bulk_visc_j, particlepair.absdist_, specificcoefficient_ij,
                  particlepair.e_ij_, acc_i);
            }

            // apply transport velocity formulation
            if (applytransportvelocity_)
            {
              // get background pressure
              const double backgroundpressure = material_i->backgroundPressure_;

              if (transportvelocityformulation_ ==
                  INPAR::PARTICLE::TransportVelocityFormulation::StandardTransportVelocity)
              {
                // evaluate background pressure (standard formulation)
                momentumformulation_->StandardBackgroundPressure(dens_i, dens_j, backgroundpressure,
                    specificcoefficient_ij, particlepair.e_ij_, mod_acc_i);
              }
              else if (transportvelocityformulation_ ==
                       INPAR::PARTICLE::TransportVelocityFormulation::GeneralizedTransportVelocity)
              {
                // modified support radius for generalized transport velocity formulation
                const double support_tilde = kernel_->SmoothingLength(rad_i[0]);

                // evaluate first derivative of kernel
                const double dWdrij_tilde = kernel_->dWdrij(particlepair.absdist_, support_tilde);

                // modified background pressure
                const double backgroundpressure_tilde =
                    std::min(std::abs(10.0 * press_i[0]), backgroundpressure);

                // evaluate background pressure (generalized formulation)
                momentumformulation_->GeneralizedBackgroundPressure(dens_i, mass_i, mass_j,
                    backgroundpressure_tilde, dWdrij_tilde, particlepair.e_ij_, mod_acc_i);
              }

              // evaluate convection of momentum with relative velocity
              if (not norelativevelocitycontribution_)
              {
                // evaluate modified velocity contribution
                if (type_j == PARTICLEENGINE::BoundaryPhase or type_j == PARTICLEENGINE::RigidPhase)
                  momentumformulation_->ModifiedVelocityContribution(dens_i, dens_j, vel_i, nullptr,
                      mod_vel_i, nullptr, specificcoefficient_ij, particlepair.e_ij_, acc_i);
                else
                  momentumformulation_->ModifiedVelocityContribution(dens_i, dens_j, vel_i, vel_j,
                      mod_vel_i, mod_vel_j, specificcoefficient_ij, particlepair.e_ij_, acc_i);
              }
            }

            // artificial viscosity
            if (evaluateviscouscontributions and material_i->artificialViscosity_ > 0.0)
            {
              // get pointer to particle states
              rad_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_j);

              // get smoothing length of particle types
              const double h_i = kernel_->SmoothingLength(rad_i[0]);
              const double h_j = kernel_->SmoothingLength(rad_j[0]);

              // get speed of sound of particle types
              const double c_i = material_i->SpeedOfSound();
              const double c_j = material_j->SpeedOfSound();

              // evaluate artificial viscosity
              artificialviscosity_->ArtificialViscosity(dens_i, dens_j, vel_i, vel_j, mass_j,
                  material_i->artificialViscosity_, particlepair.dWdrij_, h_i, h_j, c_i, c_j,
                  particlepair.absdist_, particlepair.e_ij_, acc_i);
            }
          }
        }
      }
    }
  }
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
