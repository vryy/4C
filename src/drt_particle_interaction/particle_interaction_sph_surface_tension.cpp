/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_surface_tension.cpp

\brief surface tension handler for smoothed particle hydrodynamics (SPH) interactions

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
#include "particle_interaction_sph_surface_tension.H"

#include "particle_interaction_sph_kernel.H"
#include "particle_interaction_material_handler.H"
#include "particle_interaction_sph_equationofstate.H"
#include "particle_interaction_sph_equationofstate_bundle.H"
#include "particle_interaction_sph_neighbor_pairs.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_dserror.H"

#include "../drt_lib/drt_globalproblem.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHSurfaceTensionBase::SPHSurfaceTensionBase(
    const Teuchos::ParameterList& params)
    : params_sph_(params),
      time_(0.0),
      surfacetensionrampfctnumber_(params.get<int>("SURFACETENSION_RAMP_FUNCT")),
      haveboundaryorrigidparticles_(false)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init surface tension handler                               sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHSurfaceTensionBase::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup surface tension handler                              sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHSurfaceTensionBase::Setup(
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
}

/*---------------------------------------------------------------------------*
 | write restart of surface tension handler                   sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHSurfaceTensionBase::WriteRestart(
    const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of surface tension handler                    sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHSurfaceTensionBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | set current time                                           sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHSurfaceTensionBase::SetCurrentTime(const double currenttime)
{
  time_ = currenttime;
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::SPHSurfaceTensionContinuumSurfaceForce(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::SPHSurfaceTensionBase(params),
      surfacetensioncoefficient_(params_sph_.get<double>("SURFACETENSIONCOEFFICIENT")),
      staticcontactangle_(params_sph_.get<double>("STATICCONTACTANGLE"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init surface tension handler                               sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::Init()
{
  // call base class init
  SPHSurfaceTensionBase::Init();

  // safety check
  if (not(surfacetensioncoefficient_ > 0.0))
    dserror(
        "the parameter 'SURFACETENSIONCOEFFICIENT' must be positiv when applying continuum surface "
        "force formulation!");
}

/*---------------------------------------------------------------------------*
 | insert surface tension evaluation dependent states         sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::
    InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes)
{
  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type = typeIt.first;

    // have boundary or rigid particles
    if (type == PARTICLEENGINE::BoundaryPhase or type == PARTICLEENGINE::RigidPhase)
      haveboundaryorrigidparticles_ = true;
  }

  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type = typeIt.first;

    // set of particle states for current particle type
    std::set<PARTICLEENGINE::StateEnum>& particlestates = typeIt.second;

    // no states for boundary or rigid particles
    if (type == PARTICLEENGINE::BoundaryPhase or type == PARTICLEENGINE::RigidPhase) continue;

    // states for surface tension evaluation scheme
    particlestates.insert({PARTICLEENGINE::ColorfieldGradient, PARTICLEENGINE::InterfaceNormal,
        PARTICLEENGINE::Curvature});

    if (haveboundaryorrigidparticles_)
      particlestates.insert({PARTICLEENGINE::UnitWallNormal, PARTICLEENGINE::WallDistance});
  }
}

/*---------------------------------------------------------------------------*
 | add surface tension contribution to acceleration field     sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::AddAccelerationContribution()
    const
{
  // compute colorfield gradient
  ComputeColorfieldGradient();

  if (haveboundaryorrigidparticles_)
  {
    // compute wall normal and distance
    ComputeWallNormalAndDistance();

    // refresh colorfield gradient and wall distance
    RefreshColorfieldGradientAndWallDistance();

    // extrapolate colorfield gradient at triple point
    ExtrapolateColorfieldGradientAtTriplePoint();
  }

  // compute interface normal
  ComputeInterfaceNormal();

  if (haveboundaryorrigidparticles_)
  {
    // correct normal vector of particles close to triple point
    CorrectTriplePointNormal();
  }

  // refresh colorfield gradient and interface normal
  RefreshColorfieldGradientAndInterfaceNormal();

  // compute curvature and add acceleration contribution
  ComputeCurvatureAndAddAccelerationContribution();
}

/*---------------------------------------------------------------------------*
 | compute colorfield gradient                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::ComputeColorfieldGradient() const
{
  // iterate over particle types
  for (auto& typeIt : neighborpairs_->GetRefToNeighborPairsMap())
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type_i = typeIt.first;

    // no colorfield gradient evaluation for boundary or rigid particles
    if (type_i == PARTICLEENGINE::BoundaryPhase or type_i == PARTICLEENGINE::RigidPhase) continue;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear colorfield gradient state
    container_i->ClearState(PARTICLEENGINE::ColorfieldGradient);

    // particles of current type with neighbors
    const auto& currparticles = typeIt.second;

    // iterate over particles of current type
    for (auto& particleIt : currparticles)
    {
      // get local index of particle i
      const int particle_i = particleIt.first;

      // declare pointer variables for particle i
      const double *mass_i, *dens_i;
      double* colorfieldgrad_i;

      // get pointer to particle states
      mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
      dens_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);
      colorfieldgrad_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::ColorfieldGradient, particle_i);

      // volume of particle i
      const double V_i = mass_i[0] / dens_i[0];

      // iterate over particle types of neighboring particles
      for (auto& neighborTypeIt : particleIt.second)
      {
        // get type of neighboring particles
        PARTICLEENGINE::TypeEnum type_j = neighborTypeIt.first;

        // no evaluation for particles of same type or neighboring boundary or rigid particles
        if (type_i == type_j or type_j == PARTICLEENGINE::BoundaryPhase or
            type_j == PARTICLEENGINE::RigidPhase)
          continue;

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
            const double *mass_j, *dens_j;

            // get pointer to particle states
            mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);
            dens_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j);

            // volume of particle j
            const double V_j = mass_j[0] / dens_j[0];

            const double fac = (V_i * V_i + V_j * V_j) * (dens_i[0] / (dens_i[0] + dens_j[0])) *
                               particlepair.dWdrij_;

            // sum contribution of neighbor particle j
            for (int i = 0; i < 3; ++i) colorfieldgrad_i[i] += fac * particlepair.e_ij_[i];
          }
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*
 | compute wall normal and distance                           sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::ComputeWallNormalAndDistance()
    const
{
  // iterate over particle types
  for (auto& typeIt : neighborpairs_->GetRefToNeighborPairsMap())
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type_i = typeIt.first;

    // no colorfield gradient evaluation for boundary or rigid particles
    if (type_i == PARTICLEENGINE::BoundaryPhase or type_i == PARTICLEENGINE::RigidPhase) continue;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get material for current particle type
    const MAT::PAR::ParticleMaterialBase* material_i =
        particlematerial_->GetPtrToParticleMatParameter(type_i);

    // clear unit wall normal state
    container_i->ClearState(PARTICLEENGINE::UnitWallNormal);

    // particles of current type with neighbors
    const auto& currparticles = typeIt.second;

    // iterate over particles of current type
    for (auto& particleIt : currparticles)
    {
      // get local index of particle i
      const int particle_i = particleIt.first;

      // declare pointer variables for particle i
      const double *rad_i, *mass_i, *dens_i;
      double *wallnormal_i, *walldistance_i;

      // get pointer to particle states
      rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
      mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
      dens_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);
      wallnormal_i = container_i->GetPtrToParticleState(PARTICLEENGINE::UnitWallNormal, particle_i);
      walldistance_i = container_i->GetPtrToParticleState(PARTICLEENGINE::WallDistance, particle_i);

      // volume of particle i (current volume)
      const double V_i = mass_i[0] / dens_i[0];

      // volume of boundary particle j (initial volume)
      const double V_j = mass_i[0] / material_i->initDensity_;

      // set support radius of particle i as initial wall distance
      walldistance_i[0] = rad_i[0];

      // iterate over particle types of neighboring particles
      for (auto& neighborTypeIt : particleIt.second)
      {
        // get type of neighboring particles
        PARTICLEENGINE::TypeEnum type_j = neighborTypeIt.first;

        // no evaluation for neighboring non-boundary and non-rigid particles
        if (type_j != PARTICLEENGINE::BoundaryPhase and type_j != PARTICLEENGINE::RigidPhase)
          continue;

        // iterate over particle status of neighboring particles
        for (auto& neighborStatusIt : neighborTypeIt.second)
        {
          // iterate over neighboring particles of current type and status
          for (auto& neighborParticleIt : neighborStatusIt.second)
          {
            // get reference to particle pair
            const ParticlePairSPH& particlepair = neighborParticleIt.second;

            // sum contribution of neighbor boundary particle j
            const double fac = -((V_j * V_j) / V_i) * particlepair.dWdrij_;
            for (int i = 0; i < 3; ++i) wallnormal_i[i] += fac * particlepair.e_ij_[i];
          }
        }
      }

      // norm of wall normal
      const double wallnormal_i_norm =
          std::sqrt(wallnormal_i[0] * wallnormal_i[0] + wallnormal_i[1] * wallnormal_i[1] +
                    wallnormal_i[2] * wallnormal_i[2]);

      // no interacting boundary particle
      if (not(wallnormal_i_norm > 0.0)) continue;

      // scale unit wall normal
      for (int i = 0; i < 3; ++i) wallnormal_i[i] = wallnormal_i[i] / wallnormal_i_norm;

      // iterate over particle types of neighboring particles
      for (auto& neighborTypeIt : particleIt.second)
      {
        // get type of neighboring particles
        PARTICLEENGINE::TypeEnum type_j = neighborTypeIt.first;

        // no evaluation for neighboring non-boundary and non-rigid particles
        if (type_j != PARTICLEENGINE::BoundaryPhase and type_j != PARTICLEENGINE::RigidPhase)
          continue;

        // iterate over particle status of neighboring particles
        for (auto& neighborStatusIt : neighborTypeIt.second)
        {
          // iterate over neighboring particles of current type and status
          for (auto& neighborParticleIt : neighborStatusIt.second)
          {
            // get reference to particle pair
            const ParticlePairSPH& particlepair = neighborParticleIt.second;

            // distance of particle i to neighboring boundary particle j
            double currentwalldistance =
                particlepair.absdist_ *
                (wallnormal_i[0] * particlepair.e_ij_[0] + wallnormal_i[1] * particlepair.e_ij_[1] +
                    wallnormal_i[2] * particlepair.e_ij_[2]);

            // update wall distance of particle i
            if (currentwalldistance < walldistance_i[0]) walldistance_i[0] = currentwalldistance;
          }
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*
 | refresh colorfield gradient and wall distance              sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::
    RefreshColorfieldGradientAndWallDistance() const
{
  // init map
  std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>> particlestatestotypes;

  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // no refreshing of surface tension states for boundary or rigid particles
    if (typeEnum == PARTICLEENGINE::BoundaryPhase or typeEnum == PARTICLEENGINE::RigidPhase)
      continue;

    // set state enums to map
    particlestatestotypes[typeEnum] = {
        PARTICLEENGINE::ColorfieldGradient, PARTICLEENGINE::WallDistance};
  }

  // refresh specific states of particles of specific types
  particleengineinterface_->RefreshSpecificStatesOfParticlesOfSpecificTypes(particlestatestotypes);
}

/*---------------------------------------------------------------------------*
 | extrapolate colorfield gradient at triple point            sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::
    ExtrapolateColorfieldGradientAtTriplePoint() const
{
  // map to store modified colorfield gradients of particles at triple point of all particle types
  std::map<PARTICLEENGINE::TypeEnum, std::map<int, std::vector<double>>>
      modcolorfieldgradofalltypes;

  // get kernel space dimension
  int kernelspacedim = 0;
  kernel_->KernelSpaceDimension(kernelspacedim);

  // iterate over particle types
  for (auto& typeIt : neighborpairs_->GetRefToNeighborPairsMap())
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type_i = typeIt.first;

    // no colorfield gradient extrapolation for boundary or rigid particles
    if (type_i == PARTICLEENGINE::BoundaryPhase or type_i == PARTICLEENGINE::RigidPhase) continue;

    // get reference to sub-map of modified colorfield gradients of current types
    auto& modcolorfieldgradofcurrtype = modcolorfieldgradofalltypes[type_i];

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get material for current particle type
    const MAT::PAR::ParticleMaterialBase* material_i =
        particlematerial_->GetPtrToParticleMatParameter(type_i);

    // particles of current type with neighbors
    const auto& currparticles = typeIt.second;

    // iterate over particles of current type
    for (auto& particleIt : currparticles)
    {
      // get local index of particle i
      const int particle_i = particleIt.first;

      // declare pointer variables for particle i
      const double *rad_i, *mass_i, *colorfieldgrad_i, *walldistance_i;

      // get pointer to particle states
      rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
      mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
      colorfieldgrad_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::ColorfieldGradient, particle_i);
      walldistance_i = container_i->GetPtrToParticleState(PARTICLEENGINE::WallDistance, particle_i);

      // evaluation only for particles close to boundary
      if (not(walldistance_i[0] < rad_i[0])) continue;

      // evaluation only for non-zero colorfield gradient
      const double colorfieldgrad_i_norm = std::sqrt(colorfieldgrad_i[0] * colorfieldgrad_i[0] +
                                                     colorfieldgrad_i[1] * colorfieldgrad_i[1] +
                                                     colorfieldgrad_i[2] * colorfieldgrad_i[2]);
      if (not(colorfieldgrad_i_norm > 0.0)) continue;

      // initial particle spacing
      const double initspacing_i =
          std::pow((mass_i[0] / material_i->initDensity_), (1.0 / kernelspacedim));

      // corrected wall distance and maximum correction distance
      const double dw_i = walldistance_i[0] - initspacing_i;
      const double dmax_i = kernel_->SmoothingLength(rad_i[0]);

      // determine correction factor
      double f_i = 1.0;
      if (dw_i < 0.0)
        f_i = 0.0;
      else if (dw_i < dmax_i)
        f_i = dw_i / dmax_i;

      // no extrapolation for current particle i
      if (f_i == 1.0) continue;

      // initialize sum of evaluated kernel values for particle i due to neighbor particles j
      double sumj_fj_Vj_Wij(0.0);
      double sumj_fj_Vj_Wij_CFGj[3];
      for (int i = 0; i < 3; ++i) sumj_fj_Vj_Wij_CFGj[i] = 0.0;

      // iterate over particle types of neighboring particles
      for (auto& neighborTypeIt : particleIt.second)
      {
        // get type of neighboring particles
        PARTICLEENGINE::TypeEnum type_j = neighborTypeIt.first;

        // no evaluation for neighboring boundary or rigid particles
        if (type_j == PARTICLEENGINE::BoundaryPhase or type_j == PARTICLEENGINE::RigidPhase)
          continue;

        // get material for current particle type
        const MAT::PAR::ParticleMaterialBase* material_j = NULL;
        if (type_i == type_j)
          material_j = material_i;
        else
          material_j = particlematerial_->GetPtrToParticleMatParameter(type_j);

        // change sign of colorfield gradient for different particle types
        double signfac = (type_i == type_j) ? 1.0 : -1.0;

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
            const double *rad_j, *mass_j, *dens_j, *colorfieldgrad_j, *walldistance_j;

            // get pointer to particle states
            rad_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_j);
            mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);
            dens_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j);
            colorfieldgrad_j =
                container_j->GetPtrToParticleState(PARTICLEENGINE::ColorfieldGradient, particle_j);
            walldistance_j =
                container_j->GetPtrToParticleState(PARTICLEENGINE::WallDistance, particle_j);

            // initial particle spacing
            const double initspacing_j =
                std::pow((mass_j[0] / material_j->initDensity_), (1.0 / kernelspacedim));

            // corrected wall distance and maximum correction distance
            const double dw_j = walldistance_j[0] - initspacing_j;
            const double dmax_j = kernel_->SmoothingLength(rad_j[0]);

            // determine correction factor
            double f_j = 1.0;
            if (dw_j < 0.0)
              f_j = 0.0;
            else if (dw_j < dmax_j)
              f_j = dw_j / dmax_j;

            // no contribution of current particle j
            if (f_j == 0.0) continue;

            // volume of particle j
            const double V_j = mass_j[0] / dens_j[0];

            const double fac = f_j * V_j * particlepair.Wij_;

            // initial estimate
            for (int i = 0; i < 3; ++i)
              sumj_fj_Vj_Wij_CFGj[i] += signfac * fac * colorfieldgrad_j[i];

            // correction factor
            sumj_fj_Vj_Wij += fac;
          }
        }
      }

      // evaluation only for particles with contributions from neighboring particles
      if (not(sumj_fj_Vj_Wij > 0.0)) continue;

      // get reference to modified colorfield gradient of particle i
      std::vector<double>& modcolorfieldgrad_i = modcolorfieldgradofcurrtype[particle_i];
      modcolorfieldgrad_i.resize(3);

      // determine modified colorfield gradient
      for (int i = 0; i < 3; ++i)
        modcolorfieldgrad_i[i] =
            f_i * colorfieldgrad_i[i] + (1.0 - f_i) * (sumj_fj_Vj_Wij_CFGj[i] / sumj_fj_Vj_Wij);

      // normalize modified colorfield gradient and scale with magnitude of original colorfield
      // gradient
      const double modcolorfieldgrad_i_norm =
          std::sqrt(modcolorfieldgrad_i[0] * modcolorfieldgrad_i[0] +
                    modcolorfieldgrad_i[1] * modcolorfieldgrad_i[1] +
                    modcolorfieldgrad_i[2] * modcolorfieldgrad_i[2]);
      for (int i = 0; i < 3; ++i)
        modcolorfieldgrad_i[i] =
            (colorfieldgrad_i_norm / modcolorfieldgrad_i_norm) * modcolorfieldgrad_i[i];
    }
  }

  // iterate over particle types
  for (auto& typeIt : modcolorfieldgradofalltypes)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type_i = typeIt.first;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // particles of current type with modified colorfield gradient
    const auto& currparticles = typeIt.second;

    // iterate over particles of current type
    for (auto& particleIt : currparticles)
    {
      // get local index of particle i
      const int particle_i = particleIt.first;

      // declare pointer variables for particle i
      double* colorfieldgrad_i;

      // get pointer to particle states
      colorfieldgrad_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::ColorfieldGradient, particle_i);

      // set modified colorfield gradient
      for (int i = 0; i < 3; ++i) colorfieldgrad_i[i] = particleIt.second[i];
    }
  }
}

/*---------------------------------------------------------------------------*
 | compute interface normal                                   sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::ComputeInterfaceNormal() const
{
  // iterate over particle types
  for (auto& typeIt : neighborpairs_->GetRefToNeighborPairsMap())
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type_i = typeIt.first;

    // no colorfield gradient evaluation for boundary or rigid particles
    if (type_i == PARTICLEENGINE::BoundaryPhase or type_i == PARTICLEENGINE::RigidPhase) continue;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear interface normal state
    container_i->ClearState(PARTICLEENGINE::InterfaceNormal);

    // particles of current type with neighbors
    const auto& currparticles = typeIt.second;

    // iterate over particles of current type
    for (auto& particleIt : currparticles)
    {
      // get local index of particle i
      const int particle_i = particleIt.first;

      // declare pointer variables for particle i
      const double *rad_i, *colorfieldgrad_i;
      double* interfacenormal_i;

      // get pointer to particle states
      rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
      colorfieldgrad_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::ColorfieldGradient, particle_i);
      interfacenormal_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::InterfaceNormal, particle_i);

      // norm of colorfield gradient
      const double colorfieldgrad_i_norm = std::sqrt(colorfieldgrad_i[0] * colorfieldgrad_i[0] +
                                                     colorfieldgrad_i[1] * colorfieldgrad_i[1] +
                                                     colorfieldgrad_i[2] * colorfieldgrad_i[2]);

      // scale colorfield gradient
      if (colorfieldgrad_i_norm > (1.0e-10 * rad_i[0]))
        for (int i = 0; i < 3; ++i)
          interfacenormal_i[i] = colorfieldgrad_i[i] / colorfieldgrad_i_norm;
    }
  }
}

/*---------------------------------------------------------------------------*
 | correct normal vector of particles close to triple point   sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::CorrectTriplePointNormal() const
{
  // get kernel space dimension
  int kernelspacedim = 0;
  kernel_->KernelSpaceDimension(kernelspacedim);

  // iterate over particle types
  for (auto& typeIt : neighborpairs_->GetRefToNeighborPairsMap())
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type_i = typeIt.first;

    // no curvature evaluation for boundary or rigid particles
    if (type_i == PARTICLEENGINE::BoundaryPhase or type_i == PARTICLEENGINE::RigidPhase) continue;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get material for current particle type
    const MAT::PAR::ParticleMaterialBase* material_i =
        particlematerial_->GetPtrToParticleMatParameter(type_i);

    // particles of current type with neighbors
    const auto& currparticles = typeIt.second;

    // iterate over particles of current type
    for (auto& particleIt : currparticles)
    {
      // get local index of particle i
      const int particle_i = particleIt.first;

      // declare pointer variables for particle i
      const double *rad_i, *mass_i, *wallnormal_i, *walldistance_i;
      double* interfacenormal_i;

      // get pointer to particle states
      rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
      mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
      interfacenormal_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::InterfaceNormal, particle_i);
      wallnormal_i = container_i->GetPtrToParticleState(PARTICLEENGINE::UnitWallNormal, particle_i);
      walldistance_i = container_i->GetPtrToParticleState(PARTICLEENGINE::WallDistance, particle_i);

      // evaluation only for particles close to boundary
      if (not(walldistance_i[0] < rad_i[0])) continue;

      // evaluation only for non-zero interface normal
      const double interfacenormal_i_norm = std::sqrt(interfacenormal_i[0] * interfacenormal_i[0] +
                                                      interfacenormal_i[1] * interfacenormal_i[1] +
                                                      interfacenormal_i[2] * interfacenormal_i[2]);
      if (not(interfacenormal_i_norm > 0.0)) continue;

      // initial particle spacing
      const double initspacing_i =
          std::pow((mass_i[0] / material_i->initDensity_), (1.0 / kernelspacedim));

      // corrected wall distance and maximum correction distance
      const double dw_i = walldistance_i[0] - initspacing_i;
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
      double sclrprdct = interfacenormal_i[0] * wallnormal_i[0] +
                         interfacenormal_i[1] * wallnormal_i[1] +
                         interfacenormal_i[2] * wallnormal_i[2];
      for (int i = 0; i < 3; ++i)
        walltangential_i[i] = interfacenormal_i[i] - sclrprdct * wallnormal_i[i];

      // norm of wall tangential
      const double walltangential_i_norm = std::sqrt(walltangential_i[0] * walltangential_i[0] +
                                                     walltangential_i[1] * walltangential_i[1] +
                                                     walltangential_i[2] * walltangential_i[2]);

      // scale unit wall tangential
      for (int i = 0; i < 3; ++i) walltangential_i[i] = walltangential_i[i] / walltangential_i_norm;

      // convert static contact angle in radians
      double theta_0 = staticcontactangle_ * M_PI / 180.0;
      if (type_i == PARTICLEENGINE::Phase2) theta_0 = (180 - staticcontactangle_) * M_PI / 180.0;

      // determine triple point normal
      double triplepointnormal_i[3];
      for (int i = 0; i < 3; ++i)
        triplepointnormal_i[i] =
            std::sin(theta_0) * walltangential_i[i] + std::cos(theta_0) * wallnormal_i[i];

      // determine corrected normal
      double correctednormal_i[3];
      for (int i = 0; i < 3; ++i)
        correctednormal_i[i] = f_i * interfacenormal_i[i] + (1.0 - f_i) * triplepointnormal_i[i];

      // norm of corrected normal
      const double correctednormal_i_norm = std::sqrt(correctednormal_i[0] * correctednormal_i[0] +
                                                      correctednormal_i[1] * correctednormal_i[1] +
                                                      correctednormal_i[2] * correctednormal_i[2]);

      // scale corrected normal
      for (int i = 0; i < 3; ++i)
        correctednormal_i[i] = correctednormal_i[i] / correctednormal_i_norm;

      // overwrite interface normal
      for (int i = 0; i < 3; ++i) interfacenormal_i[i] = correctednormal_i[i];
    }
  }
}

/*---------------------------------------------------------------------------*
 | refresh colorfield gradient and interface normal           sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::
    RefreshColorfieldGradientAndInterfaceNormal() const
{
  // init map
  std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>> particlestatestotypes;

  // iterate over particle types
  for (auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // no refreshing of surface tension states for boundary or rigid particles
    if (typeEnum == PARTICLEENGINE::BoundaryPhase or typeEnum == PARTICLEENGINE::RigidPhase)
      continue;

    // set state enums to map
    particlestatestotypes[typeEnum] = {
        PARTICLEENGINE::ColorfieldGradient, PARTICLEENGINE::InterfaceNormal};
  }

  // refresh specific states of particles of specific types
  particleengineinterface_->RefreshSpecificStatesOfParticlesOfSpecificTypes(particlestatestotypes);
}

/*---------------------------------------------------------------------------*
 | compute curvature and add acceleration contribution        sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHSurfaceTensionContinuumSurfaceForce::
    ComputeCurvatureAndAddAccelerationContribution() const
{
  // evaluate surface tension ramp function
  double timefac = 1.0;
  if (surfacetensionrampfctnumber_ > 0)
    timefac = DRT::Problem::Instance()->Funct(surfacetensionrampfctnumber_ - 1).EvaluateTime(time_);

  // iterate over particle types
  for (auto& typeIt : neighborpairs_->GetRefToNeighborPairsMap())
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type_i = typeIt.first;

    // no curvature evaluation for boundary or rigid particles
    if (type_i == PARTICLEENGINE::BoundaryPhase or type_i == PARTICLEENGINE::RigidPhase) continue;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // clear curvature state
    container_i->ClearState(PARTICLEENGINE::Curvature);

    // particles of current type with neighbors
    const auto& currparticles = typeIt.second;

    // iterate over particles of current type
    for (auto& particleIt : currparticles)
    {
      // get local index of particle i
      const int particle_i = particleIt.first;

      // declare pointer variables for particle i
      const double *rad_i, *mass_i, *dens_i, *colorfieldgrad_i, *interfacenormal_i;
      double *acc_i, *curvature_i;

      // get pointer to particle states
      rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
      mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
      dens_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Density, particle_i);
      interfacenormal_i =
          container_i->GetPtrToParticleState(PARTICLEENGINE::InterfaceNormal, particle_i);

      // evaluation only for non-zero interface normal
      const double interfacenormal_i_norm = std::sqrt(interfacenormal_i[0] * interfacenormal_i[0] +
                                                      interfacenormal_i[1] * interfacenormal_i[1] +
                                                      interfacenormal_i[2] * interfacenormal_i[2]);
      if (not(interfacenormal_i_norm > 0.0)) continue;

      // initialize sum of evaluated kernel values for particle i due to neighbor particles j
      double sumj_nij_Vj_eij_dWij(0.0);
      double sumj_Vj_Wij(0.0);

      // evaluate kernel
      const double Wii = kernel_->W(0.0, rad_i[0]);

      // add self-interaction
      sumj_Vj_Wij += Wii * mass_i[0] / dens_i[0];

      // iterate over particle types of neighboring particles
      for (auto& neighborTypeIt : particleIt.second)
      {
        // get type of neighboring particles
        PARTICLEENGINE::TypeEnum type_j = neighborTypeIt.first;

        // no evaluation for neighboring boundary or rigid particles
        if (type_j == PARTICLEENGINE::BoundaryPhase or type_j == PARTICLEENGINE::RigidPhase)
          continue;

        // change sign of interface normal for different particle types
        double signfac = (type_i == type_j) ? 1.0 : -1.0;

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
            const double *mass_j, *dens_j, *interfacenormal_j;

            // get pointer to particle states
            mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);
            dens_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Density, particle_j);
            interfacenormal_j =
                container_j->GetPtrToParticleState(PARTICLEENGINE::InterfaceNormal, particle_j);

            // evaluation only for non-zero interface normal
            const double interfacenormal_j_norm =
                std::sqrt(interfacenormal_j[0] * interfacenormal_j[0] +
                          interfacenormal_j[1] * interfacenormal_j[1] +
                          interfacenormal_j[2] * interfacenormal_j[2]);
            if (not(interfacenormal_j_norm > 0.0)) continue;

            // volume of particle j
            const double V_j = mass_j[0] / dens_j[0];

            double n_ij[3];
            for (int i = 0; i < 3; ++i)
              n_ij[i] = interfacenormal_i[i] - signfac * interfacenormal_j[i];

            // initial curvature estimate
            sumj_nij_Vj_eij_dWij +=
                (n_ij[0] * particlepair.e_ij_[0] + n_ij[1] * particlepair.e_ij_[1] +
                    n_ij[2] * particlepair.e_ij_[2]) *
                V_j * particlepair.dWdrij_;

            // correction factor
            sumj_Vj_Wij += V_j * particlepair.Wij_;
          }
        }
      }

      // only add meaningful contributions
      if (std::abs(sumj_Vj_Wij) > (1.0e-10 * rad_i[0]))
      {
        // get pointer to particle states
        acc_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Acceleration, particle_i);
        colorfieldgrad_i =
            container_i->GetPtrToParticleState(PARTICLEENGINE::ColorfieldGradient, particle_i);
        curvature_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Curvature, particle_i);

        // compute curvature
        curvature_i[0] = -sumj_nij_Vj_eij_dWij / sumj_Vj_Wij;

        // add contribution to acceleration
        const double fac = -timefac * surfacetensioncoefficient_ * curvature_i[0] / mass_i[0];
        for (int i = 0; i < 3; ++i) acc_i[i] += fac * colorfieldgrad_i[i];
      }
    }
  }
}
