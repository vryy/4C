/*---------------------------------------------------------------------------*/
/*!
\file particle_dirichlet_bc.cpp

\brief dirichlet boundary condition handler for particle simulations

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_dirichlet_bc.H"

#include "particle_algorithm_utils.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_container_bundle.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_globalproblem.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::DirichletBoundaryConditionHandler::DirichletBoundaryConditionHandler(
    const Teuchos::ParameterList& params)
    : params_(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init dirichlet boundary condition handler                  sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::DirichletBoundaryConditionHandler::Init()
{
  // get control parameters for initial/boundary conditions
  const Teuchos::ParameterList& params_conditions =
      params_.sublist("INITIAL AND BOUNDARY CONDITIONS");

  // read parameters relating particle types to IDs
  PARTICLEALGORITHM::UTILS::ReadParamsTypesRelatedToIDs(
      params_conditions, "DIRICHLET_BOUNDARY_CONDITION", dirichletbctypetofunctid_);

  // iterate over particle types and insert into set
  for (auto& typeIt : dirichletbctypetofunctid_) typessubjectedtodirichletbc_.insert(typeIt.first);
}

/*---------------------------------------------------------------------------*
 | setup dirichlet boundary condition handler                 sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::DirichletBoundaryConditionHandler::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;
}

/*---------------------------------------------------------------------------*
 | write restart of dirichlet boundary condition handler      sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::DirichletBoundaryConditionHandler::WriteRestart(
    const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of dirichlet boundary condition handler       sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::DirichletBoundaryConditionHandler::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | insert dbc dependent states of all particle types          sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::DirichletBoundaryConditionHandler::InsertParticleStatesOfParticleTypes(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
    const
{
  // iterate over particle types subjected to dirichlet boundary conditions
  for (auto& particleType : typessubjectedtodirichletbc_)
  {
    // insert states for types subjected to dirichlet boundary conditions
    particlestatestotypes[particleType].insert(PARTICLEENGINE::ReferencePosition);
  }
}

/*---------------------------------------------------------------------------*
 | set particle reference position                            sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::DirichletBoundaryConditionHandler::SetParticleReferencePosition() const
{
  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // iterate over particle types subjected to dirichlet boundary conditions
  for (auto& particleType : typessubjectedtodirichletbc_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container =
        particlecontainerbundle->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // set particle reference position
    container->UpdateState(0.0, PARTICLEENGINE::ReferencePosition, 1.0, PARTICLEENGINE::Position);
  }
}

/*---------------------------------------------------------------------------*
 | evaluate dirichlet boundary condition                      sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::DirichletBoundaryConditionHandler::EvaluateDirichletBoundaryCondition(
    const double& evaltime, const bool evalpos, const bool evalvel, const bool evalacc) const
{
  // degree of maximal function derivative
  int deg = 0;
  if (evalacc)
    deg = 2;
  else if (evalvel)
    deg = 1;

  // get bounding box dimensions
  LINALG::Matrix<3, 2> xaabb = particleengineinterface_->XAABB();

  // get bin size
  const double* binsize = particleengineinterface_->BinSize();

  // init vector containing evaluated function and derivatives
  std::vector<double> functtimederiv(deg + 1);

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // iterate over particle types subjected to dirichlet boundary conditions
  for (auto& typeIt : dirichletbctypetofunctid_)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum particleType = typeIt.first;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container =
        particlecontainerbundle->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get id of function
    const int functid = typeIt.second;

    // get reference to function
    DRT::UTILS::Function& function = DRT::Problem::Instance()->Funct(functid - 1);

    // declare pointer variables
    const double* refpos;
    double *pos, *vel, *acc;

    // get pointer to particle states
    refpos = container->GetPtrToParticleState(PARTICLEENGINE::ReferencePosition, 0);
    if (evalpos) pos = container->GetPtrToParticleState(PARTICLEENGINE::Position, 0);
    if (evalvel) vel = container->GetPtrToParticleState(PARTICLEENGINE::Velocity, 0);
    if (evalacc) acc = container->GetPtrToParticleState(PARTICLEENGINE::Acceleration, 0);

    // get particle state dimension
    int statedim = container->GetParticleStateDim(PARTICLEENGINE::Position);

    // safety check
    if (statedim != function.NumberComponents())
      dserror("dimension of function defining dirichlet boundary condition not correct!");

    // iterate over owned particles of current type
    for (int i = 0; i < particlestored; ++i)
    {
      // iterate over spatial dimension
      for (int dim = 0; dim < statedim; ++dim)
      {
        // evaluate function, first and second time derivative
        functtimederiv =
            function.EvaluateTimeDerivative(dim, &(refpos[statedim * i]), evaltime, deg);

        // set position state
        if (evalpos)
        {
          // check for periodic boundary condition in current spatial direction
          if (particleengineinterface_->HavePBC(dim))
          {
            // periodic length in current spatial direction
            const double pbcdelta = particleengineinterface_->PBCDelta(dim);

            // get displacement from reference position canceling out multiples of periodic length
            // in current spatial direction
            double displacement = std::fmod(functtimederiv[0], pbcdelta);

            // get new position
            double newpos = refpos[statedim * i + dim] + displacement;

            // shift by periodic length if new position is close to the periodic boundary and old
            // position is on other end domain
            if ((newpos > (xaabb(dim, 1) - binsize[dim])) and
                (std::abs(newpos - pos[statedim * i + dim]) > binsize[dim]))
              pos[statedim * i + dim] = newpos - pbcdelta;
            else if ((newpos < (xaabb(dim, 0) + binsize[dim])) and
                     (std::abs(newpos - pos[statedim * i + dim]) > binsize[dim]))
              pos[statedim * i + dim] = newpos + pbcdelta;
            else
              pos[statedim * i + dim] = newpos;
          }
          // no periodic boundary conditions in current spatial direction
          else
            pos[statedim * i + dim] = refpos[statedim * i + dim] + functtimederiv[0];
        }

        // set velocity state
        if (evalvel) vel[statedim * i + dim] = functtimederiv[1];

        // set acceleration state
        if (evalacc) acc[statedim * i + dim] = functtimederiv[2];
      }
    }
  }
}
