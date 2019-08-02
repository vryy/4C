/*---------------------------------------------------------------------------*/
/*!
\brief adhesion handler for discrete element method (DEM) interactions

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_dem_adhesion.H"

#include "particle_interaction_utils.H"

#include "particle_interaction_dem_neighbor_pairs.H"
#include "particle_interaction_dem_history_pairs.H"
#include "particle_interaction_dem_adhesion_law.H"
#include "particle_interaction_dem_adhesion_surface_energy.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_particle_wall/particle_wall_interface.H"
#include "../drt_particle_wall/particle_wall_datastate.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMAdhesion::DEMAdhesion(const Teuchos::ParameterList& params)
    : params_dem_(params), adhesion_distance_(params_dem_.get<double>("ADHESION_DISTANCE"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | destructor                                                 sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMAdhesion::~DEMAdhesion()
{
  // note: destructor declaration here since at compile-time a complete type
  // of class T as used in class member std::unique_ptr<T> ptr_T_ is required
}

/*---------------------------------------------------------------------------*
 | init adhesion handler                                      sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMAdhesion::Init()
{
  // init adhesion law handler
  InitAdhesionLawHandler();

  // init adhesion surface energy handler
  InitAdhesionSurfaceEnergyHandler();
}

/*---------------------------------------------------------------------------*
 | setup adhesion handler                                     sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMAdhesion::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface,
    const std::shared_ptr<PARTICLEINTERACTION::DEMNeighborPairs> neighborpairs,
    const std::shared_ptr<PARTICLEINTERACTION::DEMHistoryPairs> historypairs,
    const double& k_normal)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // set interface to particle wall hander
  particlewallinterface_ = particlewallinterface;

  // set neighbor pair handler
  neighborpairs_ = neighborpairs;

  // set history pair handler
  historypairs_ = historypairs;

  // setup adhesion law handler
  adhesionlaw_->Setup(k_normal);

  // setup adhesion surface energy handler
  adhesionsurfaceenergy_->Setup();

  // safety check
  if (adhesion_distance_ < 0.0) dserror("negative adhesion distance!");
}

/*---------------------------------------------------------------------------*
 | write restart of adhesion handler                          sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMAdhesion::WriteRestart(const int step, const double time) const
{
  // write restart of adhesion law handler
  adhesionlaw_->WriteRestart(step, time);

  // write restart of adhesion surface energy handler
  adhesionsurfaceenergy_->WriteRestart(step, time);
}

/*---------------------------------------------------------------------------*
 | read restart of adhesion handler                           sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMAdhesion::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // read restart of adhesion law handler
  adhesionlaw_->ReadRestart(reader);

  // read restart of adhesion surface energy handler
  adhesionsurfaceenergy_->ReadRestart(reader);
}

/*---------------------------------------------------------------------------*
 | add adhesion contribution to force field                   sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMAdhesion::AddForceContribution()
{
  // evaluate particle adhesion contribution
  EvaluateParticleAdhesion();

  // evaluate particle-wall adhesion contribution
  if (particlewallinterface_) EvaluateParticleWallAdhesion();
}

/*---------------------------------------------------------------------------*
 | init adhesion law handler                                  sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMAdhesion::InitAdhesionLawHandler()
{
  // get type of adhesion law
  INPAR::PARTICLE::AdhesionLaw adhesionlaw =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::AdhesionLaw>(params_dem_, "ADHESIONLAW");

  // create adhesion law handler
  switch (adhesionlaw)
  {
    case INPAR::PARTICLE::AdhesionVdWDMT:
    {
      adhesionlaw_ = std::unique_ptr<PARTICLEINTERACTION::DEMAdhesionLawVdWDMT>(
          new PARTICLEINTERACTION::DEMAdhesionLawVdWDMT(params_dem_));
      break;
    }
    case INPAR::PARTICLE::AdhesionRegDMT:
    {
      adhesionlaw_ = std::unique_ptr<PARTICLEINTERACTION::DEMAdhesionLawRegDMT>(
          new PARTICLEINTERACTION::DEMAdhesionLawRegDMT(params_dem_));
      break;
    }
    default:
    {
      dserror("unknown adhesion law type!");
      break;
    }
  }

  // init adhesion law handler
  adhesionlaw_->Init();
}

/*---------------------------------------------------------------------------*
 | init adhesion surface energy handler                       sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMAdhesion::InitAdhesionSurfaceEnergyHandler()
{
  // get type of adhesion surface energy distribution
  INPAR::PARTICLE::SurfaceEnergyDistribution surfaceenergydistributiontype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::SurfaceEnergyDistribution>(
          params_dem_, "ADHESION_SURFACE_ENERGY_DISTRIBUTION");

  // create adhesion surface energy handler
  switch (surfaceenergydistributiontype)
  {
    case INPAR::PARTICLE::ConstantSurfaceEnergy:
    {
      adhesionsurfaceenergy_ =
          std::unique_ptr<PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyConstant>(
              new PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyConstant(params_dem_));
      break;
    }
    case INPAR::PARTICLE::NormalSurfaceEnergyDistribution:
    {
      adhesionsurfaceenergy_ =
          std::unique_ptr<PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionNormal>(
              new PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionNormal(params_dem_));
      break;
    }
    case INPAR::PARTICLE::LogNormalSurfaceEnergyDistribution:
    {
      adhesionsurfaceenergy_ =
          std::unique_ptr<PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionLogNormal>(
              new PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionLogNormal(params_dem_));
      break;
    }
    default:
    {
      dserror("unknown adhesion surface energy distribution type!");
      break;
    }
  }

  // init adhesion surface energy handler
  adhesionsurfaceenergy_->Init();
}

/*---------------------------------------------------------------------------*
 | evaluate particle adhesion contribution                    sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMAdhesion::EvaluateParticleAdhesion()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::DEMAdhesion::EvaluateParticleAdhesion");

  // get reference to particle adhesion history pair data
  DEMHistoryPairAdhesionData& adhesionhistorydata =
      historypairs_->GetRefToParticleAdhesionHistoryData();

  // iterate over particle pairs
  for (const auto& particlepair : neighborpairs_->GetRefToParticlePairAdhesionData())
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

    // get corresponding particle containers
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    PARTICLEENGINE::ParticleContainer* container_j =
        particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

    // get global ids of particle
    const int* globalid_i = container_i->GetPtrToParticleGlobalID(particle_i);
    const int* globalid_j = container_j->GetPtrToParticleGlobalID(particle_j);

    // declare pointer variables for particle i and j
    const double *vel_i, *rad_i;
    double* force_i;

    const double *vel_j, *rad_j;
    double* force_j;

    // get pointer to particle states
    vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);
    rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
    force_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Force, particle_i);

    vel_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_j);
    rad_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_j);
    force_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Force, particle_j);

    // relative velocity in contact point c between particle i and j (neglecting angular velocity)
    double vel_rel[3];
    UTILS::vec_set(vel_rel, vel_i);
    UTILS::vec_sub(vel_rel, vel_j);

    // magnitude of relative velocity in normal direction
    const double vel_rel_normal = UTILS::vec_dot(vel_rel, particlepair.e_ji_);

    // calculate effective radius
    const double r_eff = (rad_i[0] * rad_j[0]) / (rad_i[0] + rad_j[0]);

    // get reference to touched adhesion history
    TouchedDEMHistoryPairAdhesion& touchedadhesionhistory_ij =
        adhesionhistorydata[globalid_i[0]][globalid_j[0]];

    // mark adhesion history as touched
    touchedadhesionhistory_ij.first = true;

    // get reference to adhesion history
    DEMHistoryPairAdhesion& adhesionhistory_ij = touchedadhesionhistory_ij.second;

    // calculate adhesion surface energy
    if (not(adhesionhistory_ij.surface_energy_ > 0.0))
      adhesionsurfaceenergy_->AdhesionSurfaceEnergy(adhesionhistory_ij.surface_energy_);

    // calculate adhesion force
    adhesionlaw_->AdhesionForce(particlepair.gap_, adhesionhistory_ij.surface_energy_, r_eff,
        vel_rel_normal, particlepair.m_eff_, adhesionhistory_ij.adhesion_force_);

    // copy history from interaction pair ij to ji
    if (status_j == PARTICLEENGINE::Owned)
    {
      // get reference to touched adhesion history
      TouchedDEMHistoryPairAdhesion& touchedadhesionhistory_ji =
          adhesionhistorydata[globalid_j[0]][globalid_i[0]];

      // mark adhesion history as touched
      touchedadhesionhistory_ji.first = true;

      // get reference to adhesion history
      DEMHistoryPairAdhesion& adhesionhistory_ji = touchedadhesionhistory_ji.second;

      // set adhesion surface energy and adhesion force
      adhesionhistory_ji.surface_energy_ = adhesionhistory_ij.surface_energy_;
      adhesionhistory_ji.adhesion_force_ = adhesionhistory_ij.adhesion_force_;
    }

    // add adhesion force contribution
    UTILS::vec_addscale(force_i, adhesionhistory_ij.adhesion_force_, particlepair.e_ji_);
    if (status_j == PARTICLEENGINE::Owned)
      UTILS::vec_addscale(force_j, -adhesionhistory_ij.adhesion_force_, particlepair.e_ji_);
  }
}

/*---------------------------------------------------------------------------*
 | evaluate particle-wall adhesion contribution               sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMAdhesion::EvaluateParticleWallAdhesion()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::DEMAdhesion::EvaluateParticleWallAdhesion");

  // get wall data state container
  std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate =
      particlewallinterface_->GetWallDataState();

  // get reference to particle-wall adhesion history pair data
  DEMHistoryPairAdhesionData& adhesionhistorydata =
      historypairs_->GetRefToParticleWallAdhesionHistoryData();

  // iterate over particle-wall pairs
  for (const auto& particlewallpair : neighborpairs_->GetRefToParticleWallPairAdhesionData())
  {
    // access values of local index tuple of particle i
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = particlewallpair.tuple_i_;

    // get corresponding particle container
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    // get global id of particle
    const int* globalid_i = container_i->GetPtrToParticleGlobalID(particle_i);

    // declare pointer variables for particle i
    const double *vel_i, *rad_i, *mass_i;
    double* force_i;

    // get pointer to particle states
    vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);
    rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
    mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
    force_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Force, particle_i);

    // get pointer to column wall element
    DRT::Element* ele = particlewallpair.ele_;

    // number of nodes of wall element
    const int numnodes = ele->NumNode();

    // evaluate shape functions of element at wall contact point
    Epetra_SerialDenseVector funct(numnodes);
    DRT::UTILS::shape_function_2D(
        funct, particlewallpair.elecoords_[0], particlewallpair.elecoords_[1], ele->Shape());

    // get location vector of wall element
    std::vector<int> lmele;
    lmele.reserve(numnodes * 3);
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    ele->LocationVector(*particlewallinterface_->GetWallDiscretization(), lmele, lmowner, lmstride);

    // velocity of wall contact point j
    double vel_j[3] = {0.0};

    if (walldatastate->GetVelCol() != Teuchos::null)
    {
      // get nodal velocities
      std::vector<double> nodal_vel(numnodes * 3);
      DRT::UTILS::ExtractMyValues(*walldatastate->GetVelCol(), nodal_vel, lmele);

      // determine velocity of wall contact point j
      for (int node = 0; node < numnodes; ++node)
        for (int dim = 0; dim < 3; ++dim) vel_j[dim] += funct[node] * nodal_vel[node * 3 + dim];
    }

    // relative velocity in wall contact point j (neglecting angular velocity)
    double vel_rel[3];
    UTILS::vec_set(vel_rel, vel_i);
    UTILS::vec_sub(vel_rel, vel_j);

    // magnitude of relative velocity in normal direction
    const double vel_rel_normal = UTILS::vec_dot(vel_rel, particlewallpair.e_ji_);

    // get reference to touched adhesion history
    TouchedDEMHistoryPairAdhesion& touchedadhesionhistory_ij =
        adhesionhistorydata[globalid_i[0]][ele->Id()];

    // mark adhesion history as touched
    touchedadhesionhistory_ij.first = true;

    // get reference to adhesion history
    DEMHistoryPairAdhesion& adhesionhistory_ij = touchedadhesionhistory_ij.second;

    // calculate adhesion surface energy
    if (not(adhesionhistory_ij.surface_energy_ > 0.0))
      adhesionsurfaceenergy_->AdhesionSurfaceEnergy(adhesionhistory_ij.surface_energy_);

    // calculate adhesion force
    adhesionlaw_->AdhesionForce(particlewallpair.gap_, adhesionhistory_ij.surface_energy_, rad_i[0],
        vel_rel_normal, mass_i[0], adhesionhistory_ij.adhesion_force_);

    // add adhesion force contribution
    UTILS::vec_addscale(force_i, adhesionhistory_ij.adhesion_force_, particlewallpair.e_ji_);

    // copy history to relevant wall elements in penetration volume
    for (int histele : particlewallpair.histeles_)
      adhesionhistorydata[globalid_i[0]][histele] = touchedadhesionhistory_ij;

    // assemble adhesion force acting on wall element
    if (walldatastate->GetForceCol() != Teuchos::null)
    {
      // wall contact force
      double wallcontactforce[3];
      UTILS::vec_setscale(
          wallcontactforce, -adhesionhistory_ij.adhesion_force_, particlewallpair.e_ji_);

      // determine nodal forces
      double nodal_force[numnodes * 3];
      for (int node = 0; node < numnodes; ++node)
        for (int dim = 0; dim < 3; ++dim)
          nodal_force[node * 3 + dim] = funct[node] * wallcontactforce[dim];

      // assemble nodal forces
      const int err = walldatastate->GetMutableForceCol()->SumIntoGlobalValues(
          numnodes * 3, &nodal_force[0], &(lmele)[0]);
      if (err < 0) dserror("sum into Epetra_Vector failed!");
    }
  }
}
