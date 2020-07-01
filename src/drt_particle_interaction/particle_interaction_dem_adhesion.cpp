/*---------------------------------------------------------------------------*/
/*! \file
\brief adhesion handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_dem_adhesion.H"

#include "particle_interaction_runtime_vtp_writer.H"
#include "particle_interaction_utils.H"

#include "particle_interaction_dem_neighbor_pairs.H"
#include "particle_interaction_dem_history_pairs.H"
#include "particle_interaction_dem_adhesion_law.H"
#include "particle_interaction_dem_adhesion_surface_energy.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_particle_wall/particle_wall_interface.H"
#include "../drt_particle_wall/particle_wall_datastate.H"

#include "../drt_mat/particle_wall_material_dem.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"

#include "../drt_io/runtime_vtp_writer.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMAdhesion::DEMAdhesion(const Teuchos::ParameterList& params)
    : params_dem_(params),
      adhesion_distance_(params_dem_.get<double>("ADHESION_DISTANCE")),
      writeparticlewallinteraction_(
          DRT::INPUT::IntegralValue<int>(params_dem_, "WRITE_PARTICLE_WALL_INTERACTION"))
{
  // empty constructor
}

PARTICLEINTERACTION::DEMAdhesion::~DEMAdhesion()
{
  // note: destructor declaration here since at compile-time a complete type
  // of class T as used in class member std::unique_ptr<T> ptr_T_ is required
}

void PARTICLEINTERACTION::DEMAdhesion::Init()
{
  // init adhesion law handler
  InitAdhesionLawHandler();

  // init adhesion surface energy handler
  InitAdhesionSurfaceEnergyHandler();
}

void PARTICLEINTERACTION::DEMAdhesion::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface,
    const std::shared_ptr<PARTICLEINTERACTION::InteractionWriter> particleinteractionwriter,
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

  // set particle interaction writer
  particleinteractionwriter_ = particleinteractionwriter;

  // setup particle interaction writer
  SetupParticleInteractionWriter();

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

void PARTICLEINTERACTION::DEMAdhesion::WriteRestart(const int step, const double time) const
{
  // write restart of adhesion law handler
  adhesionlaw_->WriteRestart(step, time);

  // write restart of adhesion surface energy handler
  adhesionsurfaceenergy_->WriteRestart(step, time);
}

void PARTICLEINTERACTION::DEMAdhesion::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // read restart of adhesion law handler
  adhesionlaw_->ReadRestart(reader);

  // read restart of adhesion surface energy handler
  adhesionsurfaceenergy_->ReadRestart(reader);
}

void PARTICLEINTERACTION::DEMAdhesion::AddForceContribution()
{
  // evaluate particle adhesion contribution
  EvaluateParticleAdhesion();

  // evaluate particle-wall adhesion contribution
  if (particlewallinterface_) EvaluateParticleWallAdhesion();
}

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

void PARTICLEINTERACTION::DEMAdhesion::SetupParticleInteractionWriter()
{
  // register specific runtime vtp writer
  if (writeparticlewallinteraction_)
    particleinteractionwriter_->RegisterSpecificRuntimeVtpWriter("particle-wall-adhesion");
}

void PARTICLEINTERACTION::DEMAdhesion::EvaluateParticleAdhesion()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::DEMAdhesion::EvaluateParticleAdhesion");

  // get reference to particle adhesion history pair data
  DEMHistoryPairAdhesionData& adhesionhistorydata =
      historypairs_->GetRefToParticleAdhesionHistoryData();

  // adhesion surface energy
  const double surface_energy = params_dem_.get<double>("ADHESION_SURFACE_ENERGY");

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
      adhesionsurfaceenergy_->AdhesionSurfaceEnergy(
          surface_energy, adhesionhistory_ij.surface_energy_);

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

void PARTICLEINTERACTION::DEMAdhesion::EvaluateParticleWallAdhesion()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::DEMAdhesion::EvaluateParticleWallAdhesion");

  // get wall data state container
  std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate =
      particlewallinterface_->GetWallDataState();

  // get reference to particle-wall pair data
  const DEMParticleWallPairData& particlewallpairdata =
      neighborpairs_->GetRefToParticleWallPairAdhesionData();

  // get reference to particle-wall adhesion history pair data
  DEMHistoryPairAdhesionData& adhesionhistorydata =
      historypairs_->GetRefToParticleWallAdhesionHistoryData();

  // write interaction output
  const bool writeinteractionoutput =
      particleinteractionwriter_->GetCurrentWriteResultFlag() and writeparticlewallinteraction_;

  // init storage for interaction output
  std::vector<double> attackpoints;
  std::vector<double> adhesionforces;
  std::vector<double> normaldirection;
  std::vector<double> surfaceenergy;

  // prepare storage for interaction output
  if (writeinteractionoutput)
  {
    const int numparticlewallpairs = particlewallpairdata.size();

    attackpoints.reserve(3 * numparticlewallpairs);
    adhesionforces.reserve(3 * numparticlewallpairs);
    normaldirection.reserve(3 * numparticlewallpairs);
    surfaceenergy.reserve(numparticlewallpairs);
  }

  // iterate over particle-wall pairs
  for (const auto& particlewallpair : particlewallpairdata)
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
    const double *pos_i, *vel_i, *rad_i, *mass_i;
    double* force_i;

    // get pointer to particle states
    pos_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i);
    vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);
    rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
    mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
    force_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Force, particle_i);

    // get pointer to column wall element
    DRT::Element* ele = particlewallpair.ele_;

    // number of nodes of wall element
    const int numnodes = ele->NumNode();

    // shape functions and location vector of wall element
    Epetra_SerialDenseVector funct(numnodes);
    std::vector<int> lmele;

    if (walldatastate->GetVelCol() != Teuchos::null or
        walldatastate->GetForceCol() != Teuchos::null)
    {
      // evaluate shape functions of element at wall contact point
      DRT::UTILS::shape_function_2D(
          funct, particlewallpair.elecoords_[0], particlewallpair.elecoords_[1], ele->Shape());

      // get location vector of wall element
      lmele.reserve(numnodes * 3);
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      ele->LocationVector(
          *particlewallinterface_->GetWallDiscretization(), lmele, lmowner, lmstride);
    }

    // adhesion surface energy
    double surface_energy = 0.0;

    // get material parameters of wall element
    {
      // cast material to particle wall material
      const Teuchos::RCP<const MAT::ParticleWallMaterialDEM>& particlewallmaterial =
          Teuchos::rcp_dynamic_cast<const MAT::ParticleWallMaterialDEM>(ele->Material());
      if (particlewallmaterial == Teuchos::null)
        dserror("cast to MAT::ParticleWallMaterialDEM failed!");

      // get adhesion surface energy
      surface_energy = particlewallmaterial->AdhesionSurfaceEnergy();
    }

    // no evaluation of adhesion contribution
    if (not(surface_energy > 0.0)) continue;

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
      adhesionsurfaceenergy_->AdhesionSurfaceEnergy(
          surface_energy, adhesionhistory_ij.surface_energy_);

    // calculate adhesion force
    adhesionlaw_->AdhesionForce(particlewallpair.gap_, adhesionhistory_ij.surface_energy_, rad_i[0],
        vel_rel_normal, mass_i[0], adhesionhistory_ij.adhesion_force_);

    // add adhesion force contribution
    UTILS::vec_addscale(force_i, adhesionhistory_ij.adhesion_force_, particlewallpair.e_ji_);

    // copy history to relevant wall elements in penetration volume
    for (int histele : particlewallpair.histeles_)
      adhesionhistorydata[globalid_i[0]][histele] = touchedadhesionhistory_ij;

    // calculation of wall adhesion force
    double walladhesionforce[3] = {0.0};
    if (writeinteractionoutput or walldatastate->GetForceCol() != Teuchos::null)
    {
      UTILS::vec_setscale(
          walladhesionforce, -adhesionhistory_ij.adhesion_force_, particlewallpair.e_ji_);
    }

    // write interaction output
    if (writeinteractionoutput)
    {
      // compute vector from particle i to wall contact point j
      double r_ji[3];
      UTILS::vec_setscale(r_ji, (rad_i[0] + particlewallpair.gap_), particlewallpair.e_ji_);

      // calculate wall contact point
      double wallcontactpoint[3];
      UTILS::vec_set(wallcontactpoint, pos_i);
      UTILS::vec_add(wallcontactpoint, r_ji);

      // set wall attack point and states
      for (int dim = 0; dim < 3; ++dim) attackpoints.push_back(wallcontactpoint[dim]);
      for (int dim = 0; dim < 3; ++dim) adhesionforces.push_back(walladhesionforce[dim]);
      for (int dim = 0; dim < 3; ++dim) normaldirection.push_back(-particlewallpair.e_ji_[dim]);
      surfaceenergy.push_back(adhesionhistory_ij.surface_energy_);
    }

    // assemble adhesion force acting on wall element
    if (walldatastate->GetForceCol() != Teuchos::null)
    {
      // determine nodal forces
      double nodal_force[numnodes * 3];
      for (int node = 0; node < numnodes; ++node)
        for (int dim = 0; dim < 3; ++dim)
          nodal_force[node * 3 + dim] = funct[node] * walladhesionforce[dim];

      // assemble nodal forces
      const int err = walldatastate->GetMutableForceCol()->SumIntoGlobalValues(
          numnodes * 3, &nodal_force[0], &(lmele)[0]);
      if (err < 0) dserror("sum into Epetra_Vector failed!");
    }
  }

  if (writeinteractionoutput)
  {
    // get specific runtime vtp writer
    RuntimeVtpWriter* runtime_vtpwriter =
        particleinteractionwriter_->GetSpecificRuntimeVtpWriter("particle-wall-adhesion");

    // set wall attack points
    runtime_vtpwriter->ResetGeometry(attackpoints);

    // append states
    runtime_vtpwriter->AppendVisualizationPointDataVector(adhesionforces, 3, "adhesion force");
    runtime_vtpwriter->AppendVisualizationPointDataVector(normaldirection, 3, "normal direction");
    runtime_vtpwriter->AppendVisualizationPointDataVector(surfaceenergy, 1, "surface energy");
  }
}
