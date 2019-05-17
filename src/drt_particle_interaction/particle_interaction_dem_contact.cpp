/*---------------------------------------------------------------------------*/
/*!

\brief contact handler for discrete element method (DEM) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_dem_contact.H"

#include "particle_interaction_material_handler.H"
#include "particle_interaction_utils.H"

#include "particle_interaction_dem_neighbor_pairs.H"
#include "particle_interaction_dem_history_pairs.H"
#include "particle_interaction_dem_contact_normal.H"
#include "particle_interaction_dem_contact_tangential.H"

#include "../drt_particle_algorithm/particle_wall_interface.H"
#include "../drt_particle_algorithm/particle_wall_datastate.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMContact::DEMContact(const Teuchos::ParameterList& params)
    : params_dem_(params), dt_(0.0)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | destructor                                                 sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMContact::~DEMContact()
{
  // note: destructor declaration here since at compile-time a complete type
  // of class T as used in class member std::unique_ptr<T> ptr_T_ is required
}

/*---------------------------------------------------------------------------*
 | init contact handler                                       sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::Init()
{
  // init normal contact handler
  InitNormalContactHandler();

  // init tangential contact handler
  InitTangentialContactHandler();
}

/*---------------------------------------------------------------------------*
 | setup contact handler                                      sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEALGORITHM::WallHandlerInterface> particlewallinterface,
    const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
    const std::shared_ptr<PARTICLEINTERACTION::DEMNeighborPairs> neighborpairs,
    const std::shared_ptr<PARTICLEINTERACTION::DEMHistoryPairs> historypairs)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // set interface to particle wall hander
  particlewallinterface_ = particlewallinterface;

  // set particle material handler
  particlematerial_ = particlematerial;

  // set neighbor pair handler
  neighborpairs_ = neighborpairs;

  // set history pair handler
  historypairs_ = historypairs;

  // get maximum density of all materials
  const double maxdensity = GetMaxDensityOfAllMaterials();

  // setup normal contact handler
  contactnormal_->Setup(maxdensity);

  // setup tangential contact handler
  if (contacttangential_) contacttangential_->Setup(contactnormal_->GetNormalContactStiffness());

  // safety check
  if (contacttangential_)
  {
    // get type of normal contact law
    INPAR::PARTICLE::NormalContact normalcontacttype =
        DRT::INPUT::IntegralValue<INPAR::PARTICLE::NormalContact>(params_dem_, "NORMALCONTACTLAW");

    if (normalcontacttype != INPAR::PARTICLE::NormalLinSpring and
        normalcontacttype != INPAR::PARTICLE::NormalLinSpringDamp)
      dserror("tangential contact law only valid with linear normal contact law!");
  }
}

/*---------------------------------------------------------------------------*
 | write restart of contact handler                           sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::WriteRestart(const int step, const double time) const
{
  // write restart of normal contact handler
  contactnormal_->WriteRestart(step, time);

  // write restart of tangential contact handler
  if (contacttangential_) contacttangential_->WriteRestart(step, time);
}

/*---------------------------------------------------------------------------*
 | read restart of contact handler                            sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // read restart of normal contact handler
  contactnormal_->ReadRestart(reader);

  // read restart of tangential contact handler
  if (contacttangential_) contacttangential_->ReadRestart(reader);
}

/*---------------------------------------------------------------------------*
 | set current step size                                      sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::SetCurrentStepSize(const double currentstepsize)
{
  dt_ = currentstepsize;

  // set current step size
  if (contacttangential_) contacttangential_->SetCurrentStepSize(currentstepsize);
}

/*---------------------------------------------------------------------------*
 | insert contact evaluation dependent states                 sfuchs 12/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::InsertParticleStatesOfParticleTypes(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
    const
{
  // iterate over particle types
  for (auto& typeIt : particlestatestotypes)
  {
    // set of particle states for current particle type
    std::set<PARTICLEENGINE::StateEnum>& particlestates = typeIt.second;

    // states for tangential contact evaluation scheme
    if (contacttangential_)
      particlestates.insert({PARTICLEENGINE::Moment, PARTICLEENGINE::AngularVelocity,
          PARTICLEENGINE::AngularAcceleration});
  }
}

/*---------------------------------------------------------------------------*
 | check critical time step (on this processor)               sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::CheckCriticalTimeStep() const
{
  // init value of minimum mass
  double minmass = 0.0;

  // iterate over particle types
  for (const auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(typeEnum, PARTICLEENGINE::Owned);

    // get minimum stored value of state
    double currminmass = container->GetMinValueOfState(PARTICLEENGINE::Mass);

    if ((not(minmass > 0.0)) or (currminmass < minmass)) minmass = currminmass;
  }

  // currently no particles in simulation domain
  if (minmass == 0.0) return;

  // get time critical contact stiffness
  const double k_tcrit = contactnormal_->GetTimeCriticalStiffness();

  // critical time step size based on particle-particle contact
  const double safety = 0.75;
  const double factor = contacttangential_ ? 0.22 : 0.34;
  const double dt_crit = safety * factor * std::sqrt(minmass / k_tcrit);

  // checks time step
  if (dt_ > dt_crit) dserror("time step %f larger than critical time step %f!", dt_, dt_crit);
}

/*---------------------------------------------------------------------------*
 | add contact contribution to force and moment field         sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::AddForceAndMomentContribution()
{
  // evaluate particle contact contribution
  EvaluateParticleContact();

  // evaluate particle-wall contact contribution
  if (particlewallinterface_) EvaluateParticleWallContact();
}

/*---------------------------------------------------------------------------*
 | init normal contact handler                                sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::InitNormalContactHandler()
{
  // get type of normal contact law
  INPAR::PARTICLE::NormalContact normalcontacttype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::NormalContact>(params_dem_, "NORMALCONTACTLAW");

  // create normal contact handler
  switch (normalcontacttype)
  {
    case INPAR::PARTICLE::NormalLinSpring:
    {
      contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLinearSpring>(
          new PARTICLEINTERACTION::DEMContactNormalLinearSpring(params_dem_));
      break;
    }
    case INPAR::PARTICLE::NormalLinSpringDamp:
    {
      contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp>(
          new PARTICLEINTERACTION::DEMContactNormalLinearSpringDamp(params_dem_));
      break;
    }
    case INPAR::PARTICLE::NormalHertz:
    {
      contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalHertz>(
          new PARTICLEINTERACTION::DEMContactNormalHertz(params_dem_));
      break;
    }
    case INPAR::PARTICLE::NormalLeeHerrmann:
    {
      contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalLeeHerrmann>(
          new PARTICLEINTERACTION::DEMContactNormalLeeHerrmann(params_dem_));
      break;
    }
    case INPAR::PARTICLE::NormalKuwabaraKono:
    {
      contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalKuwabaraKono>(
          new PARTICLEINTERACTION::DEMContactNormalKuwabaraKono(params_dem_));
      break;
    }
    case INPAR::PARTICLE::NormalTsuji:
    {
      contactnormal_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalTsuji>(
          new PARTICLEINTERACTION::DEMContactNormalTsuji(params_dem_));
      break;
    }
    default:
    {
      dserror("unknown normal contact law type!");
      break;
    }
  }

  // init normal contact handler
  contactnormal_->Init();
}

/*---------------------------------------------------------------------------*
 | init tangential contact handler                            sfuchs 12/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::InitTangentialContactHandler()
{
  // get type of tangential contact law
  INPAR::PARTICLE::TangentialContact tangentialcontacttype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::TangentialContact>(
          params_dem_, "TANGENTIALCONTACTLAW");

  // create tangential contact handler
  switch (tangentialcontacttype)
  {
    case INPAR::PARTICLE::NoTangentialContact:
    {
      contacttangential_ = std::unique_ptr<PARTICLEINTERACTION::DEMContactTangentialBase>(nullptr);
      break;
    }
    case INPAR::PARTICLE::TangentialLinSpringDamp:
    {
      contacttangential_ =
          std::unique_ptr<PARTICLEINTERACTION::DEMContactTangentialLinearSpringDamp>(
              new PARTICLEINTERACTION::DEMContactTangentialLinearSpringDamp(params_dem_));
      break;
    }
    default:
    {
      dserror("unknown tangential contact law type!");
      break;
    }
  }

  // init tangential contact handler
  if (contacttangential_) contacttangential_->Init();
}

/*---------------------------------------------------------------------------*
 | get maximum density of all materials                       sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
double PARTICLEINTERACTION::DEMContact::GetMaxDensityOfAllMaterials() const
{
  // init value of maximum density
  double maxdensity(0.0);

  // iterate over particle types
  for (const auto& typeEnum : particlecontainerbundle_->GetParticleTypes())
  {
    // get material for current particle type
    const MAT::PAR::ParticleMaterialBase* material =
        particlematerial_->GetPtrToParticleMatParameter(typeEnum);

    // compare to current maximum
    maxdensity = std::max(maxdensity, material->initDensity_);
  }

  return maxdensity;
}

/*---------------------------------------------------------------------------*
 | evaluate particle contact contribution                     sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::EvaluateParticleContact()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::DEMContact::EvaluateParticleContact");

  // get reference to particle tangential history pair data
  DEMHistoryPairTangentialData& tangentialhistorydata =
      historypairs_->GetRefToParticleTangentialHistoryData();

  // iterate over particle pairs
  for (const auto& particlepair : neighborpairs_->GetRefToParticlePairData())
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
    const double *vel_i, *rad_i, *angvel_i;
    double *force_i, *moment_i;

    const double *vel_j, *rad_j, *angvel_j;
    double *force_j, *moment_j;

    // get pointer to particle states
    vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);
    rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
    force_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Force, particle_i);

    vel_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_j);
    rad_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_j);
    force_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Force, particle_j);

    if (contacttangential_)
    {
      angvel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::AngularVelocity, particle_i);
      moment_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Moment, particle_i);

      angvel_j = container_j->GetPtrToParticleState(PARTICLEENGINE::AngularVelocity, particle_j);
      moment_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Moment, particle_j);
    }

    // compute vectors from particle i and j to contact point c
    double r_ci[3], r_cj[3];
    if (contacttangential_)
    {
      UTILS::vec_setscale(r_ci, (rad_i[0] + 0.5 * particlepair.gap_), particlepair.e_ji_);
      UTILS::vec_setscale(r_cj, -(rad_j[0] + 0.5 * particlepair.gap_), particlepair.e_ji_);
    }

    // relative velocity in contact point c between particle i and j
    double vel_rel[3];
    UTILS::vec_set(vel_rel, vel_i);
    UTILS::vec_sub(vel_rel, vel_j);
    if (contacttangential_)
    {
      UTILS::vec_addcross(vel_rel, angvel_i, r_ci);
      UTILS::vec_addcross(vel_rel, r_cj, angvel_j);
    }

    // magnitude of relative velocity in normal direction
    const double vel_rel_normal = UTILS::vec_dot(vel_rel, particlepair.e_ji_);

    // calculate normal contact force
    double normalcontactforce(0.0);
    contactnormal_->NormalContactForce(
        particlepair.gap_, rad_i, rad_j, vel_rel_normal, particlepair.m_eff_, normalcontactforce);

    // add normal contact force contribution
    UTILS::vec_addscale(force_i, normalcontactforce, particlepair.e_ji_);
    if (status_j == PARTICLEENGINE::Owned)
      UTILS::vec_addscale(force_j, -normalcontactforce, particlepair.e_ji_);

    // calculation of tangential contact force
    if (contacttangential_)
    {
      // get reference to touched tangential history
      TouchedDEMHistoryPairTangential& touchedtangentialhistory_ij =
          tangentialhistorydata[globalid_i[0]][globalid_j[0]];

      // mark tangential history as touched
      touchedtangentialhistory_ij.first = true;

      // get reference to tangential history
      DEMHistoryPairTangential& tangentialhistory_ij = touchedtangentialhistory_ij.second;

      // relative velocity in tangential direction
      double vel_rel_tangential[3];
      UTILS::vec_set(vel_rel_tangential, vel_rel);
      UTILS::vec_addscale(vel_rel_tangential, -vel_rel_normal, particlepair.e_ji_);

      // calculate tangential contact force
      double tangentialcontactforce[3];
      contacttangential_->TangentialContactForce(tangentialhistory_ij.gap_t_,
          tangentialhistory_ij.stick_, particlepair.e_ji_, vel_rel_tangential, particlepair.m_eff_,
          normalcontactforce, tangentialcontactforce);

      // copy history from interaction pair ij to ji
      if (status_j == PARTICLEENGINE::Owned)
      {
        // get reference to touched tangential history
        TouchedDEMHistoryPairTangential& touchedtangentialhistory_ji =
            tangentialhistorydata[globalid_j[0]][globalid_i[0]];

        // mark tangential history as touched
        touchedtangentialhistory_ji.first = true;

        // get reference to tangential history
        DEMHistoryPairTangential& tangentialhistory_ji = touchedtangentialhistory_ji.second;

        // set tangential gap and tangential stick flag
        UTILS::vec_setscale(tangentialhistory_ji.gap_t_, -1.0, tangentialhistory_ij.gap_t_);
        tangentialhistory_ji.stick_ = tangentialhistory_ij.stick_;
      }

      // add tangential contact force contribution
      UTILS::vec_add(force_i, tangentialcontactforce);
      if (status_j == PARTICLEENGINE::Owned) UTILS::vec_sub(force_j, tangentialcontactforce);

      // add tangential contact moment contribution
      UTILS::vec_addcross(moment_i, r_ci, tangentialcontactforce);
      if (status_j == PARTICLEENGINE::Owned)
        UTILS::vec_addcross(moment_j, tangentialcontactforce, r_cj);
    }
  }
}

/*---------------------------------------------------------------------------*
 | evaluate particle-wall contact contribution                sfuchs 05/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMContact::EvaluateParticleWallContact()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::DEMContact::EvaluateParticleWallContact");

  // get wall data state container
  std::shared_ptr<PARTICLEALGORITHM::WallDataState> walldatastate =
      particlewallinterface_->GetWallDataState();

  // get reference to particle-wall tangential history pair data
  DEMHistoryPairTangentialData& tangentialhistorydata =
      historypairs_->GetRefToParticleWallTangentialHistoryData();

  // iterate over particle-wall pairs
  for (const auto& particlewallpair : neighborpairs_->GetRefToParticleWallPairData())
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
    const double *vel_i, *rad_i, *mass_i, *angvel_i;
    double *force_i, *moment_i;

    // get pointer to particle states
    vel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_i);
    rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
    mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);
    force_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Force, particle_i);

    if (contacttangential_)
    {
      angvel_i = container_i->GetPtrToParticleState(PARTICLEENGINE::AngularVelocity, particle_i);
      moment_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Moment, particle_i);
    }

    // get pointer to column wall element
    DRT::Element* ele = particlewallpair.ele_;

    // declare zero radius for wall contact point j
    const double rad_j = 0.0;

    // velocity of wall contact point j
    double vel_j[3] = {0.0};

    if (walldatastate->GetVelCol() != Teuchos::null)
    {
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
      ele->LocationVector(
          *particlewallinterface_->GetWallDiscretization(), lmele, lmowner, lmstride);

      // get nodal velocities
      std::vector<double> nodal_vel(numnodes * 3);
      DRT::UTILS::ExtractMyValues(*walldatastate->GetVelCol(), nodal_vel, lmele);

      // determine velocity of wall contact point j
      for (int node = 0; node < numnodes; ++node)
        for (int dim = 0; dim < 3; ++dim) vel_j[dim] += funct[node] * nodal_vel[node * 3 + dim];
    }

    // compute vector from particle i to wall contact point j
    double r_ji[3];
    if (contacttangential_)
      UTILS::vec_setscale(r_ji, (rad_i[0] + particlewallpair.gap_), particlewallpair.e_ji_);

    // relative velocity in wall contact point j
    double vel_rel[3];
    UTILS::vec_set(vel_rel, vel_i);
    UTILS::vec_sub(vel_rel, vel_j);
    if (contacttangential_) UTILS::vec_addcross(vel_rel, angvel_i, r_ji);

    // magnitude of relative velocity in normal direction
    const double vel_rel_normal = UTILS::vec_dot(vel_rel, particlewallpair.e_ji_);

    // calculate normal contact force
    double normalcontactforce(0.0);
    contactnormal_->NormalContactForce(
        particlewallpair.gap_, rad_i, &rad_j, vel_rel_normal, mass_i[0], normalcontactforce);

    // add normal contact force contribution
    UTILS::vec_addscale(force_i, normalcontactforce, particlewallpair.e_ji_);

    // calculation of tangential contact force
    if (contacttangential_)
    {
      // get reference to touched tangential history
      TouchedDEMHistoryPairTangential& touchedtangentialhistory_ij =
          tangentialhistorydata[globalid_i[0]][ele->Id()];

      // mark tangential history as touched
      touchedtangentialhistory_ij.first = true;

      // get reference to tangential history
      DEMHistoryPairTangential& tangentialhistory_ij = touchedtangentialhistory_ij.second;

      // relative velocity in tangential direction
      double vel_rel_tangential[3];
      UTILS::vec_set(vel_rel_tangential, vel_rel);
      UTILS::vec_addscale(vel_rel_tangential, -vel_rel_normal, particlewallpair.e_ji_);

      // calculate tangential contact force
      double tangentialcontactforce[3];
      contacttangential_->TangentialContactForce(tangentialhistory_ij.gap_t_,
          tangentialhistory_ij.stick_, particlewallpair.e_ji_, vel_rel_tangential, mass_i[0],
          normalcontactforce, tangentialcontactforce);

      // add tangential contact force contribution
      UTILS::vec_add(force_i, tangentialcontactforce);

      // add tangential contact moment contribution
      UTILS::vec_addcross(moment_i, r_ji, tangentialcontactforce);

      // copy history to relevant wall elements in penetration volume
      for (int histele : particlewallpair.histeles_)
        tangentialhistorydata[globalid_i[0]][histele] = touchedtangentialhistory_ij;
    }
  }
}
