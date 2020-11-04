/*---------------------------------------------------------------------------*/
/*! \file
\brief virtual wall particle handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_virtual_wall_particle.H"

#include "particle_interaction_sph_kernel.H"
#include "particle_interaction_sph_neighbor_pairs.H"

#include "particle_interaction_utils.H"

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
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHVirtualWallParticle::SPHVirtualWallParticle(
    const Teuchos::ParameterList& params)
    : params_sph_(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHVirtualWallParticle::Init()
{
  // init with potential fluid particle types
  allfluidtypes_ = {PARTICLEENGINE::Phase1, PARTICLEENGINE::Phase2, PARTICLEENGINE::DirichletPhase,
      PARTICLEENGINE::NeumannPhase};
  intfluidtypes_ = {PARTICLEENGINE::Phase1, PARTICLEENGINE::Phase2, PARTICLEENGINE::NeumannPhase};
}

void PARTICLEINTERACTION::SPHVirtualWallParticle::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface,
    const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
    const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // set interface to particle wall handler
  particlewallinterface_ = particlewallinterface;

  // set kernel handler
  kernel_ = kernel;

  // set neighbor pair handler
  neighborpairs_ = neighborpairs;

  // safety check
  if (not particlewallinterface_)
    dserror("interface to particle wall handler required in virtual wall particle handler!");

  // update with actual fluid particle types
  const auto allfluidtypes = allfluidtypes_;
  for (const auto& type_i : allfluidtypes)
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i))
      allfluidtypes_.erase(type_i);

  const auto intfluidtypes = intfluidtypes_;
  for (const auto& type_i : intfluidtypes)
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i))
      intfluidtypes_.erase(type_i);
}

void PARTICLEINTERACTION::SPHVirtualWallParticle::InitRelativePositionsOfVirtualParticles(
    const double maxinteractiondistance)
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::SPHVirtualWallParticle::InitRelativePositionsOfVirtualParticles");

  // clear relative positions of virtual particles
  virtualparticles_.clear();

  // get initial particle spacing
  const double initialparticlespacing = params_sph_.get<double>("INITIALPARTICLESPACING");

  // number of particles per direction
  const int numparticleperdir = std::round(maxinteractiondistance / initialparticlespacing);

  // iterate over virtual particles
  for (int r = 0; r < numparticleperdir; ++r)
  {
    for (int s = (-numparticleperdir + 1); s < numparticleperdir; ++s)
    {
      for (int t = (-numparticleperdir + 1); t < numparticleperdir; ++t)
      {
        // relative position of current virtual particle
        std::vector<double> currvirtualparticle(3);

        currvirtualparticle[0] = (0.5 + r) * initialparticlespacing;
        currvirtualparticle[1] = s * initialparticlespacing;
        currvirtualparticle[2] = t * initialparticlespacing;

        // current virtual particle within support radius
        if (UTILS::vec_norm2(&currvirtualparticle[0]) > maxinteractiondistance) continue;

        // add to relative positions of virtual particles
        virtualparticles_.push_back(currvirtualparticle);
      }
    }
  }
}

void PARTICLEINTERACTION::SPHVirtualWallParticle::InitStatesAtWallContactPoints(
    std::vector<double>& gravity)
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::SPHVirtualWallParticle::InitStatesAtWallContactPoints");

  // get wall data state container
  std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate =
      particlewallinterface_->GetWallDataState();

  // get reference to particle-wall pair data
  const SPHParticleWallPairData& particlewallpairdata =
      neighborpairs_->GetRefToParticleWallPairData();

  // get number of particle-wall pairs
  const int numparticlewallpairs = particlewallpairdata.size();

  // allocate memory
  weightedpressure_.assign(numparticlewallpairs, 0.0);
  weightedpressuregradient_.assign(numparticlewallpairs, std::vector<double>(3, 0.0));
  weighteddistancevector_.assign(numparticlewallpairs, std::vector<double>(3, 0.0));
  weightedvelocity_.assign(numparticlewallpairs, std::vector<double>(3, 0.0));

  // get relevant particle wall pair indices for specific particle types
  std::vector<int> relindices;
  neighborpairs_->GetRelevantParticleWallPairIndices(intfluidtypes_, relindices);

  // iterate over relevant particle-wall pairs
  for (const int particlewallpairindex : relindices)
  {
    const SPHParticleWallPair& particlewallpair = particlewallpairdata[particlewallpairindex];

    // access values of local index tuple of particle i
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = particlewallpair.tuple_i_;

    // get corresponding particle container
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    // get pointer to particle states
    const double* pos_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i);

    // get pointer to wall contact point states
    const double* rad_j = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);

    // get pointer to column wall element
    DRT::Element* ele = particlewallpair.ele_;

    // number of nodes of wall element
    const int numnodes = ele->NumNode();

    // shape functions and location vector of wall element
    Epetra_SerialDenseVector funct(numnodes);
    std::vector<int> lmele;

    if (walldatastate->GetVelCol() != Teuchos::null or walldatastate->GetAccCol() != Teuchos::null)
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

    // acceleration of wall contact point j
    double acc_j[3] = {0.0};

    if (walldatastate->GetAccCol() != Teuchos::null)
    {
      // get nodal accelerations
      std::vector<double> nodal_acc(numnodes * 3);
      DRT::UTILS::ExtractMyValues(*walldatastate->GetAccCol(), nodal_acc, lmele);

      // determine acceleration of wall contact point j
      for (int node = 0; node < numnodes; ++node)
        for (int dim = 0; dim < 3; ++dim) acc_j[dim] += funct[node] * nodal_acc[node * 3 + dim];
    }

    // compute vector from wall contact point j to particle i
    double r_ij[3];
    UTILS::vec_setscale(r_ij, particlewallpair.absdist_, particlewallpair.e_ij_);

    // compute position of wall contact point j
    double pos_j[3];
    UTILS::vec_set(pos_j, pos_i);
    UTILS::vec_sub(pos_j, r_ij);

    // get particles within radius
    std::vector<PARTICLEENGINE::LocalIndexTuple> neighboringparticles;
    particleengineinterface_->GetParticlesWithinRadius(&pos_j[0], rad_j[0], neighboringparticles);

#ifdef DEBUG
    if (not(neighboringparticles.size() > 0))
      dserror("expected at least one neighboring particle for wall contact point!");
#endif

    double sumk_Wjk = 0.0;
    double sumk_press_k_Wjk = 0.0;
    double sumk_dens_k_Wjk = 0.0;
    double sumk_r_jk_Wjk[3] = {0.0};
    double sumk_vel_k_Wjk[3] = {0.0};

    // iterate over neighboring particles
    for (const auto& neighboringparticle : neighboringparticles)
    {
      // access values of local index tuple of particle k
      PARTICLEENGINE::TypeEnum type_k;
      PARTICLEENGINE::StatusEnum status_k;
      int particle_k;
      std::tie(type_k, status_k, particle_k) = neighboringparticle;

      // evaluation only for fluid particles
      if (not allfluidtypes_.count(type_k)) continue;

      // get container of particles of current particle type
      PARTICLEENGINE::ParticleContainer* container_k =
          particlecontainerbundle_->GetSpecificContainer(type_k, status_k);

      // get pointer to particle states
      const double* pos_k =
          container_k->GetPtrToParticleState(PARTICLEENGINE::Position, particle_k);
      const double* vel_k =
          container_k->GetPtrToParticleState(PARTICLEENGINE::Velocity, particle_k);
      const double* dens_k =
          container_k->GetPtrToParticleState(PARTICLEENGINE::Density, particle_k);
      const double* press_k =
          container_k->GetPtrToParticleState(PARTICLEENGINE::Pressure, particle_k);

      // vector from particle k to wall contact point j
      double r_jk[3];

      // distance between wall contact point and particle considering periodic boundaries
      particleengineinterface_->DistanceBetweenParticles(pos_k, pos_j, r_jk);

      // absolute distance between particles
      const double absdist = UTILS::vec_norm2(r_jk);

      // evaluate kernel
      const double Wjk = kernel_->W(absdist, rad_j[0]);

      // sum contribution of neighboring particle k
      sumk_Wjk += Wjk;
      sumk_press_k_Wjk += press_k[0] * Wjk;
      sumk_dens_k_Wjk += dens_k[0] * Wjk;
      UTILS::vec_addscale(sumk_r_jk_Wjk, Wjk, r_jk);
      UTILS::vec_addscale(sumk_vel_k_Wjk, Wjk, vel_k);
    }

#ifdef DEBUG
    if (not(sumk_Wjk > 0))
      dserror("expected at least one neighboring particle for wall contact point!");
#endif

    const double inv_sumk_Wjk = 1.0 / sumk_Wjk;

    // compute relative acceleration of wall contact point
    double relacc[3];
    UTILS::vec_set(relacc, &gravity[0]);
    UTILS::vec_sub(relacc, acc_j);

    // set weighted fluid particle pressure
    weightedpressure_[particlewallpairindex] = sumk_press_k_Wjk * inv_sumk_Wjk;

    // set weighted fluid particle pressure gradient
    UTILS::vec_setscale(&weightedpressuregradient_[particlewallpairindex][0],
        sumk_dens_k_Wjk * inv_sumk_Wjk, relacc);

    // set weighted fluid particle distance vector
    UTILS::vec_setscale(
        &weighteddistancevector_[particlewallpairindex][0], inv_sumk_Wjk, sumk_r_jk_Wjk);

    // set weighted fluid particle velocity
    UTILS::vec_setscale(&weightedvelocity_[particlewallpairindex][0], inv_sumk_Wjk, sumk_vel_k_Wjk);
  }
}
