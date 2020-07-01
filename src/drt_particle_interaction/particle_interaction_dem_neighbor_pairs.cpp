/*---------------------------------------------------------------------------*/
/*! \file
\brief neighbor pair handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_dem_neighbor_pairs.H"

#include "particle_interaction_utils.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_particle_wall/particle_wall_interface.H"

#include "../drt_mat/particle_wall_material_dem.H"

#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/element_coordtrafo.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMNeighborPairs::DEMNeighborPairs()
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMNeighborPairs::Init()
{
  // nothing to do
}

void PARTICLEINTERACTION::DEMNeighborPairs::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // set interface to particle wall handler
  particlewallinterface_ = particlewallinterface;
}

void PARTICLEINTERACTION::DEMNeighborPairs::WriteRestart(const int step, const double time) const
{
  // nothing to do
}

void PARTICLEINTERACTION::DEMNeighborPairs::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

void PARTICLEINTERACTION::DEMNeighborPairs::EvaluateNeighborPairs()
{
  // evaluate particle pairs
  EvaluateParticlePairs();

  // evaluate particle-wall pairs
  if (particlewallinterface_) EvaluateParticleWallPairs();
}

void PARTICLEINTERACTION::DEMNeighborPairs::EvaluateNeighborPairsAdhesion(
    const double& adhesion_distance)
{
  // evaluate adhesion particle pairs
  EvaluateParticlePairsAdhesion(adhesion_distance);

  // evaluate adhesion particle-wall pairs
  if (particlewallinterface_) EvaluateParticleWallPairsAdhesion(adhesion_distance);
}

void PARTICLEINTERACTION::DEMNeighborPairs::EvaluateParticlePairs()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::DEMNeighborPairs::EvaluateParticlePairs");

  // clear particle pair data
  particlepairdata_.clear();

  // iterate over potential particle neighbors
  for (const auto& potentialneighbors : particleengineinterface_->GetPotentialParticleNeighbors())
  {
    // access values of local index tuples of particle i and j
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = potentialneighbors.first;

    PARTICLEENGINE::TypeEnum type_j;
    PARTICLEENGINE::StatusEnum status_j;
    int particle_j;
    std::tie(type_j, status_j, particle_j) = potentialneighbors.second;

    // get corresponding particle containers
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    PARTICLEENGINE::ParticleContainer* container_j =
        particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

    // declare pointer variables for particle i and j
    const double *pos_i, *rad_i, *mass_i;
    const double *pos_j, *rad_j, *mass_j;

    // get pointer to particle states
    pos_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i);
    rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
    mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);

    pos_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Position, particle_j);
    rad_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_j);
    mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);

    // vector from particle i to j
    double r_ji[3];

    // distance between particles considering periodic boundaries
    particleengineinterface_->DistanceBetweenParticles(pos_i, pos_j, r_ji);

    // absolute distance between particles
    const double absdist = UTILS::vec_norm2(r_ji);

    // gap between particles
    const double gap = absdist - rad_i[0] - rad_j[0];

    // neighboring particles within interaction distance
    if (gap < 0.0)
    {
      // initialize particle pair
      particlepairdata_.push_back(DEMParticlePair());

      // get reference to current particle pair
      DEMParticlePair& particlepair = particlepairdata_.back();

      // set local index tuple of particles i and j
      particlepair.tuple_i_ = potentialneighbors.first;
      particlepair.tuple_j_ = potentialneighbors.second;

      // set gap between particles
      particlepair.gap_ = gap;

      // versor from particle i to j
      UTILS::vec_setscale(particlepair.e_ji_, (1.0 / absdist), r_ji);

      // set effective mass of particles i and j
      particlepair.m_eff_ = mass_i[0] * mass_j[0] / (mass_i[0] + mass_j[0]);
    }
  }
}

void PARTICLEINTERACTION::DEMNeighborPairs::EvaluateParticleWallPairs()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::DEMNeighborPairs::EvaluateParticleWallPairs");

  // clear particle-wall pair data
  particlewallpairdata_.clear();

  // relate particles to index of particle-wall pairs (considering object type of contact point)
  std::unordered_map<int, std::vector<std::pair<GEO::ObjectType, int>>>
      particletoindexofparticlewallpairs;

  // index of particle-wall pairs
  int particlewallpairindex = 0;

  // iterate over potential wall neighbors
  for (const auto& potentialneighbors : particlewallinterface_->GetPotentialWallNeighbors())
  {
    // access values of local index tuple of particle i
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = potentialneighbors.first;

    // get corresponding particle container
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    // get global id of particle i
    const int* globalid_i = container_i->GetPtrToParticleGlobalID(particle_i);

    // declare pointer variables for particle i
    const double* rad_i;

    // get pointer to particle states
    rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);

    // get position of particle i
    const LINALG::Matrix<3, 1> pos_i(
        container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i));

    // get pointer to column wall element
    DRT::Element* ele = potentialneighbors.second;

    // determine nodal positions of column wall element
    std::map<int, LINALG::Matrix<3, 1>> colelenodalpos;
    particlewallinterface_->DetermineColWallEleNodalPos(ele, colelenodalpos);

    // get coordinates of closest point on current column wall element to particle
    LINALG::Matrix<3, 1> closestpos;
    GEO::ObjectType objecttype =
        GEO::nearest3DObjectOnElement(ele, colelenodalpos, pos_i, closestpos);

    // vector from particle i to wall contact point j
    double r_ji[3];
    for (int i = 0; i < 3; i++) r_ji[i] = closestpos(i) - pos_i(i);

    // absolute distance between particle and wall contact point
    const double absdist = UTILS::vec_norm2(r_ji);

    // gap between particle and wall contact point
    const double gap = absdist - rad_i[0];

    // neighboring particle and wall element within interaction distance
    if (gap < 0.0)
    {
      // initialize particle-wall pair
      particlewallpairdata_.push_back(DEMParticleWallPair());

      // get reference to current particle-wall pair
      DEMParticleWallPair& particlewallpair = particlewallpairdata_[particlewallpairindex];

      // store index of particle-wall pair
      particletoindexofparticlewallpairs[globalid_i[0]].push_back(
          std::make_pair(objecttype, particlewallpairindex));

      // increase index
      ++particlewallpairindex;

      // set local index tuple of particle i
      particlewallpair.tuple_i_ = potentialneighbors.first;

      // set pointer to column wall element
      particlewallpair.ele_ = potentialneighbors.second;

      // set gap between particle and wall contact point
      particlewallpair.gap_ = gap;

      // versor from particle i to wall contact point j
      UTILS::vec_setscale(particlewallpair.e_ji_, (1.0 / absdist), r_ji);

      // get coordinates of wall contact point in element parameter space
      LINALG::Matrix<2, 1> elecoords(true);
      const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(ele, colelenodalpos));
      GEO::CurrentToSurfaceElementCoordinates(ele->Shape(), xyze, closestpos, elecoords);

      // set parameter space coordinates of wall contact point
      particlewallpair.elecoords_[0] = elecoords(0, 0);
      particlewallpair.elecoords_[1] = elecoords(1, 0);
    }
  }

  // set of particle-wall pairs to remove
  std::set<int> particlewallpairstoremove;

  // iterate over particles with neighboring wall contact points
  for (auto& particleIt : particletoindexofparticlewallpairs)
  {
    // get reference to index of particle-wall pairs for current particle
    std::vector<std::pair<GEO::ObjectType, int>>& indexofparticlewallpairs = particleIt.second;

    // only one particle-wall pair for current particle
    if (indexofparticlewallpairs.size() == 1) continue;

    // get local index tuple of current particle
    PARTICLEENGINE::LocalIndexTuple tuple_i =
        particlewallpairdata_[indexofparticlewallpairs[0].second].tuple_i_;

    // access values of local index tuple of particle i
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = tuple_i;

    // get corresponding particle container
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    // declare pointer variables for particle i
    const double* rad_i;

    // get pointer to particle states
    rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);

    // define tolerance dependent on the particle radius
    const double adaptedtol = 1.0e-7 * rad_i[0];

    // iterate over particle-wall pairs (master)
    for (std::pair<GEO::ObjectType, int>& master : indexofparticlewallpairs)
    {
      // get reference to particle-wall pair (master)
      DEMParticleWallPair& masterpair = particlewallpairdata_[master.second];

      // intersection radius of particle with column wall element in wall contact point
      const double intersectionradius =
          std::sqrt(UTILS::pow<2>(rad_i[0]) - UTILS::pow<2>(rad_i[0] + masterpair.gap_));

      // check with other particle-wall pairs (slave)
      for (std::pair<GEO::ObjectType, int>& slave : indexofparticlewallpairs)
      {
        // no-self checking
        if (master.second == slave.second) continue;

        // get reference to particle-wall pair (slave)
        DEMParticleWallPair& slavepair = particlewallpairdata_[slave.second];

        // vector between detected wall contact points
        double dist[3];
        UTILS::vec_setscale(dist, (rad_i[0] + masterpair.gap_), masterpair.e_ji_);
        UTILS::vec_addscale(dist, -(rad_i[0] + slavepair.gap_), slavepair.e_ji_);

        // absolute distance between wall contact points
        const double absdist = UTILS::vec_norm2(dist);

        bool removeslavepair = false;

        // check for coincident contact points of same type (e.g. on line between two surfaces)
        if (master.first == slave.first)
        {
          // contact point already detected (e.g. on line between two surfaces)
          if (absdist <= adaptedtol)
          {
            if (master.second < slave.second) removeslavepair = true;
          }
        }
        // check for line/node contact points within penetration volume of a surface contact point
        else if (master.first == GEO::SURFACE_OBJECT)
        {
          if (absdist <= intersectionradius) removeslavepair = true;
        }
        // check for node contact points within penetration volume of a line contact point
        else if (master.first == GEO::LINE_OBJECT and slave.first == GEO::NODE_OBJECT)
        {
          if (absdist <= intersectionradius) removeslavepair = true;
        }

        if (removeslavepair)
        {
          // mark particle-wall pair (slave) to be removed
          particlewallpairstoremove.insert(slave.second);
          // add global id of slave wall element for interaction history
          masterpair.histeles_.insert(slavepair.ele_->Id());
        }
      }
    }
  }

  // erase particle-wall pairs to be removed
  {
    int numparticlewallpairs = particlewallpairdata_.size();

    std::set<int>::reverse_iterator rit;
    for (rit = particlewallpairstoremove.rbegin(); rit != particlewallpairstoremove.rend(); ++rit)
      particlewallpairdata_[*rit] = particlewallpairdata_[--numparticlewallpairs];

    particlewallpairdata_.resize(numparticlewallpairs);
  }
}

void PARTICLEINTERACTION::DEMNeighborPairs::EvaluateParticlePairsAdhesion(
    const double& adhesion_distance)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::DEMNeighborPairs::EvaluateParticlePairsAdhesion");

  // clear adhesion particle pair data
  particlepairadhesiondata_.clear();

  // iterate over potential particle neighbors
  for (const auto& potentialneighbors : particleengineinterface_->GetPotentialParticleNeighbors())
  {
    // access values of local index tuples of particle i and j
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = potentialneighbors.first;

    PARTICLEENGINE::TypeEnum type_j;
    PARTICLEENGINE::StatusEnum status_j;
    int particle_j;
    std::tie(type_j, status_j, particle_j) = potentialneighbors.second;

    // get corresponding particle containers
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    PARTICLEENGINE::ParticleContainer* container_j =
        particlecontainerbundle_->GetSpecificContainer(type_j, status_j);

    // declare pointer variables for particle i and j
    const double *pos_i, *rad_i, *mass_i;
    const double *pos_j, *rad_j, *mass_j;

    // get pointer to particle states
    pos_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i);
    rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);
    mass_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_i);

    pos_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Position, particle_j);
    rad_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_j);
    mass_j = container_j->GetPtrToParticleState(PARTICLEENGINE::Mass, particle_j);

    // vector from particle i to j
    double r_ji[3];

    // distance between particles considering periodic boundaries
    particleengineinterface_->DistanceBetweenParticles(pos_i, pos_j, r_ji);

    // absolute distance between particles
    const double absdist = UTILS::vec_norm2(r_ji);

    // gap between particles
    const double gap = absdist - rad_i[0] - rad_j[0];

    // neighboring particles within adhesion distance
    if (gap < adhesion_distance)
    {
      // initialize particle pair
      particlepairadhesiondata_.push_back(DEMParticlePair());

      // get reference to current particle pair
      DEMParticlePair& particlepair = particlepairadhesiondata_.back();

      // set local index tuple of particles i and j
      particlepair.tuple_i_ = potentialneighbors.first;
      particlepair.tuple_j_ = potentialneighbors.second;

      // set gap between particles
      particlepair.gap_ = gap;

      // versor from particle i to j
      UTILS::vec_setscale(particlepair.e_ji_, (1.0 / absdist), r_ji);

      // set effective mass of particles i and j
      particlepair.m_eff_ = mass_i[0] * mass_j[0] / (mass_i[0] + mass_j[0]);
    }
  }
}

void PARTICLEINTERACTION::DEMNeighborPairs::EvaluateParticleWallPairsAdhesion(
    const double& adhesion_distance)
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::DEMNeighborPairs::EvaluateParticleWallPairsAdhesion");

  // clear particle-wall pair data
  particlewallpairadhesiondata_.clear();

  // relate particles to index of particle-wall pairs (considering object type of contact point)
  std::unordered_map<int, std::vector<std::pair<GEO::ObjectType, int>>>
      particletoindexofparticlewallpairs;

  // index of particle-wall pairs
  int particlewallpairindex = 0;

  // iterate over potential wall neighbors
  for (const auto& potentialneighbors : particlewallinterface_->GetPotentialWallNeighbors())
  {
    // access values of local index tuple of particle i
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = potentialneighbors.first;

    // get corresponding particle container
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    // get global id of particle i
    const int* globalid_i = container_i->GetPtrToParticleGlobalID(particle_i);

    // declare pointer variables for particle i
    const double* rad_i;

    // get pointer to particle states
    rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);

    // get position of particle i
    const LINALG::Matrix<3, 1> pos_i(
        container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i));

    // get pointer to column wall element
    DRT::Element* ele = potentialneighbors.second;

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

    // determine nodal positions of column wall element
    std::map<int, LINALG::Matrix<3, 1>> colelenodalpos;
    particlewallinterface_->DetermineColWallEleNodalPos(ele, colelenodalpos);

    // get coordinates of closest point on current column wall element to particle
    LINALG::Matrix<3, 1> closestpos;
    GEO::ObjectType objecttype =
        GEO::nearest3DObjectOnElement(ele, colelenodalpos, pos_i, closestpos);

    // vector from particle i to wall contact point j
    double r_ji[3];
    for (int i = 0; i < 3; i++) r_ji[i] = closestpos(i) - pos_i(i);

    // absolute distance between particle and wall contact point
    const double absdist = UTILS::vec_norm2(r_ji);

    // gap between particle and wall contact point
    const double gap = absdist - rad_i[0];

    // neighboring particle and wall element within adhesion distance
    if (gap < adhesion_distance)
    {
      // initialize particle-wall pair
      particlewallpairadhesiondata_.push_back(DEMParticleWallPair());

      // get reference to current particle-wall pair
      DEMParticleWallPair& particlewallpair = particlewallpairadhesiondata_[particlewallpairindex];

      // store index of particle-wall pair
      particletoindexofparticlewallpairs[globalid_i[0]].push_back(
          std::make_pair(objecttype, particlewallpairindex));

      // increase index
      ++particlewallpairindex;

      // set local index tuple of particle i
      particlewallpair.tuple_i_ = potentialneighbors.first;

      // set pointer to column wall element
      particlewallpair.ele_ = potentialneighbors.second;

      // set gap between particle and wall contact point
      particlewallpair.gap_ = gap;

      // versor from particle i to wall contact point j
      UTILS::vec_setscale(particlewallpair.e_ji_, (1.0 / absdist), r_ji);

      // get coordinates of wall contact point in element parameter space
      LINALG::Matrix<2, 1> elecoords(true);
      const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(ele, colelenodalpos));
      GEO::CurrentToSurfaceElementCoordinates(ele->Shape(), xyze, closestpos, elecoords);

      // set parameter space coordinates of wall contact point
      particlewallpair.elecoords_[0] = elecoords(0, 0);
      particlewallpair.elecoords_[1] = elecoords(1, 0);
    }
  }

  // set of particle-wall pairs to remove
  std::set<int> particlewallpairstoremove;

  // iterate over particles with neighboring wall contact points
  for (auto& particleIt : particletoindexofparticlewallpairs)
  {
    // get reference to index of particle-wall pairs for current particle
    std::vector<std::pair<GEO::ObjectType, int>>& indexofparticlewallpairs = particleIt.second;

    // only one particle-wall pair for current particle
    if (indexofparticlewallpairs.size() == 1) continue;

    // get local index tuple of current particle
    PARTICLEENGINE::LocalIndexTuple tuple_i =
        particlewallpairadhesiondata_[indexofparticlewallpairs[0].second].tuple_i_;

    // access values of local index tuple of particle i
    PARTICLEENGINE::TypeEnum type_i;
    PARTICLEENGINE::StatusEnum status_i;
    int particle_i;
    std::tie(type_i, status_i, particle_i) = tuple_i;

    // get corresponding particle container
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, status_i);

    // declare pointer variables for particle i
    const double* rad_i;

    // get pointer to particle states
    rad_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Radius, particle_i);

    // define tolerance dependent on the particle radius
    const double adaptedtol = 1.0e-7 * rad_i[0];

    // iterate over particle-wall pairs (master)
    for (std::pair<GEO::ObjectType, int>& master : indexofparticlewallpairs)
    {
      // get reference to particle-wall pair (master)
      DEMParticleWallPair& masterpair = particlewallpairadhesiondata_[master.second];

      // intersection radius of particle with column wall element in wall contact point
      const double intersectionradius = std::sqrt(
          UTILS::pow<2>(rad_i[0] + adhesion_distance) - UTILS::pow<2>(rad_i[0] + masterpair.gap_));

      // check with other particle-wall pairs (slave)
      for (std::pair<GEO::ObjectType, int>& slave : indexofparticlewallpairs)
      {
        // no-self checking
        if (master.second == slave.second) continue;

        // get reference to particle-wall pair (slave)
        DEMParticleWallPair& slavepair = particlewallpairadhesiondata_[slave.second];

        // vector between detected wall contact points
        double dist[3];
        UTILS::vec_setscale(dist, (rad_i[0] + masterpair.gap_), masterpair.e_ji_);
        UTILS::vec_addscale(dist, -(rad_i[0] + slavepair.gap_), slavepair.e_ji_);

        // absolute distance between wall contact points
        const double absdist = UTILS::vec_norm2(dist);

        bool removeslavepair = false;

        // check for coincident contact points of same type (e.g. on line between two surfaces)
        if (master.first == slave.first)
        {
          // contact point already detected (e.g. on line between two surfaces)
          if (absdist <= adaptedtol)
          {
            if (master.second < slave.second) removeslavepair = true;
          }
        }
        // check for line/node contact points within penetration volume of a surface contact point
        else if (master.first == GEO::SURFACE_OBJECT)
        {
          if (absdist <= intersectionradius) removeslavepair = true;
        }
        // check for node contact points within penetration volume of a line contact point
        else if (master.first == GEO::LINE_OBJECT and slave.first == GEO::NODE_OBJECT)
        {
          if (absdist <= intersectionradius) removeslavepair = true;
        }

        if (removeslavepair)
        {
          // mark particle-wall pair (slave) to be removed
          particlewallpairstoremove.insert(slave.second);
          // add global id of slave wall element for interaction history
          masterpair.histeles_.insert(slavepair.ele_->Id());
        }
      }
    }
  }

  // erase particle-wall pairs to be removed
  {
    int numparticlewallpairs = particlewallpairadhesiondata_.size();

    std::set<int>::reverse_iterator rit;
    for (rit = particlewallpairstoremove.rbegin(); rit != particlewallpairstoremove.rend(); ++rit)
      particlewallpairadhesiondata_[*rit] = particlewallpairadhesiondata_[--numparticlewallpairs];

    particlewallpairadhesiondata_.resize(numparticlewallpairs);
  }
}
