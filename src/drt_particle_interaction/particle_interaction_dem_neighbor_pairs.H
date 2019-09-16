/*---------------------------------------------------------------------------*/
/*! \file
\brief neighbor pair handler for discrete element method (DEM) interactions

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_DEM_NEIGHBOR_PAIRS_H
#define PARTICLE_INTERACTION_DEM_NEIGHBOR_PAIRS_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include "particle_interaction_dem_neighbor_pair_struct.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
  class ParticleContainerBundle;
}  // namespace PARTICLEENGINE

namespace PARTICLEWALL
{
  class WallHandlerInterface;
}

/*---------------------------------------------------------------------------*
 | type definitions                                           sfuchs 01/2019 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  typedef std::vector<PARTICLEINTERACTION::DEMParticlePair> DEMParticlePairData;
  typedef std::vector<PARTICLEINTERACTION::DEMParticleWallPair> DEMParticleWallPairData;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMNeighborPairs final
  {
   public:
    //! constructor
    explicit DEMNeighborPairs();

    //! destructor
    ~DEMNeighborPairs(){};

    //! init neighbor pair handler
    void Init();

    //! setup neighbor pair handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface);

    //! write restart of neighbor pair handler
    void WriteRestart(const int step, const double time) const;

    //! read restart of neighbor pair handler
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! get reference to particle pair data
    inline const DEMParticlePairData& GetRefToParticlePairData() const
    {
      return particlepairdata_;
    };

    //! get reference to particle-wall pair data
    inline const DEMParticleWallPairData& GetRefToParticleWallPairData() const
    {
      return particlewallpairdata_;
    };

    //! get reference to adhesion particle pair data
    inline const DEMParticlePairData& GetRefToParticlePairAdhesionData() const
    {
      return particlepairadhesiondata_;
    };

    //! get reference to adhesion particle-wall pair data
    inline const DEMParticleWallPairData& GetRefToParticleWallPairAdhesionData() const
    {
      return particlewallpairadhesiondata_;
    };

    //! evaluate neighbor pairs
    void EvaluateNeighborPairs();

    //! evaluate adhesion neighbor pairs
    void EvaluateNeighborPairsAdhesion(const double& adhesion_distance);

   private:
    //! evaluate particle pairs
    void EvaluateParticlePairs();

    //! evaluate particle-wall pairs
    void EvaluateParticleWallPairs();

    //! evaluate adhesion particle pairs
    void EvaluateParticlePairsAdhesion(const double& adhesion_distance);

    //! evaluate adhesion particle-wall pairs
    void EvaluateParticleWallPairsAdhesion(const double& adhesion_distance);

    //! particle pair data with evaluated quantities
    DEMParticlePairData particlepairdata_;

    //! particle-wall pair data with evaluated quantities
    DEMParticleWallPairData particlewallpairdata_;

    //! adhesion particle pair data with evaluated quantities
    DEMParticlePairData particlepairadhesiondata_;

    //! adhesion particle-wall pair data with evaluated quantities
    DEMParticleWallPairData particlewallpairadhesiondata_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! interface to particle wall handler
    std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
