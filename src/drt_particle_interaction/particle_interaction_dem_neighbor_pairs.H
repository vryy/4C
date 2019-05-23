/*---------------------------------------------------------------------------*/
/*!

\brief neighbor pair handler for discrete element method (DEM) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

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

namespace PARTICLEALGORITHM
{
  class WallHandlerInterface;
}

namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
  class ParticleContainerBundle;
}  // namespace PARTICLEENGINE

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
  class DEMNeighborPairs
  {
   public:
    //! constructor
    explicit DEMNeighborPairs();

    //! virtual destructor
    ~DEMNeighborPairs(){};

    //! init neighbor pair handler
    void Init();

    //! setup neighbor pair handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEALGORITHM::WallHandlerInterface> particlewallinterface);

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

    //! evaluate neighbor pairs
    void EvaluateNeighborPairs();

   private:
    //! evaluate particle pairs
    void EvaluateParticlePairs();

    //! evaluate particle-wall pairs
    void EvaluateParticleWallPairs();

    //! particle pair data with evaluated quantities
    DEMParticlePairData particlepairdata_;

    //! particle-wall pair data with evaluated quantities
    DEMParticleWallPairData particlewallpairdata_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! interface to particle wall handler
    std::shared_ptr<PARTICLEALGORITHM::WallHandlerInterface> particlewallinterface_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
