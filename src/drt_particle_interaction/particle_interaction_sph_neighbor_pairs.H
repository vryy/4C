/*---------------------------------------------------------------------------*/
/*! \file
\brief neighbor pair handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_NEIGHBOR_PAIRS_H
#define PARTICLE_INTERACTION_SPH_NEIGHBOR_PAIRS_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include "particle_interaction_sph_neighbor_pair_struct.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
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

namespace PARTICLEINTERACTION
{
  class SPHKernelBase;
}

/*---------------------------------------------------------------------------*
 | type definitions                                           sfuchs 01/2019 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  typedef std::vector<PARTICLEINTERACTION::SPHParticlePair> SPHParticlePairData;
  typedef std::vector<PARTICLEINTERACTION::SPHParticleWallPair> SPHParticleWallPairData;
  typedef std::vector<std::vector<std::vector<int>>> SPHIndexOfParticlePairs;
  typedef std::vector<std::vector<int>> SPHIndexOfParticleWallPairs;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHNeighborPairs final
  {
   public:
    //! constructor
    explicit SPHNeighborPairs();

    //! destructor
    ~SPHNeighborPairs() = default;

    //! init neighbor pair handler
    void Init();

    //! setup neighbor pair handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel);

    //! get reference to particle pair data
    inline const SPHParticlePairData& GetRefToParticlePairData() const
    {
      return particlepairdata_;
    };

    //! get reference to particle-wall pair data
    inline const SPHParticleWallPairData& GetRefToParticleWallPairData() const
    {
      return particlewallpairdata_;
    };

    //! get relevant particle pair indices for disjoint combination of particle types
    void GetRelevantParticlePairIndicesForDisjointCombination(
        const std::set<PARTICLEENGINE::TypeEnum>& types_a,
        const std::set<PARTICLEENGINE::TypeEnum>& types_b, std::vector<int>& relindices) const;

    //! get relevant particle pair indices for equal combination of particle types
    void GetRelevantParticlePairIndicesForEqualCombination(
        const std::set<PARTICLEENGINE::TypeEnum>& types_a, std::vector<int>& relindices) const;

    //! get relevant particle wall pair indices for specific particle types
    void GetRelevantParticleWallPairIndices(
        const std::set<PARTICLEENGINE::TypeEnum>& types_a, std::vector<int>& relindices) const;

    //! evaluate neighbor pairs
    void EvaluateNeighborPairs();

   private:
    //! evaluate particle pairs
    void EvaluateParticlePairs();

    //! evaluate particle-wall pairs
    void EvaluateParticleWallPairs();

    //! particle pair data with evaluated quantities
    SPHParticlePairData particlepairdata_;

    //! particle-wall pair data with evaluated quantities
    SPHParticleWallPairData particlewallpairdata_;

    //! index of particle pairs for each type
    SPHIndexOfParticlePairs indexofparticlepairs_;

    //! index of particle-wall pairs for each type
    SPHIndexOfParticleWallPairs indexofparticlewallpairs_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! interface to particle wall handler
    std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface_;

    //! kernel handler
    std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
