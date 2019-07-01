/*---------------------------------------------------------------------------*/
/*!
\brief neighbor pair handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_NEIGHBOR_PAIRS_H
#define PARTICLE_INTERACTION_SPH_NEIGHBOR_PAIRS_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include "particle_interaction_sph_neighbor_pair_struct.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 08/2018 |
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
  typedef std::vector<std::vector<int>> SPHIndexOfParticlePairs;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHNeighborPairs
  {
   public:
    //! constructor
    explicit SPHNeighborPairs();

    //! virtual destructor
    ~SPHNeighborPairs(){};

    //! init neighbor pair handler
    void Init();

    //! setup neighbor pair handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel);

    //! write restart of neighbor pair handler
    void WriteRestart(const int step, const double time) const;

    //! read restart of neighbor pair handler
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! get reference to particle pair data
    inline const SPHParticlePairData& GetRefToParticlePairData() const
    {
      return particlepairdata_;
    };

    //! get relevant particle pair indices for particle types
    void GetRelevantParticlePairIndices(
        const std::set<PARTICLEENGINE::TypeEnum>& reltypes, std::vector<int>& relindices) const;

    //! evaluate neighbor pairs
    void EvaluateNeighborPairs();

   private:
    //! evaluate particle pairs
    void EvaluateParticlePairs();

    //! particle pair data with evaluated quantities
    SPHParticlePairData particlepairdata_;

    //! index of particle pairs for each type
    SPHIndexOfParticlePairs indexofparticlepairs_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! kernel handler
    std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
