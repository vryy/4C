/*---------------------------------------------------------------------------*/
/*! \file
\brief neighbor pair handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_INTERACTION_SPH_NEIGHBOR_PAIRS_HPP
#define FOUR_C_PARTICLE_INTERACTION_SPH_NEIGHBOR_PAIRS_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_particle_engine_enums.hpp"
#include "4C_particle_engine_typedefs.hpp"
#include "4C_particle_interaction_sph_neighbor_pair_struct.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
  class ParticleContainerBundle;
}  // namespace PARTICLEENGINE

namespace PARTICLEWALL
{
  class WallHandlerInterface;
}

namespace ParticleInteraction
{
  class SPHKernelBase;
}

/*---------------------------------------------------------------------------*
 | type definitions                                                          |
 *---------------------------------------------------------------------------*/
namespace ParticleInteraction
{
  using SPHParticlePairData = std::vector<ParticleInteraction::SPHParticlePair>;
  using SPHParticleWallPairData = std::vector<ParticleInteraction::SPHParticleWallPair>;
  using SPHIndexOfParticlePairs = std::vector<std::vector<std::vector<int>>>;
  using SPHIndexOfParticleWallPairs = std::vector<std::vector<int>>;
}  // namespace ParticleInteraction

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace ParticleInteraction
{
  class SPHNeighborPairs final
  {
   public:
    //! constructor
    explicit SPHNeighborPairs();

    //! init neighbor pair handler
    void Init();

    //! setup neighbor pair handler
    void setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface,
        const std::shared_ptr<ParticleInteraction::SPHKernelBase> kernel);

    //! get reference to particle pair data
    inline const SPHParticlePairData& get_ref_to_particle_pair_data() const
    {
      return particlepairdata_;
    };

    //! get reference to particle-wall pair data
    inline const SPHParticleWallPairData& get_ref_to_particle_wall_pair_data() const
    {
      return particlewallpairdata_;
    };

    //! get relevant particle pair indices for disjoint combination of particle types
    void get_relevant_particle_pair_indices_for_disjoint_combination(
        const std::set<PARTICLEENGINE::TypeEnum>& types_a,
        const std::set<PARTICLEENGINE::TypeEnum>& types_b, std::vector<int>& relindices) const;

    //! get relevant particle pair indices for equal combination of particle types
    void get_relevant_particle_pair_indices_for_equal_combination(
        const std::set<PARTICLEENGINE::TypeEnum>& types_a, std::vector<int>& relindices) const;

    //! get relevant particle wall pair indices for specific particle types
    void get_relevant_particle_wall_pair_indices(
        const std::set<PARTICLEENGINE::TypeEnum>& types_a, std::vector<int>& relindices) const;

    //! evaluate neighbor pairs
    void evaluate_neighbor_pairs();

   private:
    //! evaluate particle pairs
    void evaluate_particle_pairs();

    //! evaluate particle-wall pairs
    void evaluate_particle_wall_pairs();

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
    std::shared_ptr<ParticleInteraction::SPHKernelBase> kernel_;
  };

}  // namespace ParticleInteraction

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
