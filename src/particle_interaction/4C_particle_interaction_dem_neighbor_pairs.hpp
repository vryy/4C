/*---------------------------------------------------------------------------*/
/*! \file
\brief neighbor pair handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_INTERACTION_DEM_NEIGHBOR_PAIRS_HPP
#define FOUR_C_PARTICLE_INTERACTION_DEM_NEIGHBOR_PAIRS_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_particle_engine_enums.hpp"
#include "4C_particle_engine_typedefs.hpp"
#include "4C_particle_interaction_dem_neighbor_pair_struct.hpp"

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

/*---------------------------------------------------------------------------*
 | type definitions                                                          |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  using DEMParticlePairData = std::vector<PARTICLEINTERACTION::DEMParticlePair>;
  using DEMParticleWallPairData = std::vector<PARTICLEINTERACTION::DEMParticleWallPair>;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMNeighborPairs final
  {
   public:
    //! constructor
    explicit DEMNeighborPairs();

    //! init neighbor pair handler
    void Init();

    //! setup neighbor pair handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface);

    //! get reference to particle pair data
    inline const DEMParticlePairData& get_ref_to_particle_pair_data() const
    {
      return particlepairdata_;
    };

    //! get reference to particle-wall pair data
    inline const DEMParticleWallPairData& get_ref_to_particle_wall_pair_data() const
    {
      return particlewallpairdata_;
    };

    //! get reference to adhesion particle pair data
    inline const DEMParticlePairData& get_ref_to_particle_pair_adhesion_data() const
    {
      return particlepairadhesiondata_;
    };

    //! get reference to adhesion particle-wall pair data
    inline const DEMParticleWallPairData& get_ref_to_particle_wall_pair_adhesion_data() const
    {
      return particlewallpairadhesiondata_;
    };

    //! evaluate neighbor pairs
    void evaluate_neighbor_pairs();

    //! evaluate adhesion neighbor pairs
    void evaluate_neighbor_pairs_adhesion(const double& adhesion_distance);

   private:
    //! evaluate particle pairs
    void evaluate_particle_pairs();

    //! evaluate particle-wall pairs
    void evaluate_particle_wall_pairs();

    //! evaluate adhesion particle pairs
    void evaluate_particle_pairs_adhesion(const double& adhesion_distance);

    //! evaluate adhesion particle-wall pairs
    void evaluate_particle_wall_pairs_adhesion(const double& adhesion_distance);

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
FOUR_C_NAMESPACE_CLOSE

#endif
