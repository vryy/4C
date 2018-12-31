/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_dem_neighbor_pairs.H

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
#include "particle_interaction_typedefs.H"

#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include "particle_interaction_dem_particle_pair.H"

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
    virtual ~DEMNeighborPairs(){};

    //! init neighbor pair handler
    virtual void Init();

    //! setup neighbor pair handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);

    //! write restart of neighbor pair handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of neighbor pair handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! get reference to map holding neighbor pairs
    const DEMNeighborPairsMap& GetRefToNeighborPairsMap() const { return neighborpairsmap_; };

    //! evaluate neighbor pairs
    virtual void EvaluateNeighborPairs();

   private:
    //! map of neighbor pairs with evaluated quantities
    PARTICLEINTERACTION::DEMNeighborPairsMap neighborpairsmap_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
