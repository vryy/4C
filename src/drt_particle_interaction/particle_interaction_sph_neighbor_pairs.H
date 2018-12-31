/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_neighbor_pairs.H

\brief neighbor pair handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

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
#include "particle_interaction_typedefs.H"

#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include "particle_interaction_sph_particle_pair.H"

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
    virtual ~SPHNeighborPairs(){};

    //! init neighbor pair handler
    virtual void Init();

    //! setup neighbor pair handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel);

    //! write restart of neighbor pair handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of neighbor pair handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! get reference to map holding neighbor pairs
    const SPHNeighborPairsMap& GetRefToNeighborPairsMap() const { return neighborpairsmap_; };

    //! evaluate neighbor pairs
    virtual void EvaluateNeighborPairs();

   private:
    //! map of neighbor pairs with evaluated absolute distance, versor and kernel
    PARTICLEINTERACTION::SPHNeighborPairsMap neighborpairsmap_;

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
