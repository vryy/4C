/*---------------------------------------------------------------------------*/
/*!
\file particle_engine_interface.H

\brief pure virtual interface to provide access to particle engine class

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_ENGINE_INTERFACE_H
#define PARTICLE_ENGINE_INTERFACE_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_enums.H"
#include "particle_container_bundle.H"

#include "../linalg/linalg_fixedsizematrix.H"

#include <unordered_map>

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleEngineInterface
  {
   public:
    //! virtual destructor
    virtual ~ParticleEngineInterface(){};

    //! refresh particles of specific states and types
    virtual void RefreshParticlesOfSpecificStatesAndTypes(
        const StatesOfTypesToRefresh& particlestatestotypes) const = 0;

    //! change type of particles
    virtual void TypeChangeParticles(std::vector<std::set<int>>& particlestoremove,
        std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoinsert) = 0;

    //! get particle container bundle
    virtual const ParticleContainerBundleShrdPtr GetParticleContainerBundle() const = 0;

    //! get reference to potential particle neighbors
    virtual const PotentialParticleNeighbors& GetPotentialParticleNeighbors() const = 0;

    //! get reference to particles being communicated to target processors
    virtual const std::vector<std::vector<int>>& GetCommunicatedParticleTargets() const = 0;

    //! get local index in specific particle container
    virtual const LocalIndexTupleShrdPtr GetLocalIndexInSpecificContainer(int globalid) const = 0;

    //! relate all particles to all processors
    virtual void RelateAllParticlesToAllProcs(std::vector<int>& particlestoproc) const = 0;

    //! return bin size
    virtual const double* BinSize() const = 0;

    //! return minimum relevant bin size
    virtual double MinBinSize() const = 0;

    //! return flag whether pbcs are applied to the problem
    virtual bool HavePBC(const int dim = -1) const = 0;

    //! return delta for pbc in x, y, or z direction
    virtual double PBCDelta(const int dim) const = 0;

    //! get bounding box dimensions
    virtual LINALG::Matrix<3, 2>& XAABB() const = 0;

    //! get distance between particles considering periodic boundaries
    virtual void DistanceBetweenParticles(
        const double* pos_i, const double* pos_j, double* r_ji) const = 0;

    //! get number of particles on this processors
    virtual int GetNumberOfParticles() const = 0;

    //! get number of particles on this processor of specific type
    virtual int GetNumberOfParticlesOfSpecificType(const TypeEnum particleType) const = 0;
  };

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
