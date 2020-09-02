/*---------------------------------------------------------------------------*/
/*! \file
\brief particle type definitions
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_TYPEDEFS_H
#define PARTICLE_TYPEDEFS_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_enums.H"

#include <map>
#include <set>
#include <tuple>
#include <vector>
#include <memory>

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleContainer;
  class ParticleContainerBundle;
  class ParticleObject;
}  // namespace PARTICLEENGINE

namespace DRT
{
  class Element;
}

/*---------------------------------------------------------------------------*
 | type definitions                                                          |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  //! particle type enum
  typedef enum PARTICLEENGINE::ParticleType TypeEnum;

  //! particle status enum
  typedef enum PARTICLEENGINE::ParticleStatus StatusEnum;

  //! particle state enum
  typedef enum PARTICLEENGINE::ParticleState StateEnum;

  //! states of particle indexed by particle state enum
  typedef std::vector<std::vector<double>> ParticleStates;

  //! shared pointer to particle object
  typedef std::shared_ptr<PARTICLEENGINE::ParticleObject> ParticleObjShrdPtr;

  //! collection of particle containers indexed by particle type enum and particle status enum
  typedef std::vector<std::vector<std::shared_ptr<PARTICLEENGINE::ParticleContainer>>>
      TypeStatusContainers;

  //! shared pointer to particle container bundle
  typedef std::shared_ptr<PARTICLEENGINE::ParticleContainerBundle> ParticleContainerBundleShrdPtr;

  /*!
   * \brief local index tuple of a particle
   *
   * Local index tuple of a particle consisting of particle type enum, particle status enum, and
   * index of particle in container.
   */
  typedef std::tuple<TypeEnum, StatusEnum, int> LocalIndexTuple;

  //! shared pointer to local index tuple of a particle
  typedef std::shared_ptr<LocalIndexTuple> LocalIndexTupleShrdPtr;

  /*!
   * \brief relate particles to corresponding bins being located in
   *
   * Relate all particles (consisting of particle type enum and index of particle in container) to
   * their respective corresponding bin the particles are located in. The bins are indexed by the
   * local id of the bin column map.
   *
   * \note The information of the particle status (owned or ghosted) can be deduced directly from
   *       the status of the bin.
   */
  typedef std::vector<std::vector<std::pair<TypeEnum, int>>> ParticlesToBins;

  //! potential particle neighbor pairs
  typedef std::vector<std::pair<LocalIndexTuple, LocalIndexTuple>> PotentialParticleNeighbors;

  //! relate particle type enums and particle state enums to be refreshed
  typedef std::vector<std::pair<TypeEnum, std::vector<StateEnum>>> StatesOfTypesToRefresh;

  /*!
   * \brief relate range of owned bins to corresponding column wall element
   *
   * Relate the range of owned bins (global id of bin) to the column wall element occupying that
   * bins. The column wall elements are indexed by the local id of their column map.
   */
  typedef std::vector<std::vector<int>> BinsToColWallEles;

  //! potential particle wall neighbor pairs
  typedef std::vector<std::pair<LocalIndexTuple, DRT::Element*>> PotentialWallNeighbors;

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
