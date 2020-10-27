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
  using TypeEnum = PARTICLEENGINE::ParticleType;

  //! particle status enum
  using StatusEnum = PARTICLEENGINE::ParticleStatus;

  //! particle state enum
  using StateEnum = PARTICLEENGINE::ParticleState;

  //! states of particle indexed by particle state enum
  using ParticleStates = std::vector<std::vector<double>>;

  //! shared pointer to particle object
  using ParticleObjShrdPtr = std::shared_ptr<PARTICLEENGINE::ParticleObject>;

  //! collection of particle containers indexed by particle type enum and particle status enum
  using TypeStatusContainers =
      std::vector<std::vector<std::shared_ptr<PARTICLEENGINE::ParticleContainer>>>;

  //! shared pointer to particle container bundle
  using ParticleContainerBundleShrdPtr = std::shared_ptr<PARTICLEENGINE::ParticleContainerBundle>;

  /*!
   * \brief local index tuple of a particle
   *
   * Local index tuple of a particle consisting of particle type enum, particle status enum, and
   * index of particle in container.
   */
  using LocalIndexTuple = std::tuple<TypeEnum, StatusEnum, int>;

  //! shared pointer to local index tuple of a particle
  using LocalIndexTupleShrdPtr = std::shared_ptr<LocalIndexTuple>;

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
  using ParticlesToBins = std::vector<std::vector<std::pair<TypeEnum, int>>>;

  //! potential particle neighbor pairs
  using PotentialParticleNeighbors = std::vector<std::pair<LocalIndexTuple, LocalIndexTuple>>;

  //! relate particle type enums and particle state enums to be refreshed
  using StatesOfTypesToRefresh = std::vector<std::pair<TypeEnum, std::vector<StateEnum>>>;

  /*!
   * \brief relate range of owned bins to corresponding column wall element
   *
   * Relate the range of owned bins (global id of bin) to the column wall element occupying that
   * bins. The column wall elements are indexed by the local id of their column map.
   */
  using BinsToColWallEles = std::vector<std::vector<int>>;

  //! potential particle wall neighbor pairs
  using PotentialWallNeighbors = std::vector<std::pair<LocalIndexTuple, DRT::Element*>>;

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
