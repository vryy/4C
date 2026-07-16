// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_PARTICLE_ENGINE_REFRESH_ENTRY_HPP
#define FOUR_C_PARTICLE_ENGINE_REFRESH_ENTRY_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_particle_engine_typedefs.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace Core::Communication
{
  class PackBuffer;
  class UnpackBuffer;
}  // namespace Core::Communication

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace Particle
{
  /*!
   * \brief entry for refreshing the states of a ghosted particle on another processor
   *
   * pack() and unpack() together define the wire format used by refresh_particles() and
   * refresh_particles_of_specific_states_and_types() to update states of ghosted particles.
   * Both must be kept in sync: unpack() extracts fields in exactly the order pack() writes
   * them.
   */
  struct ParticleRefreshEntry
  {
    //! particle type
    ParticleType type;

    //! local index of the particle in the container of ghosted particles on the target processor
    int ghostedindex;

    //! states of the particle
    ParticleStates states;

    /*!
     * \brief pack header (type, ghosted index) and append pre-packed states to a send buffer
     *
     * The states are passed in already packed since, for a single owned particle that is
     * ghosted on multiple processors, the same packed bytes are appended to every target
     * processor's send buffer without re-serializing the states for each target.
     *
     * \param[in]  type             particle type
     * \param[in]  ghostedindex     local index of the particle in the container of ghosted
     *                              particles on the target processor
     * \param[in]  prepacked_states pre-packed particle states
     * \param[out] sendbuffer       send buffer of the target processor to append to
     */
    static void pack(ParticleType type, int ghostedindex,
        const Core::Communication::PackBuffer& prepacked_states, std::vector<char>& sendbuffer);

    /*!
     * \brief unpack one refresh entry previously packed by pack()
     *
     *
     * \param[in,out] buffer buffer to unpack from
     * \return unpacked refresh entry
     */
    static ParticleRefreshEntry unpack(Core::Communication::UnpackBuffer& buffer);
  };

}  // namespace Particle

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
