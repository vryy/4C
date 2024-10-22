// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_PARTICLE_RIGIDBODY_INITIAL_FIELD_HPP
#define FOUR_C_PARTICLE_RIGIDBODY_INITIAL_FIELD_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_particle_engine_typedefs.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

FOUR_C_NAMESPACE_OPEN

namespace ParticleRigidBody
{
  /*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
  class RigidBodyDataState;

  /*!
   * \brief set initial fields
   *
   */
  void set_initial_fields(const Teuchos::ParameterList& params,
      const std::vector<int>& ownedrigidbodies,
      ParticleRigidBody::RigidBodyDataState& rigidbodydatastates);

}  // namespace ParticleRigidBody

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
