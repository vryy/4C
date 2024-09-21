/*---------------------------------------------------------------------------*/
/*! \file
\brief initial field handler for rigid bodies
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
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
