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

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace PARTICLERIGIDBODY
{
  /*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
  class RigidBodyDataState;

  /*!
   * \brief set initial fields
   *
   */
  void SetInitialFields(const Teuchos::ParameterList& params,
      const std::vector<int>& ownedrigidbodies,
      PARTICLERIGIDBODY::RigidBodyDataState& rigidbodydatastates);

}  // namespace PARTICLERIGIDBODY

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
