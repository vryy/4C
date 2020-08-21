/*---------------------------------------------------------------------------*/
/*! \file
\brief interface to provide restricted access to rigid body handler
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_RIGIDBODY_INTERFACE_H
#define PARTICLE_RIGIDBODY_INTERFACE_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include <memory>

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLERIGIDBODY
{
  class RigidBodyDataState;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLERIGIDBODY
{
  /*!
   * \brief interface to provide restricted access to rigid body handler
   *
   * \author Sebastian Fuchs \date 09/2020
   */
  class RigidBodyHandlerInterface
  {
   public:
    //! virtual destructor
    virtual ~RigidBodyHandlerInterface() = default;
  };

}  // namespace PARTICLERIGIDBODY

/*---------------------------------------------------------------------------*/
#endif
