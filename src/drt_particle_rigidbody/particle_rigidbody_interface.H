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
#include <vector>
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

    /*!
     * \brief update positions with given time increment
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] timeincrement time increment
     */
    virtual void UpdatePositions(const double timeincrement) = 0;

    /*!
     * \brief update velocities with given time increment
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] timeincrement time increment
     */
    virtual void UpdateVelocities(const double timeincrement) = 0;

    /*!
     * \brief clear accelerations
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    virtual void ClearAccelerations() = 0;

    /*!
     * \brief get rigid body data state container
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \return rigid body data state container
     */
    virtual std::shared_ptr<PARTICLERIGIDBODY::RigidBodyDataState> GetRigidBodyDataState()
        const = 0;

    /*!
     * \brief get owned rigid bodies by this processor
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \return owned rigid bodies by this processor
     */
    virtual const std::vector<int>& GetOwnedRigidBodies() const = 0;
  };

}  // namespace PARTICLERIGIDBODY

/*---------------------------------------------------------------------------*/
#endif
