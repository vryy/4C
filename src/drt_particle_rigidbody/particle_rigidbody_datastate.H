/*---------------------------------------------------------------------------*/
/*! \file
\brief rigid body data state container
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_RIGIDBODY_DATASTATE_H
#define PARTICLE_RIGIDBODY_DATASTATE_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include <vector>
#include <memory>

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLERIGIDBODY
{
  /*!
   * \brief rigid body data state container
   *
   * \author Sebastian Fuchs \date 08/2020
   */
  class RigidBodyDataState final
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    explicit RigidBodyDataState() = default;

    //! destructor
    ~RigidBodyDataState() = default;

    /*!
     * \brief init rigid body data state container
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void Init();

    /*!
     * \brief setup rigid body data state container
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void Setup();

    /*!
     * \brief allocate stored states
     *
     * \author Sebastian Fuchs \date 08/2020
     *
     * \param[in] numrigidbodies number of rigid bodies
     */
    void AllocateStoredStates(const int numrigidbodies);

    //! @name get states (read only access)
    //! @{

    //! get mass of rigid bodies
    inline const std::vector<double>& GetRefMass() const { return mass_; };

    //! get mass moment of inertia of rigid bodies
    inline const std::vector<std::vector<double>>& GetRefInertia() const { return inertia_; };

    //! get position of rigid bodies
    inline const std::vector<std::vector<double>>& GetRefPosition() const { return position_; };

    //! get velocity of rigid bodies
    inline const std::vector<std::vector<double>>& GetRefVelocity() const { return velocity_; };

    //! get acceleration of rigid bodies
    inline const std::vector<std::vector<double>>& GetRefAcceleration() const
    {
      return acceleration_;
    };

    //! get rotation of rigid bodies
    inline const std::vector<std::vector<double>>& GetRefRotation() const { return rotation_; };

    //! get angular velocity of rigid bodies
    inline const std::vector<std::vector<double>>& GetRefAngularVelocity() const
    {
      return angularvelocity_;
    };

    //! get angular acceleration of rigid bodies
    inline const std::vector<std::vector<double>>& GetRefAngularAcceleration() const
    {
      return angularacceleration_;
    };

    //! get force of rigid bodies
    inline const std::vector<std::vector<double>>& GetRefForce() const { return force_; };

    //! get torque of rigid bodies
    inline const std::vector<std::vector<double>>& GetRefTorque() const { return torque_; };

    //! @}

    //! @name get states (read and write access)
    //! @{

    //! get mass of rigid bodies
    inline std::vector<double>& GetRefMutableMass() { return mass_; };

    //! get mass moment of inertia of rigid bodies
    inline std::vector<std::vector<double>>& GetRefMutableInertia() { return inertia_; };

    //! get position of rigid bodies
    inline std::vector<std::vector<double>>& GetRefMutablePosition() { return position_; };

    //! get velocity of rigid bodies
    inline std::vector<std::vector<double>>& GetRefMutableVelocity() { return velocity_; };

    //! get acceleration of rigid bodies
    inline std::vector<std::vector<double>>& GetRefMutableAcceleration() { return acceleration_; };

    //! get rotation of rigid bodies
    inline std::vector<std::vector<double>>& GetRefMutableRotation() { return rotation_; };

    //! get angular velocity of rigid bodies
    inline std::vector<std::vector<double>>& GetRefMutableAngularVelocity()
    {
      return angularvelocity_;
    };

    //! get angular acceleration of rigid bodies
    inline std::vector<std::vector<double>>& GetRefMutableAngularAcceleration()
    {
      return angularacceleration_;
    };

    //! get force of rigid bodies
    inline std::vector<std::vector<double>>& GetRefMutableForce() { return force_; };

    //! get torque of rigid bodies
    inline std::vector<std::vector<double>>& GetRefMutableTorque() { return torque_; };

    //! @}

   private:
    //! @name stored states
    //! @{

    //! mass of rigid bodies
    std::vector<double> mass_;

    //! mass moment of inertia of rigid bodies
    std::vector<std::vector<double>> inertia_;

    //! position of rigid bodies
    std::vector<std::vector<double>> position_;

    //! velocity of rigid bodies
    std::vector<std::vector<double>> velocity_;

    //! acceleration of rigid bodies
    std::vector<std::vector<double>> acceleration_;

    //! rotation of rigid bodies
    std::vector<std::vector<double>> rotation_;

    //! angular velocity of rigid bodies
    std::vector<std::vector<double>> angularvelocity_;

    //! angular acceleration of rigid bodies
    std::vector<std::vector<double>> angularacceleration_;

    //! force of rigid bodies
    std::vector<std::vector<double>> force_;

    //! torque of rigid bodies
    std::vector<std::vector<double>> torque_;

    //! @}
  };
}  // namespace PARTICLERIGIDBODY

/*---------------------------------------------------------------------------*/
#endif
