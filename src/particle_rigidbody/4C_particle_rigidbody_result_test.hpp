/*---------------------------------------------------------------------------*/
/*! \file
\brief rigid body result test for particle simulations
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_RIGIDBODY_RESULT_TEST_HPP
#define FOUR_C_PARTICLE_RIGIDBODY_RESULT_TEST_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_utils_result_test.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace ParticleRigidBody
{
  class RigidBodyHandlerInterface;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace ParticleRigidBody
{
  /*!
   * \brief rigid body result test handler
   *
   * \author Sebastian Fuchs \date 09/2020
   */
  class RigidBodyResultTest final : public Core::UTILS::ResultTest
  {
   public:
    //! constructor
    explicit RigidBodyResultTest();

    /*!
     * \brief init rigid body result test
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void init();

    /*!
     * \brief setup rigid body result test
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] particlerigidbodyinterface interface to rigid body handler
     */
    void setup(const std::shared_ptr<ParticleRigidBody::RigidBodyHandlerInterface>
            particlerigidbodyinterface);

    /*!
     * \brief test special quantity
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in]  res        result line definition
     * \param[out] nerr       number of tests with errors
     * \param[out] test_count number of tests performed
     */
    void test_special(Input::LineDefinition& res, int& nerr, int& test_count) override;

   private:
    //! interface to rigid body handler
    std::shared_ptr<ParticleRigidBody::RigidBodyHandlerInterface> particlerigidbodyinterface_;
  };

}  // namespace ParticleRigidBody

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
