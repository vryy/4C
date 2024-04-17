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
#include "baci_config.hpp"

#include "baci_lib_resulttest.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLERIGIDBODY
{
  class RigidBodyHandlerInterface;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLERIGIDBODY
{
  /*!
   * \brief rigid body result test handler
   *
   * \author Sebastian Fuchs \date 09/2020
   */
  class RigidBodyResultTest final : public DRT::ResultTest
  {
   public:
    //! constructor
    explicit RigidBodyResultTest();

    /*!
     * \brief init rigid body result test
     *
     * \author Sebastian Fuchs \date 09/2020
     */
    void Init();

    /*!
     * \brief setup rigid body result test
     *
     * \author Sebastian Fuchs \date 09/2020
     *
     * \param[in] particlerigidbodyinterface interface to rigid body handler
     */
    void Setup(const std::shared_ptr<PARTICLERIGIDBODY::RigidBodyHandlerInterface>
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
    void TestSpecial(INPUT::LineDefinition& res, int& nerr, int& test_count) override;

   private:
    //! interface to rigid body handler
    std::shared_ptr<PARTICLERIGIDBODY::RigidBodyHandlerInterface> particlerigidbodyinterface_;
  };

}  // namespace PARTICLERIGIDBODY

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
