/*---------------------------------------------------------------------------*/
/*! \file
\brief particle result test for particle simulations
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_ALGORITHM_RESULT_TEST_HPP
#define FOUR_C_PARTICLE_ALGORITHM_RESULT_TEST_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_lib_resulttest.hpp"
#include "baci_particle_engine_typedefs.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  /*!
   * \brief particle field result test handler
   *
   * \author Sebastian Fuchs \date 07/2018
   */
  class ParticleResultTest final : public DRT::ResultTest
  {
   public:
    //! constructor
    explicit ParticleResultTest();

    /*!
     * \brief init particle result test
     *
     * \author Sebastian Fuchs \date 07/2018
     */
    void Init();

    /*!
     * \brief setup particle result test
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \param[in] particleengineinterface interface to particle engine
     */
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);

    /*!
     * \brief test special quantity
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \param[in]  res        result line definition
     * \param[out] nerr       number of tests with errors
     * \param[out] test_count number of tests performed
     */
    void TestSpecial(INPUT::LineDefinition& res, int& nerr, int& test_count) override;

   private:
    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
