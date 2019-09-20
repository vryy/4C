/*---------------------------------------------------------------------------*/
/*! \file
\brief particle result test for particle simulations

\level 1

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_RESULT_TEST_H
#define PARTICLE_RESULT_TEST_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_typedefs.H"

#include "../drt_lib/drt_resulttest.H"

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

    //! destructor
    ~ParticleResultTest() override{};

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
    void TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count) override;

   private:
    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
