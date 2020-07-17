/*---------------------------------------------------------------------------*/
/*! \file
\brief particle wall result test for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_WALL_RESULT_TEST_H
#define PARTICLE_WALL_RESULT_TEST_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_typedefs.H"

#include "../drt_lib/drt_resulttest.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLEWALL
{
  class WallHandlerInterface;
}

namespace DRT
{
  class Discretization;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEWALL
{
  /*!
   * \brief particle wall result test handler
   *
   * \author Sebastian Fuchs \date 03/2019
   */
  class WallResultTest final : public DRT::ResultTest
  {
   public:
    //! constructor
    explicit WallResultTest();

    /*!
     * \brief init wall result test
     *
     * \author Sebastian Fuchs \date 03/2019
     */
    void Init();

    /*!
     * \brief setup wall result test
     *
     * \author Sebastian Fuchs \date 03/2019
     *
     * \param[in] particleengineinterface interface to particle engine
     */
    void Setup(const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface);

    /*!
     * \brief test node value
     *
     * \author Sebastian Fuchs \date 03/2019
     *
     * \param[in]  res        result line definition
     * \param[out] nerr       number of tests with errors
     * \param[out] test_count number of tests performed
     */
    void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count) override;

    /*!
     * \brief test special quantity
     *
     * \author Sebastian Fuchs \date 03/2019
     *
     * \param[in]  res        result line definition
     * \param[out] nerr       number of tests with errors
     * \param[out] test_count number of tests performed
     */
    void TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count) override;

   private:
    //! interface to particle wall handler
    std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface_;

    //! wall discretization
    Teuchos::RCP<const DRT::Discretization> walldiscretization_;
  };

}  // namespace PARTICLEWALL

/*---------------------------------------------------------------------------*/
#endif
