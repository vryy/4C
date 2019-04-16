/*---------------------------------------------------------------------------*/
/*!
\file particle_wall_result_test.H

\brief particle wall result test for particle simulations

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 03/2019 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_WALL_RESULT_TEST_H
#define PARTICLE_WALL_RESULT_TEST_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2019 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_typedefs.H"

#include "../drt_lib/drt_resulttest.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 03/2019 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  class WallHandlerInterface;
}

namespace DRT
{
  class Discretization;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 03/2019 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  class WallResultTest : public DRT::ResultTest
  {
   public:
    //! constructor
    explicit WallResultTest();

    //! virtual destructor
    virtual ~WallResultTest(){};

    //! init wall result test
    virtual void Init();

    //! setup wall result test
    virtual void Setup(
        const std::shared_ptr<PARTICLEALGORITHM::WallHandlerInterface> particlewallinterface);

    //! test node value
    virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count) override;

    //! test special quantity
    virtual void TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count) override;

   private:
    //! interface to particle wall handler
    std::shared_ptr<PARTICLEALGORITHM::WallHandlerInterface> particlewallinterface_;

    //! wall discretization
    Teuchos::RCP<const DRT::Discretization> walldiscretization_;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
