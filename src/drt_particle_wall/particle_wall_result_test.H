/*---------------------------------------------------------------------------*/
/*! \file
\brief particle wall result test for particle simulations

\level 3

\maintainer  Sebastian Fuchs
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
namespace PARTICLEWALL
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
namespace PARTICLEWALL
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
        const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface);

    //! test node value
    virtual void TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count) override;

    //! test special quantity
    virtual void TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count) override;

   private:
    //! interface to particle wall handler
    std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface_;

    //! wall discretization
    Teuchos::RCP<const DRT::Discretization> walldiscretization_;
  };

}  // namespace PARTICLEWALL

/*---------------------------------------------------------------------------*/
#endif
