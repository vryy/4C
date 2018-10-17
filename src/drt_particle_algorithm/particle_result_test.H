/*---------------------------------------------------------------------------*/
/*!
\file particle_result_test.H

\brief particle result test for particle simulations

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_RESULT_TEST_H
#define PARTICLE_RESULT_TEST_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_typedefs.H"

#include "../drt_lib/drt_resulttest.H"

#include <Epetra_Comm.h>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  class ResultTest : public DRT::ResultTest
  {
   public:
    //! constructor
    explicit ResultTest(const Epetra_Comm& comm);

    //! virtual destructor
    virtual ~ResultTest(){};

    //! init result test
    virtual void Init();

    //! setup result test
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);

    //! test special quantity
    void TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count);

   protected:
    //! communication
    const Epetra_Comm& comm_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
