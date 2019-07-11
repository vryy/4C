/*----------------------------------------------------------------------*/
/*!

\brief

\level 2

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_UNIT_ELASTHYPER_H
#define BACI_UNIT_ELASTHYPER_H

#include <cxxtest/TestSuite.h>
#include <Teuchos_RCPDecl.hpp>
#include <src/drt_mat/material_service.H>
#include "src/drt_mat/elasthyper.H"
#include "src/linalg/linalg_fixedsizematrix.H"

namespace MAT
{
  class ElastHyper_TestSuite;
}

// class implementation
class MAT::ElastHyper_TestSuite : public CxxTest::TestSuite
{
 public:
  /// Setup Testsuite
  void setUp() override
  {
    glstrain_ = Teuchos::rcp(new LINALG::Matrix<6, 1>);

    // create dummy Right Cauchy-Green tensor for testing of the other methods
    (*glstrain_)(0, 0) = 0.1;
    (*glstrain_)(1, 0) = 0.2;
    (*glstrain_)(2, 0) = 0.3;
    (*glstrain_)(3, 0) = 0.01;
    (*glstrain_)(4, 0) = 0.02;
    (*glstrain_)(5, 0) = 0.03;
  }

  /// Tear down test suite
  void tearDown() override {}

  void test_EvaluateKinQuant()
  {
    LINALG::Matrix<6, 1> C_strain;
    LINALG::Matrix<6, 1> iC_strain;
    LINALG::Matrix<3, 1> prinv;

    EvaluateRightCauchyGreenStrainLikeVoigt(*glstrain_, C_strain);
    InvariantsPrincipal<VoigtNotation::strain>(prinv, C_strain);
    VStrainUtils::InverseTensor(C_strain, iC_strain);

    // test right cauchy green in strain like Voigt notation
    TS_ASSERT_DELTA(C_strain(0), 1.2, 1e-10);
    TS_ASSERT_DELTA(C_strain(1), 1.4, 1e-10);
    TS_ASSERT_DELTA(C_strain(2), 1.6, 1e-10);
    TS_ASSERT_DELTA(C_strain(3), 0.02, 1e-10);
    TS_ASSERT_DELTA(C_strain(4), 0.04, 1e-10);
    TS_ASSERT_DELTA(C_strain(5), 0.06, 1e-10);

    // test inverse right cauchy green in stress like Voigt notation
    TS_ASSERT_DELTA(iC_strain(0), 0.83377015, 1e-7);
    TS_ASSERT_DELTA(iC_strain(1), 0.71445271, 1e-7);
    TS_ASSERT_DELTA(iC_strain(2), 0.62540207, 1e-7);
    TS_ASSERT_DELTA(iC_strain(3), 2 * -0.00573319, 1e-7);
    TS_ASSERT_DELTA(iC_strain(4), 2 * -0.00882316, 1e-7);
    TS_ASSERT_DELTA(iC_strain(5), 2 * -0.01556153, 1e-7);

    // test principal invariants
    TS_ASSERT_DELTA(prinv(0), 4.1999999999999993, 1e-10)
    TS_ASSERT_DELTA(prinv(1), 5.8386000000000005, 1e-10)
    TS_ASSERT_DELTA(prinv(2), 2.6861119999999996, 1e-10)
  }

 private:
  Teuchos::RCP<LINALG::Matrix<6, 1>> glstrain_;
};

#endif  // BACI_UNIT_ELASTHYPER_H
