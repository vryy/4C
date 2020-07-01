/*----------------------------------------------------------------------*/
/*! \file

\brief Testcases for the material_service functions

\level 2


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_UNIT_MATERIAL_SERVICE_H
#define BACI_UNIT_MATERIAL_SERVICE_H

#include <cxxtest/TestSuite.h>
#include <src/linalg/linalg_fixedsizematrix.H>
#include <Teuchos_RCPDecl.hpp>
#include "src/drt_mat/material_service.H"

namespace MAT
{
  class MaterialService_TestSuite;
}

// class implementation
class MAT::MaterialService_TestSuite : public CxxTest::TestSuite
{
 public:
  /// Setup Testsuite
  void setUp() override
  {
    tens = Teuchos::rcp(new LINALG::Matrix<3, 3>);

    (*tens)(0, 0) = 1.1;
    (*tens)(1, 1) = 1.2;
    (*tens)(2, 2) = 1.3;
    (*tens)(0, 1) = (*tens)(1, 0) = 0.01;
    (*tens)(1, 2) = (*tens)(2, 1) = 0.02;
    (*tens)(0, 2) = (*tens)(2, 0) = 0.03;
  }

  /// Tear down test suite
  void tearDown() override {}

  void test_invariants_principal()
  {
    LINALG::Matrix<3, 1> prinv(false);
    MAT::InvariantsPrincipal(prinv, *tens);
    TS_ASSERT_DELTA(prinv(0), 3.5999999999999996, 1e-10);
    TS_ASSERT_DELTA(prinv(1), 4.3085999999999984, 1e-10);
    TS_ASSERT_DELTA(prinv(2), 1.7143620000000002, 1e-10);
  }

 private:
  Teuchos::RCP<LINALG::Matrix<3, 3>> tens;
};

#endif  // BACI_UNIT_MATERIAL_SERVICE_H
