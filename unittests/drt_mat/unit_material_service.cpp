/*----------------------------------------------------------------------*/
/*! \file
\brief Testcases for the material_service functions
\level 2

*/
/*----------------------------------------------------------------------*/

#include "gtest/gtest.h"

#include "src/drt_mat/material_service.H"

#include "src/linalg/linalg_fixedsizematrix.H"

#include "unittests/common/assertions.h"

namespace
{
  TEST(MaterialServiceTest, TestInvariantsPrincipal)
  {
    LINALG::Matrix<3, 3> sym_tensor(false);
    sym_tensor(0, 0) = 1.1;
    sym_tensor(1, 1) = 1.2;
    sym_tensor(2, 2) = 1.3;
    sym_tensor(0, 1) = sym_tensor(1, 0) = 0.01;
    sym_tensor(1, 2) = sym_tensor(2, 1) = 0.02;
    sym_tensor(0, 2) = sym_tensor(2, 0) = 0.03;

    LINALG::Matrix<3, 1> prinv(false);
    MAT::InvariantsPrincipal(prinv, sym_tensor);

    LINALG::Matrix<3, 1> prinv_reference(false);
    prinv_reference(0) = 3.5999999999999996;
    prinv_reference(1) = 4.3085999999999984;
    prinv_reference(2) = 1.7143620000000002;

    BACI_EXPECT_NEAR(prinv, prinv_reference, 1.0e-10);
  }
}  // namespace