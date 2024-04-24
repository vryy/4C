/*----------------------------------------------------------------------*/
/*! \file
\brief Testcases for the material_service functions
\level 2

*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_mat_service.hpp"
#include "baci_unittest_utils_assertions_test.hpp"

namespace
{
  using namespace FourC;

  TEST(MaterialServiceTest, TestInvariantsPrincipal)
  {
    CORE::LINALG::Matrix<3, 3> sym_tensor(false);
    sym_tensor(0, 0) = 1.1;
    sym_tensor(1, 1) = 1.2;
    sym_tensor(2, 2) = 1.3;
    sym_tensor(0, 1) = sym_tensor(1, 0) = 0.01;
    sym_tensor(1, 2) = sym_tensor(2, 1) = 0.02;
    sym_tensor(0, 2) = sym_tensor(2, 0) = 0.03;

    CORE::LINALG::Matrix<3, 1> prinv(false);
    MAT::InvariantsPrincipal(prinv, sym_tensor);

    CORE::LINALG::Matrix<3, 1> prinv_reference(false);
    prinv_reference(0) = 3.5999999999999996;
    prinv_reference(1) = 4.3085999999999984;
    prinv_reference(2) = 1.7143620000000002;

    FOUR_C_EXPECT_NEAR(prinv, prinv_reference, 1.0e-10);
  }

  TEST(MaterialServiceTest, TestAddDerivInvABInvBProduct)
  {
    CORE::LINALG::Matrix<6, 1> A(false);
    A(0) = 0.5;
    A(1) = 0.3;
    A(2) = 0.6;
    A(3) = 0.2;
    A(4) = 0.4;
    A(5) = 0.9;

    CORE::LINALG::Matrix<6, 1> InvABInvB(false);
    InvABInvB(0) = 1.72;
    InvABInvB(1) = 1.65;
    InvABInvB(2) = 1.13;
    InvABInvB(3) = 1.27;
    InvABInvB(4) = 1.46;
    InvABInvB(5) = 1.23;

    CORE::LINALG::Matrix<6, 6> Result(true);
    double scalar = 0.5;

    // result_ijkl = A_ik InvABInvB_jl +  A_il InvABInvB_jk + A_jk InvABInvB_il + A_jl InvABInvB_ik
    MAT::AddDerivInvABInvBProduct(scalar, A, InvABInvB, Result);

    CORE::LINALG::Matrix<6, 6> Result_reference(false);
    Result_reference(0, 0) = -0.86;
    Result_reference(0, 1) = -0.254;
    Result_reference(0, 2) = -1.107;
    Result_reference(0, 3) = -0.4895;
    Result_reference(0, 4) = -0.6945;
    Result_reference(0, 5) = -1.0815;
    Result_reference(1, 0) = -0.254;
    Result_reference(1, 1) = -0.495;
    Result_reference(1, 2) = -0.584;
    Result_reference(1, 3) = -0.3555;
    Result_reference(1, 4) = -0.549;
    Result_reference(1, 5) = -0.40;
    Result_reference(2, 0) = -1.107;
    Result_reference(2, 1) = -0.584;
    Result_reference(2, 2) = -0.678;
    Result_reference(2, 3) = -0.903;
    Result_reference(2, 4) = -0.664;
    Result_reference(2, 5) = -0.8775;
    Result_reference(3, 0) = -0.4895;
    Result_reference(3, 1) = -0.3555;
    Result_reference(3, 2) = -0.903;
    Result_reference(3, 3) = -0.46225;
    Result_reference(3, 4) = -0.6635;
    Result_reference(3, 5) = -0.70175;
    Result_reference(4, 0) = -0.6945;
    Result_reference(4, 1) = -0.549;
    Result_reference(4, 2) = -0.664;
    Result_reference(4, 3) = -0.6635;
    Result_reference(4, 4) = -0.62425;
    Result_reference(4, 5) = -0.6985;
    Result_reference(5, 0) = -1.0815;
    Result_reference(5, 1) = -0.40;
    Result_reference(5, 2) = -0.8775;
    Result_reference(5, 3) = -0.70175;
    Result_reference(5, 4) = -0.6985;
    Result_reference(5, 5) = -0.95275;

    FOUR_C_EXPECT_NEAR(Result, Result_reference, 1.0e-10);
  }

}  // namespace