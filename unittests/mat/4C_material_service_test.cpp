/*----------------------------------------------------------------------*/
/*! \file
\brief Testcases for the Material_service functions
\level 2

*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_four_tensor.hpp"
#include "4C_mat_service.hpp"
#include "4C_unittest_utils_assertions_test.hpp"

namespace
{
  using namespace FourC;

  TEST(MaterialServiceTest, TestInvariantsPrincipal)
  {
    Core::LinAlg::Matrix<3, 3> sym_tensor(false);
    sym_tensor(0, 0) = 1.1;
    sym_tensor(1, 1) = 1.2;
    sym_tensor(2, 2) = 1.3;
    sym_tensor(0, 1) = sym_tensor(1, 0) = 0.01;
    sym_tensor(1, 2) = sym_tensor(2, 1) = 0.02;
    sym_tensor(0, 2) = sym_tensor(2, 0) = 0.03;

    Core::LinAlg::Matrix<3, 1> prinv(false);
    Mat::invariants_principal(prinv, sym_tensor);

    Core::LinAlg::Matrix<3, 1> prinv_reference(false);
    prinv_reference(0) = 3.5999999999999996;
    prinv_reference(1) = 4.3085999999999984;
    prinv_reference(2) = 1.7143620000000002;

    FOUR_C_EXPECT_NEAR(prinv, prinv_reference, 1.0e-10);
  }

  TEST(MaterialServiceTest, Testadd_derivative_of_inva_b_inva_product)
  {
    Core::LinAlg::Matrix<6, 1> A(false);
    A(0) = 0.5;
    A(1) = 0.3;
    A(2) = 0.6;
    A(3) = 0.2;
    A(4) = 0.4;
    A(5) = 0.9;

    Core::LinAlg::Matrix<6, 1> InvABInvB(false);
    InvABInvB(0) = 1.72;
    InvABInvB(1) = 1.65;
    InvABInvB(2) = 1.13;
    InvABInvB(3) = 1.27;
    InvABInvB(4) = 1.46;
    InvABInvB(5) = 1.23;

    Core::LinAlg::Matrix<6, 6> Result(true);
    double scalar = 0.5;

    // result_ijkl = A_ik InvABInvB_jl +  A_il InvABInvB_jk + A_jk InvABInvB_il + A_jl InvABInvB_ik
    Mat::add_derivative_of_inva_b_inva_product(scalar, A, InvABInvB, Result);

    Core::LinAlg::Matrix<6, 6> Result_reference(false);
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

  TEST(MaterialServiceTest, TestComputeJ2)
  {
    // test the calculation of J2 invariant
    Core::LinAlg::Matrix<3, 3> stress(false);
    stress(0, 0) = 1.0;
    stress(1, 1) = 2.0;
    stress(2, 2) = 3.0;
    stress(0, 1) = stress(1, 0) = 0.1;
    stress(1, 2) = stress(2, 1) = 0.0;
    stress(0, 2) = stress(2, 0) = 0.0;
    const double j2 = Mat::second_invariant_of_deviatoric_stress(stress);
    const double j2_ref = 0.5 * (2 + 2 * 0.1 * 0.1);
    EXPECT_NEAR(j2, j2_ref, 1.0e-10);
  }

  TEST(MaterialServiceTest, TestElasticTensor)
  {
    // test the calculation of fourth order elastic tensor and subroutine to transform 4th order
    // tensor to Matrix
    Core::LinAlg::FourTensor<3> Ce;
    const double E = 2.0;
    const double NU = 0.3;
    Mat::setup_linear_isotropic_elastic_tensor(Ce, E, NU);

    Core::LinAlg::Matrix<6, 6> De;
    Mat::four_tensor_to_matrix(Ce, De);

    Core::LinAlg::Matrix<6, 6> De_ref;
    const double c1 = E / ((1.00 + NU) * (1 - 2 * NU));
    const double c2 = c1 * (1 - NU);
    const double c3 = c1 * NU;
    const double c4 = c1 * 0.5 * (1 - 2 * NU);

    De_ref(0, 0) = c2;
    De_ref(0, 1) = c3;
    De_ref(0, 2) = c3;
    De_ref(0, 3) = 0.;
    De_ref(0, 4) = 0.;
    De_ref(0, 5) = 0.;
    De_ref(1, 0) = c3;
    De_ref(1, 1) = c2;
    De_ref(1, 2) = c3;
    De_ref(1, 3) = 0.;
    De_ref(1, 4) = 0.;
    De_ref(1, 5) = 0.;
    De_ref(2, 0) = c3;
    De_ref(2, 1) = c3;
    De_ref(2, 2) = c2;
    De_ref(2, 3) = 0.;
    De_ref(2, 4) = 0.;
    De_ref(2, 5) = 0.;
    De_ref(3, 0) = 0.;
    De_ref(3, 1) = 0.;
    De_ref(3, 2) = 0.;
    De_ref(3, 3) = c4;
    De_ref(3, 4) = 0.;
    De_ref(3, 5) = 0.;
    De_ref(4, 0) = 0.;
    De_ref(4, 1) = 0.;
    De_ref(4, 2) = 0.;
    De_ref(4, 3) = 0.;
    De_ref(4, 4) = c4;
    De_ref(4, 5) = 0.;
    De_ref(5, 0) = 0.;
    De_ref(5, 1) = 0.;
    De_ref(5, 2) = 0.;
    De_ref(5, 3) = 0.;
    De_ref(5, 4) = 0.;
    De_ref(5, 5) = c4;

    FOUR_C_EXPECT_NEAR(De, De_ref, 1.0e-10);
  }

  TEST(MaterialServiceTest, TestDeviatoricTensor)
  {
    // test the calculation of fourth order deviatoric tensor and contraction of fourth order tensor
    // and Matrix
    Core::LinAlg::FourTensor<3> Id;
    Mat::setup_deviatoric_projection_tensor(Id);

    Core::LinAlg::Matrix<3, 3> stress(false);
    stress(0, 0) = 1.0;
    stress(1, 1) = 2.0;
    stress(2, 2) = 3.0;
    stress(0, 1) = stress(1, 0) = 0.4;
    stress(1, 2) = stress(2, 1) = 0.5;
    stress(0, 2) = stress(2, 0) = 0.6;

    Core::LinAlg::Matrix<3, 3> s_ref(false);
    s_ref(0, 0) = -1.0;
    s_ref(1, 1) = 0.0;
    s_ref(2, 2) = 1.0;
    s_ref(0, 1) = s_ref(1, 0) = 0.4;
    s_ref(1, 2) = s_ref(2, 1) = 0.5;
    s_ref(0, 2) = s_ref(2, 0) = 0.6;

    Core::LinAlg::Matrix<3, 3> s;
    Mat::add_contraction_matrix_four_tensor(s, 1.0, Id, stress);

    FOUR_C_EXPECT_NEAR(s, s_ref, 1.0e-10);
  }

  TEST(MaterialServiceTest, TestKroneckerProduct)
  {
    // test the calculation of the fourth order kronecker product of two second-order tensors
    Core::LinAlg::Matrix<3, 3> a(false);
    a(0, 0) = 1.0000000000;
    a(0, 1) = 2.0000000000;
    a(0, 2) = 3.0000000000;
    a(1, 0) = 0.0000000000;
    a(1, 1) = 5.0000000000;
    a(1, 2) = 1.2000000000;
    a(2, 0) = 1.0000000000;
    a(2, 1) = 1.0000000000;
    a(2, 2) = 1.0000000000;

    Core::LinAlg::Matrix<3, 3> b(false);
    b(0, 0) = 1.0000000000;
    b(0, 1) = 0.0000000000;
    b(0, 2) = 5.0000000000;
    b(1, 0) = 0.0000000000;
    b(1, 1) = 0.0000000000;
    b(1, 2) = 7.0000000000;
    b(2, 0) = 1.3000000000;
    b(2, 1) = 1.2000000000;
    b(2, 2) = 1.0000000000;

    Core::LinAlg::Matrix<6, 6> a_kron_b_ref(true);
    a_kron_b_ref(0, 0) = 1.0000000000;
    a_kron_b_ref(0, 1) = 0.0000000000;
    a_kron_b_ref(0, 2) = 15.0000000000;
    a_kron_b_ref(0, 3) = 1.0000000000;
    a_kron_b_ref(0, 4) = 5.0000000000;
    a_kron_b_ref(0, 5) = 4.0000000000;
    a_kron_b_ref(1, 0) = 0.0000000000;
    a_kron_b_ref(1, 1) = 0.0000000000;
    a_kron_b_ref(1, 2) = 8.4000000000;
    a_kron_b_ref(1, 3) = 0.0000000000;
    a_kron_b_ref(1, 4) = 17.5000000000;
    a_kron_b_ref(1, 5) = 0.0000000000;
    a_kron_b_ref(2, 0) = 1.3000000000;
    a_kron_b_ref(2, 1) = 1.2000000000;
    a_kron_b_ref(2, 2) = 1.0000000000;
    a_kron_b_ref(2, 3) = 1.2500000000;
    a_kron_b_ref(2, 4) = 1.1000000000;
    a_kron_b_ref(2, 5) = 1.1500000000;
    a_kron_b_ref(3, 0) = 0.0000000000;
    a_kron_b_ref(3, 1) = 0.0000000000;
    a_kron_b_ref(3, 2) = 21.0000000000;
    a_kron_b_ref(3, 3) = 0.0000000000;
    a_kron_b_ref(3, 4) = 7.0000000000;
    a_kron_b_ref(3, 5) = 3.5000000000;
    a_kron_b_ref(4, 0) = 0.0000000000;
    a_kron_b_ref(4, 1) = 6.0000000000;
    a_kron_b_ref(4, 2) = 1.2000000000;
    a_kron_b_ref(4, 3) = 3.2500000000;
    a_kron_b_ref(4, 4) = 3.2200000000;
    a_kron_b_ref(4, 5) = 0.7800000000;
    a_kron_b_ref(5, 0) = 1.3000000000;
    a_kron_b_ref(5, 1) = 2.4000000000;
    a_kron_b_ref(5, 2) = 3.0000000000;
    a_kron_b_ref(5, 3) = 1.9000000000;
    a_kron_b_ref(5, 4) = 2.8000000000;
    a_kron_b_ref(5, 5) = 2.4500000000;

    Core::LinAlg::Matrix<6, 6> a_kron_b(true);
    Mat::add_kronecker_tensor_product(a_kron_b, 1.0, a, b, 0.0);

    FOUR_C_EXPECT_NEAR(a_kron_b, a_kron_b_ref, 1.0e-10);
  }


  TEST(MaterialServiceTest, TestDyadicProductMatrixMatrix)
  {
    // test the dyadic product of two Matrices
    Core::LinAlg::Matrix<3, 3> stress(false);
    stress(0, 0) = 1.0;
    stress(1, 1) = 2.0;
    stress(2, 2) = 3.0;
    stress(0, 1) = stress(1, 0) = 0.4;
    stress(1, 2) = stress(2, 1) = 0.5;
    stress(0, 2) = stress(2, 0) = 0.6;

    Core::LinAlg::FourTensor<3> T;
    Mat::add_dyadic_product_matrix_matrix(T, 1.0, stress, stress);

    Core::LinAlg::Matrix<3, 3> result;
    Mat::add_contraction_matrix_four_tensor(result, 1.0, T, stress);

    Core::LinAlg::Matrix<3, 3> result_ref = stress;
    result_ref.scale(std::pow(stress.norm2(), 2));

    FOUR_C_EXPECT_NEAR(result, result_ref, 1.0e-10);
  }

}  // namespace
