/*---------------------------------------------------------------------------*/
/*! \file

\brief Unittests for utilities concerning voigt notation

\level 3
*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>
#include "baci_linalg_fixedsizematrix.H"
#include "baci_lib_voigt_notation.H"
#include "baci_unittest_utils_assertions.h"

namespace
{
  class VoigtNotationTest : public testing::Test
  {
   public:
    CORE::LINALG::Matrix<3, 3> tens;
    CORE::LINALG::Matrix<3, 3> itens;
    CORE::LINALG::Matrix<6, 1> tens_strain;
    CORE::LINALG::Matrix<6, 1> itens_strain;
    CORE::LINALG::Matrix<6, 1> tens_stress;
    CORE::LINALG::Matrix<6, 1> itens_stress;

    VoigtNotationTest()
    {
      tens(0, 0) = 1.1;
      tens(1, 1) = 1.2;
      tens(2, 2) = 1.3;
      tens(0, 1) = tens(1, 0) = 0.01;
      tens(1, 2) = tens(2, 1) = 0.02;
      tens(0, 2) = tens(2, 0) = 0.03;

      itens(0, 0) = 0.90972618;
      itens(1, 1) = 0.83360457;
      itens(2, 2) = 0.76990741;
      itens(0, 1) = itens(1, 0) = -0.00723301;
      itens(1, 2) = itens(2, 1) = -0.01265777;
      itens(0, 2) = itens(2, 0) = -0.0208824;

      tens_strain(0) = tens_stress(0) = 1.1;
      tens_strain(1) = tens_stress(1) = 1.2;
      tens_strain(2) = tens_stress(2) = 1.3;

      tens_stress(3) = 0.01;
      tens_stress(4) = 0.02;
      tens_stress(5) = 0.03;

      itens_stress(0) = 0.90972618;
      itens_stress(1) = 0.83360457;
      itens_stress(2) = 0.76990741;
      itens_stress(3) = -0.00723301;
      itens_stress(4) = -0.01265777;
      itens_stress(5) = -0.0208824;

      tens_strain(3) = 2 * 0.01;
      tens_strain(4) = 2 * 0.02;
      tens_strain(5) = 2 * 0.03;

      itens_strain(0) = 0.90972618;
      itens_strain(1) = 0.83360457;
      itens_strain(2) = 0.76990741;
      itens_strain(3) = 2 * -0.00723301;
      itens_strain(4) = 2 * -0.01265777;
      itens_strain(5) = 2 * -0.0208824;
    };
  };

  TEST_F(VoigtNotationTest, MatrixToVectorStressLike)
  {
    CORE::LINALG::Matrix<6, 1> cmp_stress(false);

    UTILS::VOIGT::Stresses::MatrixToVector(tens, cmp_stress);

    BACI_EXPECT_NEAR(cmp_stress, tens_stress, 1e-10);
  }

  TEST_F(VoigtNotationTest, MatrixToVectorStrainLike)
  {
    CORE::LINALG::Matrix<6, 1> cmp_strain(false);

    UTILS::VOIGT::Strains::MatrixToVector(tens, cmp_strain);

    BACI_EXPECT_NEAR(cmp_strain, tens_strain, 1e-10);
  }

  TEST_F(VoigtNotationTest, DeterminantStressLike)
  {
    EXPECT_NEAR(UTILS::VOIGT::Stresses::Determinant(tens_stress), 1.7143620000000002, 1e-10);
  }

  TEST_F(VoigtNotationTest, DeterminantStrainLike)
  {
    EXPECT_NEAR(UTILS::VOIGT::Strains::Determinant(tens_strain), 1.7143620000000002, 1e-10);
  }

  TEST_F(VoigtNotationTest, InvariantsPrincipalStressLike)
  {
    CORE::LINALG::Matrix<3, 1> prinv(false);
    UTILS::VOIGT::Stresses::InvariantsPrincipal(prinv, tens_stress);
    EXPECT_NEAR(prinv(0), 3.5999999999999996, 1e-10);
    EXPECT_NEAR(prinv(1), 4.3085999999999984, 1e-10);
    EXPECT_NEAR(prinv(2), 1.7143620000000002, 1e-10);
  }

  TEST_F(VoigtNotationTest, InvariantsPrincipalStrainLike)
  {
    CORE::LINALG::Matrix<3, 1> prinv(false);
    UTILS::VOIGT::Strains::InvariantsPrincipal(prinv, tens_strain);
    EXPECT_NEAR(prinv(0), 3.5999999999999996, 1e-10);
    EXPECT_NEAR(prinv(1), 4.3085999999999984, 1e-10);
    EXPECT_NEAR(prinv(2), 1.7143620000000002, 1e-10);
  }

  TEST_F(VoigtNotationTest, InverseStressLike)
  {
    CORE::LINALG::Matrix<6, 1> itens_stress_result(false);
    UTILS::VOIGT::Stresses::InverseTensor(tens_stress, itens_stress_result);

    BACI_EXPECT_NEAR(itens_stress_result, itens_stress, 1e-5);
  }

  TEST_F(VoigtNotationTest, InverseStrainLike)
  {
    CORE::LINALG::Matrix<6, 1> itens_strain_result(false);
    UTILS::VOIGT::Strains::InverseTensor(tens_strain, itens_strain_result);

    BACI_EXPECT_NEAR(itens_strain_result, itens_strain, 1e-5);
  }

  TEST_F(VoigtNotationTest, ToStressLike)
  {
    CORE::LINALG::Matrix<6, 1> strain_to_stress(false);
    CORE::LINALG::Matrix<6, 1> stress_to_stress(false);
    UTILS::VOIGT::Strains::ToStressLike(tens_strain, strain_to_stress);
    UTILS::VOIGT::Stresses::ToStressLike(tens_stress, stress_to_stress);

    BACI_EXPECT_NEAR(strain_to_stress, tens_stress, 1e-5);
    BACI_EXPECT_NEAR(stress_to_stress, stress_to_stress, 1e-5);
  }

  TEST_F(VoigtNotationTest, ToStrainLike)
  {
    CORE::LINALG::Matrix<6, 1> strain_to_strain(false);
    CORE::LINALG::Matrix<6, 1> stress_to_strain(false);
    UTILS::VOIGT::Strains::ToStrainLike(tens_strain, strain_to_strain);
    UTILS::VOIGT::Stresses::ToStrainLike(tens_stress, stress_to_strain);


    BACI_EXPECT_NEAR(strain_to_strain, tens_strain, 1e-5);
    BACI_EXPECT_NEAR(stress_to_strain, tens_strain, 1e-5);
  }

  TEST_F(VoigtNotationTest, IdentityMatrix)
  {
    CORE::LINALG::Matrix<6, 1> id(false);

    UTILS::VOIGT::IdentityMatrix(id);

    EXPECT_NEAR(id(0), 1.0, 1e-10);
    EXPECT_NEAR(id(1), 1.0, 1e-10);
    EXPECT_NEAR(id(2), 1.0, 1e-10);
    EXPECT_NEAR(id(3), 0.0, 1e-10);
    EXPECT_NEAR(id(4), 0.0, 1e-10);
    EXPECT_NEAR(id(5), 0.0, 1e-10);
  }

  TEST_F(VoigtNotationTest, StrainLikeVectorToMatrix)
  {
    CORE::LINALG::Matrix<3, 3> matrix(false);
    UTILS::VOIGT::Strains::VectorToMatrix(tens_strain, matrix);
    EXPECT_EQ(matrix, tens);
  }

  TEST_F(VoigtNotationTest, StressLikeVectorToMatrix)
  {
    CORE::LINALG::Matrix<3, 3> matrix(false);
    UTILS::VOIGT::Stresses::VectorToMatrix(tens_stress, matrix);
    EXPECT_EQ(matrix, tens);
  }
}  // namespace