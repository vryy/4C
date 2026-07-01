// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_config.hpp"

#include "4C_linalg_tensor_conversion.hpp"

#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_fixedsizematrix_tensor_products.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_symmetric_tensor_eigen.hpp"
#include "4C_linalg_tensor_generators.hpp"
#include "4C_unittest_utils_assertions_test.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  TEST(TensorConversionTest, makeTensor)
  {
    Core::LinAlg::Matrix<3, 2> matrix;
    matrix(0, 0) = 1.0;
    matrix(0, 1) = 2.0;
    matrix(1, 0) = 3.0;
    matrix(1, 1) = 4.0;
    matrix(2, 0) = 5.0;
    matrix(2, 1) = 6.0;

    Core::LinAlg::Tensor<double, 3, 2> tensor = Core::LinAlg::make_tensor(matrix);
    EXPECT_DOUBLE_EQ(tensor(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(tensor(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(tensor(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(tensor(1, 1), 4.0);
    EXPECT_DOUBLE_EQ(tensor(2, 0), 5.0);
    EXPECT_DOUBLE_EQ(tensor(2, 1), 6.0);
  }

  TEST(TensorConversionTest, makeTensorView)
  {
    Core::LinAlg::Matrix<3, 2> matrix;
    matrix(0, 0) = 1.0;
    matrix(0, 1) = 2.0;
    matrix(1, 0) = 3.0;
    matrix(1, 1) = 4.0;
    matrix(2, 0) = 5.0;
    matrix(2, 1) = 6.0;

    Core::LinAlg::TensorView<double, 3, 2> tensor = Core::LinAlg::make_tensor_view(matrix);
    EXPECT_DOUBLE_EQ(tensor(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(tensor(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(tensor(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(tensor(1, 1), 4.0);
    EXPECT_DOUBLE_EQ(tensor(2, 0), 5.0);
    EXPECT_DOUBLE_EQ(tensor(2, 1), 6.0);

    tensor(2, 0) = 7.0;
    EXPECT_DOUBLE_EQ(matrix(2, 0), 7.0);
  }

  TEST(TensorConversionTest, makeConstTensorView)
  {
    const Core::LinAlg::Matrix<3, 2> matrix = []()
    {
      Core::LinAlg::Matrix<3, 2> matrix;
      matrix(0, 0) = 1.0;
      matrix(0, 1) = 2.0;
      matrix(1, 0) = 3.0;
      matrix(1, 1) = 4.0;
      matrix(2, 0) = 5.0;
      matrix(2, 1) = 6.0;
      return matrix;
    }();

    Core::LinAlg::TensorView<const double, 3, 2> tensor = Core::LinAlg::make_tensor_view(matrix);
    EXPECT_DOUBLE_EQ(tensor(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(tensor(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(tensor(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(tensor(1, 1), 4.0);
    EXPECT_DOUBLE_EQ(tensor(2, 0), 5.0);
    EXPECT_DOUBLE_EQ(tensor(2, 1), 6.0);

    EXPECT_EQ(tensor.data(), matrix.data());
  }


  TEST(TensorConversionTest, MakeSymmetricTensorFromStressLikeVoigtMatrix)
  {
    Core::LinAlg::Matrix<6, 1> stress_like_voigt;
    stress_like_voigt(0, 0) = 1.0;
    stress_like_voigt(1, 0) = 2.0;
    stress_like_voigt(2, 0) = 3.0;
    stress_like_voigt(3, 0) = 4.0;
    stress_like_voigt(4, 0) = 5.0;
    stress_like_voigt(5, 0) = 6.0;

    Core::LinAlg::SymmetricTensor<double, 3, 3> symmetric_tensor =
        Core::LinAlg::make_symmetric_tensor_from_stress_like_voigt_matrix(stress_like_voigt);

    EXPECT_DOUBLE_EQ(symmetric_tensor(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(symmetric_tensor(1, 1), 2.0);
    EXPECT_DOUBLE_EQ(symmetric_tensor(2, 2), 3.0);
    EXPECT_DOUBLE_EQ(symmetric_tensor(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(symmetric_tensor(1, 2), 5.0);
    EXPECT_DOUBLE_EQ(symmetric_tensor(0, 2), 6.0);
  }

  TEST(TensorConversionTest, reinterpretAsTensor)
  {
    Core::LinAlg::Matrix<3, 1> matrix;
    matrix(0, 0) = 1.0;
    matrix(1, 0) = 2.0;
    matrix(2, 0) = 3.0;

    Core::LinAlg::Tensor<double, 3> tensor = Core::LinAlg::reinterpret_as_tensor<3>(matrix);
    EXPECT_DOUBLE_EQ(tensor(0), 1.0);
    EXPECT_DOUBLE_EQ(tensor(1), 2.0);
    EXPECT_DOUBLE_EQ(tensor(2), 3.0);
  }

  TEST(TensorConversionTest, reinterpretAsTensorView)
  {
    const Core::LinAlg::Matrix<3, 1> matrix = []()
    {
      Core::LinAlg::Matrix<3, 1> matrix;
      matrix(0, 0) = 1.0;
      matrix(1, 0) = 2.0;
      matrix(2, 0) = 3.0;
      return matrix;
    }();

    Core::LinAlg::TensorView<const double, 3> tensor =
        Core::LinAlg::reinterpret_as_tensor_view<3>(matrix);
    EXPECT_DOUBLE_EQ(tensor(0), 1.0);
    EXPECT_DOUBLE_EQ(tensor(1), 2.0);
    EXPECT_DOUBLE_EQ(tensor(2), 3.0);

    EXPECT_EQ(tensor.data(), matrix.data());
  }

  TEST(TensorConversionTest, reinterpretAsConstTensorView)
  {
    Core::LinAlg::Matrix<3, 1> matrix;
    matrix(0, 0) = 1.0;
    matrix(1, 0) = 2.0;
    matrix(2, 0) = 3.0;

    Core::LinAlg::TensorView<double, 3> tensor =
        Core::LinAlg::reinterpret_as_tensor_view<3>(matrix);
    EXPECT_DOUBLE_EQ(tensor(0), 1.0);
    EXPECT_DOUBLE_EQ(tensor(1), 2.0);
    EXPECT_DOUBLE_EQ(tensor(2), 3.0);

    tensor(2) = 7.0;
    EXPECT_DOUBLE_EQ(matrix(2, 0), 7.0);
  }

  TEST(TensorConversionTest, makeMatrixView)
  {
    Core::LinAlg::Tensor<double, 3, 2> tensor = {{{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}}};

    Core::LinAlg::Matrix<3, 2> matrix_view = Core::LinAlg::make_matrix_view(tensor);
    EXPECT_DOUBLE_EQ(matrix_view(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(matrix_view(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(matrix_view(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(matrix_view(1, 1), 4.0);
    EXPECT_DOUBLE_EQ(matrix_view(2, 0), 5.0);
    EXPECT_DOUBLE_EQ(matrix_view(2, 1), 6.0);

    matrix_view(2, 0) = 7.0;
    EXPECT_DOUBLE_EQ(tensor(2, 0), 7.0);
  }

  TEST(TensorConversionTest, makeMatrixViewFromArray)
  {
    std::array<Core::LinAlg::Tensor<double, 2>, 3> array_of_tensors;
    array_of_tensors[0] = {{1.0, 2.0}};
    array_of_tensors[1] = {{3.0, 4.0}};
    array_of_tensors[2] = {{5.0, 6.0}};

    Core::LinAlg::Matrix<2, 3> matrix_view = Core::LinAlg::make_matrix_view(array_of_tensors);
    EXPECT_DOUBLE_EQ(matrix_view(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(matrix_view(0, 1), 3.0);
    EXPECT_DOUBLE_EQ(matrix_view(0, 2), 5.0);
    EXPECT_DOUBLE_EQ(matrix_view(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(matrix_view(1, 1), 4.0);
    EXPECT_DOUBLE_EQ(matrix_view(1, 2), 6.0);

    matrix_view(1, 2) = 7.0;
    EXPECT_DOUBLE_EQ(array_of_tensors[2](1), 7.0);
  }

  TEST(TensorConversionTest, makeMatrixViewReinterpretation)
  {
    Core::LinAlg::Tensor<double, 3> tensor = {{1.0, 2.0, 3.0}};

    Core::LinAlg::Matrix<3, 1> matrix_view = Core::LinAlg::make_matrix_view<3, 1>(tensor);
    EXPECT_DOUBLE_EQ(matrix_view(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(matrix_view(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(matrix_view(2, 0), 3.0);

    matrix_view(2, 0) = 7.0;
    EXPECT_DOUBLE_EQ(tensor(2), 7.0);
  }

  TEST(TensorConversionTest, makeMatrix)
  {
    Core::LinAlg::Tensor<double, 3, 2> tensor = {{{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}}};

    Core::LinAlg::Matrix<3, 2> matrix = Core::LinAlg::make_matrix(tensor);
    EXPECT_DOUBLE_EQ(matrix(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(matrix(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(matrix(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(matrix(1, 1), 4.0);
    EXPECT_DOUBLE_EQ(matrix(2, 0), 5.0);
    EXPECT_DOUBLE_EQ(matrix(2, 1), 6.0);

    matrix(2, 0) = 7.0;
    EXPECT_DOUBLE_EQ(tensor(2, 0), 5.0);
  }

  TEST(TensorConversionTest, makeMatrixReinterpretation)
  {
    Core::LinAlg::Tensor<double, 3> tensor = {{1.0, 2.0, 3.0}};

    Core::LinAlg::Matrix<3, 1> matrix = Core::LinAlg::make_matrix<3, 1>(tensor);
    EXPECT_DOUBLE_EQ(matrix(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(matrix(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(matrix(2, 0), 3.0);

    matrix(2, 0) = 7.0;
    EXPECT_DOUBLE_EQ(tensor(2), 3.0);
  }

  TEST(TensorConversionTest, MakeStressLikeVoigtView)
  {
    Core::LinAlg::SymmetricTensor<double, 3, 3> tensor =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 3, 3>{{
            {1.0, 4.0, 6.0},
            {4.0, 2.0, 5.0},
            {6.0, 5.0, 3.0},
        }});

    Core::LinAlg::Matrix<6, 1> voigt_view = Core::LinAlg::make_stress_like_voigt_view(tensor);
    EXPECT_DOUBLE_EQ(voigt_view(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(voigt_view(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(voigt_view(2, 0), 3.0);
    EXPECT_DOUBLE_EQ(voigt_view(3, 0), 4.0);
    EXPECT_DOUBLE_EQ(voigt_view(4, 0), 5.0);
    EXPECT_DOUBLE_EQ(voigt_view(5, 0), 6.0);

    voigt_view(3, 0) = 7.0;
    EXPECT_DOUBLE_EQ(tensor(0, 1), 7.0);
  }

  TEST(TensorConversionTest, MakeStrainLikeVoigtMatrix)
  {
    Core::LinAlg::SymmetricTensor<double, 3, 3> tensor =
        Core::LinAlg::assume_symmetry(Core::LinAlg::Tensor<double, 3, 3>{{
            {1.0, 4.0, 6.0},
            {4.0, 2.0, 5.0},
            {6.0, 5.0, 3.0},
        }});

    Core::LinAlg::Matrix<6, 1> voigt_strain = Core::LinAlg::make_strain_like_voigt_matrix(tensor);
    EXPECT_DOUBLE_EQ(voigt_strain(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(voigt_strain(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(voigt_strain(2, 0), 3.0);
    EXPECT_DOUBLE_EQ(voigt_strain(3, 0), 8.0);
    EXPECT_DOUBLE_EQ(voigt_strain(4, 0), 10.0);
    EXPECT_DOUBLE_EQ(voigt_strain(5, 0), 12.0);

    EXPECT_DOUBLE_EQ(tensor(0, 1), 4.0);
  }

  TEST(TensorConversionTest, MakeTensorFromNestedArray)
  {
    std::array<std::array<double, 2>, 3> nested_array = {{
        {1.0, 2.0},
        {3.0, 4.0},
        {5.0, 6.0},
    }};
    Core::LinAlg::Tensor<double, 3, 2> tensor =
        Core::LinAlg::make_tensor_from_nested_array<double, 3, 2>(nested_array);


    EXPECT_DOUBLE_EQ(tensor(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(tensor(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(tensor(2, 0), 5.0);
    EXPECT_DOUBLE_EQ(tensor(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(tensor(1, 1), 4.0);
    EXPECT_DOUBLE_EQ(tensor(2, 1), 6.0);
  }

  TEST(TensorConversionTest, MakeNestedArrayFromTensor)
  {
    Core::LinAlg::Tensor<double, 3, 2> tensor = {{
        {1.0, 2.0},
        {3.0, 4.0},
        {5.0, 6.0},
    }};
    std::array<std::array<double, 2>, 3> nested_array =
        Core::LinAlg::make_nested_array_from_tensor(tensor);


    EXPECT_DOUBLE_EQ(nested_array[0][0], 1.0);
    EXPECT_DOUBLE_EQ(nested_array[1][0], 3.0);
    EXPECT_DOUBLE_EQ(nested_array[2][0], 5.0);
    EXPECT_DOUBLE_EQ(nested_array[0][1], 2.0);
    EXPECT_DOUBLE_EQ(nested_array[1][1], 4.0);
    EXPECT_DOUBLE_EQ(nested_array[2][1], 6.0);
  }
  TEST(TensorConversionTest, Make6x9VoigtMatrixFromTensorTest)
  {
    // arbitrary (non-symmetric) 4th order tensor; both paths apply the same 1st-pair symmetrization
    Core::LinAlg::Tensor<double, 3, 3, 3, 3> T{};
    Core::LinAlg::FourTensor<3> ft(true);
    int v = 1;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
          for (int l = 0; l < 3; ++l)
          {
            const double val = 0.1 * (v++);
            T(i, j, k, l) = val;
            ft(i, j, k, l) = val;
          }

    const Core::LinAlg::Matrix<6, 9> new_voigt = Core::LinAlg::make_6x9_voigt_matrix_from_tensor(T);
    Core::LinAlg::Matrix<6, 9> old_voigt(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Voigt::setup_6x9_voigt_matrix_from_four_tensor(old_voigt, ft);

    FOUR_C_EXPECT_NEAR(new_voigt, old_voigt, 1e-12);
  }

  TEST(TensorConversionTest, Make9x6VoigtMatrixFromTensorTest)
  {
    // arbitrary (non-symmetric) 4th order tensor; both paths apply the same 2nd-pair symmetrization
    Core::LinAlg::Tensor<double, 3, 3, 3, 3> T{};
    Core::LinAlg::FourTensor<3> ft(true);
    int v = 1;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
          for (int l = 0; l < 3; ++l)
          {
            const double val = 0.1 * (v++);
            T(i, j, k, l) = val;
            ft(i, j, k, l) = val;
          }

    const Core::LinAlg::Matrix<9, 6> new_voigt = Core::LinAlg::make_9x6_voigt_matrix_from_tensor(T);
    Core::LinAlg::Matrix<9, 6> old_voigt(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Voigt::setup_9x6_voigt_matrix_from_four_tensor(old_voigt, ft);

    FOUR_C_EXPECT_NEAR(new_voigt, old_voigt, 1e-12);
  }

  TEST(TensorConversionTest, Make6x6StressLikeVoigtMatrixFromTensorTest)
  {
    Core::LinAlg::Matrix<3, 3> A;
    A(0, 0) = 1.0;
    A(0, 1) = 2.0;
    A(0, 2) = 3.0;
    A(1, 0) = 0.5;
    A(1, 1) = 4.0;
    A(1, 2) = 1.2;
    A(2, 0) = 7.0;
    A(2, 1) = 0.3;
    A(2, 2) = 6.0;

    Core::LinAlg::Tensor<double, 3, 3, 3, 3> C_full{};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
          for (int l = 0; l < 3; ++l)
            C_full(i, j, k, l) = 0.5 * (A(i, k) * A(j, l) + A(i, l) * A(j, k));
    const Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> C =
        Core::LinAlg::assume_symmetry(C_full);

    const Core::LinAlg::Matrix<6, 6> result = Core::LinAlg::make_6x6_voigt_matrix_from_tensor(C);

    // reference: add_kronecker_tensor_product called on non-symmetric matrix
    Core::LinAlg::Matrix<6, 6> ref(Core::LinAlg::Initialization::zero);
    Core::LinAlg::FourTensorOperations::add_kronecker_tensor_product(ref, 1.0, A, A, 0.0);

    FOUR_C_EXPECT_NEAR(result, ref, 1e-12);
  }
}  // namespace

FOUR_C_NAMESPACE_CLOSE