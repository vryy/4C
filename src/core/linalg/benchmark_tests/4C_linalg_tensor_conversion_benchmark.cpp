// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor.hpp"
#include "4C_linalg_tensor_conversion.hpp"

#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstdlib>

FOUR_C_NAMESPACE_OPEN
template <std::size_t... size>
Core::LinAlg::Tensor<double, size...> get_random_tensor()
{
  Core::LinAlg::Tensor<double, size...> tensor;
  std::generate(tensor.container().begin(), tensor.container().end(),
      []() { return static_cast<double>(std::rand()) / RAND_MAX; });
  return tensor;
}

static void symmetric_four_tensor_to_6x6_voigt_matrix(benchmark::State& state)
{
  Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> A =
      Core::LinAlg::assume_symmetry(get_random_tensor<3, 3, 3, 3>());
  for (auto _ : state)
  {
    Core::LinAlg::Matrix<6, 6> result = Core::LinAlg::make_6x6_voigt_matrix_from_tensor(A);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(symmetric_four_tensor_to_6x6_voigt_matrix);

static void four_tensor_to_6x9_voigt_matrix(benchmark::State& state)
{
  Core::LinAlg::Tensor<double, 3, 3, 3, 3> A = get_random_tensor<3, 3, 3, 3>();
  for (auto _ : state)
  {
    Core::LinAlg::Matrix<6, 9> result = Core::LinAlg::make_6x9_voigt_matrix_from_tensor(A);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(four_tensor_to_6x9_voigt_matrix);

FOUR_C_NAMESPACE_CLOSE
