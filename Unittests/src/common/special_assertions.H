/*----------------------------------------------------------------------*/
/*! \file

\brief Special assertion functions that are often needed

\level 0

*/

#ifndef UNIT_SPECIAL_ASSERTIONS_H_
#define UNIT_SPECIAL_ASSERTIONS_H_

#include <src/linalg/linalg_fixedsizematrix.H>

namespace TESTING
{
  /*!
   * \brief Compares each entry of two matrices with a specified tolerance
   *
   * \tparam M rows of the matrices
   * \tparam N columns of the matrices
   * \tparam T type of the matrix entries
   * \param mat1 Matrix 1
   * \param mat2 Matrix 2
   * \param delta Maximum allowed tolerance
   */
  template <unsigned int M, unsigned int N, typename T>
  inline void AssertDelta(
      const LINALG::Matrix<M, N, T>& mat1, const LINALG::Matrix<M, N, T>& mat2, const T delta)
  {
    for (unsigned row = 0; row < M; ++row)
    {
      for (unsigned col = 0; col < N; ++col)
      {
        TS_ASSERT_DELTA(mat1(row, col), mat2(row, col), delta);
      }
    }
  }

  /*!
   * \brief Compares each entry of a std::vector with a specified tolerance
   *
   * \tparam T type of each entry
   * \param vector1 Vector 1
   * \param vector2 Vector 2
   * \param delta Maximum allowed tolerance
   */
  template <typename T>
  void AssertDelta(const std::vector<T>& vector1, const std::vector<T>& vector2, const T delta)
  {
    TS_ASSERT_EQUALS(vector1.size(), vector2.size());

    auto min_size = std::min(vector1.size(), vector2.size());

    for (auto i = 0U; i < min_size; ++i)
    {
      TS_ASSERT_DELTA(vector1[i], vector2[i], delta);
    }
  }
}  // namespace TESTING

#endif  // UNIT_SPECIAL_ASSERTIONS_H_