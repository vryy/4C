/*----------------------------------------------------------------------*/
/*! \file
\brief special assertions for BACI code
\level 1
*----------------------------------------------------------------------*/
#ifndef FOUR_C_UNITTEST_UTILS_ASSERTIONS_TEST_HPP
#define FOUR_C_UNITTEST_UTILS_ASSERTIONS_TEST_HPP

#include <gtest/gtest.h>

#include <gmock/gmock.h>

#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_linalg_serialdensematrix.hpp"

#include <type_traits>

namespace TESTING::INTERNAL
{
  namespace
  {
    using namespace FourC;

    //! Determine a number of decimal digits for printing to a given tolerance.
    template <typename T>
    inline int PrecisionForPrinting(T tolerance)
    {
      FOUR_C_ASSERT(tolerance > 0, "Tolerance must be positive.");
      return tolerance > 1 ? 0 : static_cast<int>(std::ceil(-1.0 * std::log10(tolerance) + 1));
    }

    //! If the @p nonMatchingEntries string is empty return success, otherwise return failure with a
    //! descriptive message.
    template <typename T>
    inline ::testing::AssertionResult ResultBasedOnNonMatchingEntries(
        const std::string& nonMatchingEntries, T tolerance, const char* expr1, const char* expr2)
    {
      if (nonMatchingEntries.empty())
        return ::testing::AssertionSuccess();
      else
      {
        return ::testing::AssertionFailure()
               << "The following entries differ: absolute difference is not within tolerance "
               << tolerance << "\n"
               << "index: " << expr1 << " vs. " << expr2 << ":" << std::endl
               << nonMatchingEntries;
      }
    }
  }  // namespace

  /**
   * Compare two iterable objects for double equality up to a tolerance. The iterables must have the
   * same underlying data types. The signature is mandated by GoogleTest's EXPECT_PRED_FORMAT4
   * macro.
   *
   * @note This function is not intended to be used directly. Use FOUR_C_EXPECT_ITERABLE_NEAR.
   */
  template <typename Iterator1, typename Iterator2>
  inline ::testing::AssertionResult AssertNear(const char* vec1Expr, const char* vec2Expr,
      const char* /*lengthExpr*/, const char* /*toleranceExpr*/, const Iterator1 iter1,
      const Iterator2 iter2, std::size_t length,
      std::decay_t<decltype(*std::declval<Iterator1>())> tolerance)
  {
    static_assert(std::is_same_v<std::decay_t<decltype(*iter1)>, std::decay_t<decltype(*iter2)>>,
        "Underlying data types are not the same.");

    const std::string nonMatchingEntries = std::invoke(
        [&]()
        {
          std::stringstream ss;
          ss << std::fixed << std::setprecision(PrecisionForPrinting(tolerance));
          for (std::size_t i = 0; i < length; ++i)
          {
            const auto value1 = *(iter1 + i);
            const auto value2 = *(iter2 + i);
            if (std::fabs(value1 - value2) > tolerance)
            {
              ss << "[" << i << "]: " << value1 << " vs. " << value2 << std::endl;
            }
          }
          return ss.str();
        });

    return ResultBasedOnNonMatchingEntries(nonMatchingEntries, tolerance, vec1Expr, vec2Expr);
  }

  /**
   * Compare two std::vector<T> objects for double equality up to a tolerance. The signature is
   * mandated by GoogleTest's EXPECT_PRED_FORMAT3 macro.
   *
   * @note This function is not intended to be used directly. Use FOUR_C_EXPECT_NEAR.
   */
  template <typename T>
  inline ::testing::AssertionResult AssertNear(const char* vec1Expr,
      const char* vec2Expr,  // NOLINT
      const char* toleranceExpr, const std::vector<T>& vec1, const std::vector<T>& vec2,
      T tolerance)
  {
    // argument is required for the EXPECT_PRED_FORMAT3 macro of GoogleTest for pretty printing
    (void)toleranceExpr;

    if (vec1.size() != vec2.size())
    {
      return ::testing::AssertionFailure()
             << "size mismatch: " << vec1Expr << " has size " << vec1.size() << " but " << vec2Expr
             << " has dimension " << vec2.size() << std::endl;
    }

    return AssertNear(vec1Expr, vec2Expr, "" /*lengthExpr*/, toleranceExpr, vec1.begin(),
        vec2.begin(), vec1.size(), tolerance);
  }

  /**
   * Compare two CORE::LINALG::Matrix objects for double equality up to a tolerance. The signature
   * is mandated by GoogleTest's EXPECT_PRED_FORMAT3 macro.
   *
   * @note This function is not intended to be used directly. Use FOUR_C_EXPECT_NEAR.
   */
  template <unsigned int M, unsigned int N, typename T>
  inline ::testing::AssertionResult AssertNear(const char* mat1Expr,
      const char* mat2Expr,  // NOLINT
      const char* toleranceExpr, const CORE::LINALG::Matrix<M, N, T>& mat1,
      const CORE::LINALG::Matrix<M, N, T>& mat2, T tolerance)
  {
    // argument is required for the EXPECT_PRED_FORMAT3 macro of GoogleTest for pretty printing
    (void)toleranceExpr;

    const std::string nonMatchingEntries = std::invoke(
        [&]()
        {
          std::stringstream ss;
          ss << std::fixed << std::setprecision(PrecisionForPrinting(tolerance));
          for (unsigned i = 0; i < M; ++i)
          {
            for (unsigned j = 0; j < N; ++j)
            {
              if (std::fabs(mat1(i, j) - mat2(i, j)) > tolerance)
              {
                ss << "(" << i << "," << j << "): " << mat1(i, j) << " vs. " << mat2(i, j)
                   << std::endl;
              }
            }
          }
          return ss.str();
        });

    return ResultBasedOnNonMatchingEntries(nonMatchingEntries, tolerance, mat1Expr, mat2Expr);
  }

  /**
   * Compare a CORE::LINALG::Matrix with a std::array for double equality up to a tolerance. The
   * entries in std::array row-major. The signature is mandated by GoogleTest's EXPECT_PRED_FORMAT3
   * macro.
   *
   * @note This function is not intended to be used directly. Use FOUR_C_EXPECT_NEAR.
   */
  template <unsigned int M, unsigned int N, typename T>
  inline ::testing::AssertionResult AssertNear(const char* mat1Expr,
      const char* mat2Expr,  // NOLINT
      const char* toleranceExpr, const CORE::LINALG::Matrix<M, N, T>& mat,
      const std::array<T, static_cast<std::size_t>(M) * N>& array, T tolerance)
  {
    // argument is required for the EXPECT_PRED_FORMAT3 macro of GoogleTest for pretty printing
    (void)toleranceExpr;

    const std::string nonMatchingEntries = std::invoke(
        [&]()
        {
          std::stringstream ss;
          ss << std::fixed << std::setprecision(PrecisionForPrinting(tolerance));
          for (unsigned i = 0; i < M; ++i)
          {
            for (unsigned j = 0; j < N; ++j)
            {
              const std::size_t arr_index = i * N + j;
              if (std::fabs(mat(i, j) - array[arr_index]) > tolerance)
              {
                ss << "(" << i << "," << j << ") vs. [" << arr_index << "]: " << mat(i, j)
                   << " vs. " << array[arr_index] << std::endl;
              }
            }
          }
          return ss.str();
        });

    return ResultBasedOnNonMatchingEntries(nonMatchingEntries, tolerance, mat1Expr, mat2Expr);
  }

  /**
   * Compare two CORE::LINALG::SerialDenseMatrix objects for double equality up to a tolerance. The
   * signature is mandated by GoogleTest's EXPECT_PRED_FORMAT3 macro.
   *
   * @note This function is not intended to be used directly. Use FOUR_C_EXPECT_NEAR.
   */
  inline ::testing::AssertionResult AssertNear(const char* mat1Expr,
      const char* mat2Expr,  // NOLINT
      const char* toleranceExpr, const CORE::LINALG::SerialDenseMatrix& mat1,
      const CORE::LINALG::SerialDenseMatrix& mat2, double tolerance)
  {
    // argument is required for the EXPECT_PRED_FORMAT3 macro of GoogleTest for pretty printing
    (void)toleranceExpr;

    const bool dimensionsMatch =
        mat1.numRows() == mat2.numRows() and mat1.numCols() == mat2.numCols();
    if (!dimensionsMatch)
    {
      return ::testing::AssertionFailure()
             << "dimension mismatch: " << mat1Expr << " has dimension " << mat1.numRows() << "x"
             << mat1.numCols() << " but " << mat2Expr << " has dimension " << mat2.numRows() << "x"
             << mat2.numCols() << std::endl;
    }

    const std::string nonMatchingEntries = std::invoke(
        [&]()
        {
          std::stringstream ss;
          ss << std::fixed << std::setprecision(PrecisionForPrinting(tolerance));
          for (int i = 0; i < mat1.numRows(); ++i)
          {
            for (int j = 0; j < mat1.numCols(); ++j)
            {
              if (std::fabs(mat1(i, j) - mat2(i, j)) > tolerance)
              {
                ss << "(" << i << "," << j << "): " << mat1(i, j) << " vs. " << mat2(i, j)
                   << std::endl;
              }
            }
          }
          return ss.str();
        });

    return ResultBasedOnNonMatchingEntries(nonMatchingEntries, tolerance, mat1Expr, mat2Expr);
  }

  /**
   * Compare two Teuchos::SerialDenseMatrix objects for double equality up to a tolerance. The
   * signature is mandated by GoogleTest's EXPECT_PRED_FORMAT3 macro.
   *
   * @note This function is not intended to be used directly. Use FOUR_C_EXPECT_NEAR.
   */
  template <typename T>
  inline ::testing::AssertionResult AssertNear(const char* mat1Expr,
      const char* mat2Expr,  // NOLINT
      const char* toleranceExpr, const Teuchos::SerialDenseMatrix<int, T>& mat1,
      const Teuchos::SerialDenseMatrix<int, T>& mat2, double tolerance)
  {
    // argument is required for the EXPECT_PRED_FORMAT3 macro of GoogleTest for pretty printing
    (void)toleranceExpr;

    const bool dimensionsMatch =
        mat1.numRows() == mat2.numRows() and mat1.numCols() == mat2.numCols();
    if (!dimensionsMatch)
    {
      return ::testing::AssertionFailure()
             << "dimension mismatch: " << mat1Expr << " has dimension " << mat1.numRows() << "x"
             << mat1.numCols() << " but " << mat2Expr << " has dimension " << mat2.numRows() << "x"
             << mat2.numCols() << std::endl;
    }

    const std::string nonMatchingEntries = std::invoke(
        [&]()
        {
          std::stringstream ss;
          ss << std::fixed << std::setprecision(PrecisionForPrinting(tolerance));
          for (int i = 0; i < mat1.numRows(); ++i)
          {
            for (int j = 0; j < mat1.numCols(); ++j)
            {
              if (std::abs(mat1(i, j) - mat2(i, j)) > tolerance)
              {
                ss << "(" << i << "," << j << "): " << mat1(i, j) << " vs. " << mat2(i, j)
                   << std::endl;
              }
            }
          }
          return ss.str();
        });

    return ResultBasedOnNonMatchingEntries(nonMatchingEntries, tolerance, mat1Expr, mat2Expr);
  }

}  // namespace TESTING::INTERNAL

/**
 * @brief Custom assertion to test for equality up to a tolerance.
 *
 * This macro tests two containers @p actual and @p expected for entry-wise equality up to an
 * absolute @p tolerance. Overloads are provided for
 * - std::vector<T>
 * - CORE::LINALG::Matrix
 * - CORE::LINALG::SerialDenseMatrix
 *
 * @note Implementation details: this and similar macros are defined to avoid writing asserts in the
 * unexpressive EXPECT_PRED_FORMATn syntax by gtest. They are all prefixed with `FOUR_C_` to easily
 * distinguish them from gtest asserts.
 */
#define FOUR_C_EXPECT_NEAR(actual, expected, tolerance) \
  EXPECT_PRED_FORMAT3(TESTING::INTERNAL::AssertNear, actual, expected, tolerance)

/**
 * @brief Custom assertion to test for equality up to a tolerance.
 *
 * This macro tests two iterables @p actual and @p expected with length @p length for entry-wise
 * equality up to an absolute @p tolerance. This works for any Iterator and also raw arrays or
 * pointers. Also, different iterables can be compared e.g., std::vector with raw array. The
 * underlying data type must be the same.
 * When testing containers you need to pass an iterator e.g. my_vector.begin()
 *
 * @note Implementation details: this and similar macros are defined to avoid writing asserts in the
 * unexpressive EXPECT_PRED_FORMATn syntax by gtest. They are all prefixed with `FOUR_C_` to easily
 * distinguish them from gtest asserts.
 */
#define FOUR_C_EXPECT_ITERABLE_NEAR(actual, expected, length, tolerance) \
  EXPECT_PRED_FORMAT4(TESTING::INTERNAL::AssertNear, actual, expected, length, tolerance)

/*!
 * Extension of EXPECT_THROW which also checks for a substring in the what() expression.
 */
#define FOUR_C_EXPECT_THROW_WITH_MESSAGE(statement, expectedException, messageSubString)        \
  {                                                                                             \
    std::function<void(void)> find_the_statement_below(                                         \
        [&]()                                                                                   \
        {                                                                                       \
          try                                                                                   \
          {                                                                                     \
            statement;                                                                          \
          }                                                                                     \
          catch (const expectedException& caughtException)                                      \
          {                                                                                     \
            using ::testing::HasSubstr;                                                         \
            EXPECT_THAT(caughtException.what(), HasSubstr(messageSubString))                    \
                << "Caught the expected exception type but message has wrong substring.";       \
            throw;                                                                              \
          }                                                                                     \
        });                                                                                     \
    EXPECT_THROW(find_the_statement_below(), expectedException) << "statement: " << #statement; \
  }

#endif
