/*----------------------------------------------------------------------*/
/*! \file
\brief special assertions for BACI code
\level 1
*----------------------------------------------------------------------*/
#ifndef UNITTESTS_COMMON_ASSERTIONS_H
#define UNITTESTS_COMMON_ASSERTIONS_H

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "src/linalg/linalg_fixedsizematrix.H"
#include "src/linalg/linalg_serialdensematrix.H"

namespace TESTING::INTERNAL
{

  //! Determine a number of decimal digits for printing to a given tolerance.
  template <typename T>
  inline int PrecisionForPrinting(T tolerance)
  {
    dsassert(tolerance > 0, "Tolerance must be positive.");
    return tolerance > 1 ? 0 : static_cast<int>(std::ceil(-1.0 * std::log10(tolerance) + 1));
  }

  /**
   * Compare two std::vector<T> objects for double equality up to a tolerance. The  signature is
   * mandated by GoogleTest's EXPECT_PRED_FORMAT3 macro.
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

    bool allEntriesMatch = true;
    std::stringstream nonMatchingEntries;
    nonMatchingEntries << std::fixed << std::setprecision(PrecisionForPrinting(tolerance));
    for (std::size_t i = 0; i < vec1.size(); ++i)
    {
      if (std::fabs(vec1[i] - vec2[i]) > tolerance)
      {
        allEntriesMatch = false;
        nonMatchingEntries << "[" << i << "]: " << vec1[i] << " vs. " << vec2[i] << std::endl;
      }
    }

    if (allEntriesMatch)
      return ::testing::AssertionSuccess();
    else
    {
      return ::testing::AssertionFailure()
             << "The following entries differ: absolute difference is not within tolerance "
             << tolerance << "\n"
             << "[i]: " << vec1Expr << " vs. " << vec2Expr << ":" << std::endl
             << nonMatchingEntries.str();
    }
  }

  /**
   * Compare two LINALG::Matrix objects for double equality up to a tolerance. The signature is
   * mandated by GoogleTest's EXPECT_PRED_FORMAT3 macro.
   */
  template <unsigned int M, unsigned int N, typename T>
  inline ::testing::AssertionResult AssertNear(const char* mat1Expr,
      const char* mat2Expr,  // NOLINT
      const char* toleranceExpr, const LINALG::Matrix<M, N, T>& mat1,
      const LINALG::Matrix<M, N, T>& mat2, T tolerance)
  {
    // argument is required for the EXPECT_PRED_FORMAT3 macro of GoogleTest for pretty printing
    (void)toleranceExpr;

    bool allEntriesMatch = true;
    std::stringstream nonMatchingEntries;
    nonMatchingEntries << std::fixed << std::setprecision(PrecisionForPrinting(tolerance));
    for (unsigned i = 0; i < M; ++i)
    {
      for (unsigned j = 0; j < N; ++j)
      {
        if (std::fabs(mat1(i, j) - mat2(i, j)) > tolerance)
        {
          allEntriesMatch = false;
          nonMatchingEntries << "(" << i << "," << j << "): " << mat1(i, j) << " vs. " << mat2(i, j)
                             << std::endl;
        }
      }
    }

    if (allEntriesMatch)
      return ::testing::AssertionSuccess();
    else
    {
      return ::testing::AssertionFailure()
             << "The following entries differ: absolute difference is not within tolerance "
             << tolerance << "\n"
             << "(i,j): " << mat1Expr << " vs. " << mat2Expr << ":" << std::endl
             << nonMatchingEntries.str();
    }
  }

  /**
   * Compare two LINALG::SerialDenseMatrix objects for double equality up to a tolerance. The
   * signature is mandated by GoogleTest's EXPECT_PRED_FORMAT3 macro.
   */
  inline ::testing::AssertionResult AssertNear(const char* mat1Expr,
      const char* mat2Expr,  // NOLINT
      const char* toleranceExpr, const LINALG::SerialDenseMatrix& mat1,
      const LINALG::SerialDenseMatrix& mat2, double tolerance)
  {
    // argument is required for the EXPECT_PRED_FORMAT3 macro of GoogleTest for pretty printing
    (void)toleranceExpr;

    const bool dimensionsMatch = mat1.RowDim() == mat2.RowDim() and mat1.ColDim() == mat2.ColDim();
    if (!dimensionsMatch)
    {
      return ::testing::AssertionFailure()
             << "dimension mismatch: " << mat1Expr << " has dimension " << mat1.RowDim() << "x"
             << mat1.ColDim() << " but " << mat2Expr << " has dimension " << mat2.RowDim() << "x"
             << mat2.ColDim() << std::endl;
    }

    bool allEntriesMatch = true;

    std::stringstream nonMatchingEntries;
    nonMatchingEntries << std::fixed << std::setprecision(PrecisionForPrinting(tolerance));
    for (int i = 0; i < mat1.RowDim(); ++i)
    {
      for (int j = 0; j < mat1.ColDim(); ++j)
      {
        if (std::fabs(mat1(i, j) - mat2(i, j)) > tolerance)
        {
          allEntriesMatch = false;
          nonMatchingEntries << "(" << i << "," << j << "): " << mat1(i, j) << " vs. " << mat2(i, j)
                             << std::endl;
        }
      }
    }

    if (allEntriesMatch)
      return ::testing::AssertionSuccess();
    else
    {
      return ::testing::AssertionFailure()
             << "The following entries differ: absolute difference is not within tolerance "
             << tolerance << "\n"
             << "(i,j): " << mat1Expr << " vs. " << mat2Expr << ":" << std::endl
             << nonMatchingEntries.str();
    }
  }
}  // namespace TESTING::INTERNAL

/**
 * Custom assertion to test for equality up to a tolerance
 *
 * Overloads are provided for
 * - std::vector<T>
 * - LINALG::Matrix
 * - LINALG::SerialDenseMatrix
 *
 * @note This and similar macros are defined to avoid writing asserts in the unexpressive
 * EXPECT_PRED_FORMATn syntax by gtest. They are all prefixed with `BACI_` to easily distinguish
 * them from gtest asserts.
 */
#define BACI_EXPECT_NEAR(actual, expected, tolerance) \
  EXPECT_PRED_FORMAT3(TESTING::INTERNAL::AssertNear, actual, expected, tolerance)

/*!
 * Extension of EXPECT_THROW which also checks for a substring in the what() expression.
 */
#define BACI_EXPECT_THROW_WITH_MESSAGE(statement, expectedException, messageSubString)  \
  std::function<void(void)> find_the_statement_below(                                   \
      [&]()                                                                             \
      {                                                                                 \
        try                                                                             \
        {                                                                               \
          statement;                                                                    \
        }                                                                               \
        catch (const expectedException& caughtException)                                \
        {                                                                               \
          using ::testing::HasSubstr;                                                   \
          EXPECT_THAT(caughtException.what(), HasSubstr(messageSubString))              \
              << "Caught the expected exception type but message has wrong substring."; \
          throw;                                                                        \
        }                                                                               \
      });                                                                               \
  EXPECT_THROW(find_the_statement_below(), expectedException) << "statement: " << #statement;

#endif  // UNITTESTS_COMMON_ASSERTIONS_H
