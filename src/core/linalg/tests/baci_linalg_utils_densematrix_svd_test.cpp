/*----------------------------------------------------------------------*/
/*! \file

\brief Unit tests for linalg SVD routines.

\level 0

*----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "baci_linalg_utils_densematrix_svd.hpp"

#include "baci_linalg_utils_densematrix_multiply.hpp"
#include "baci_unittest_utils_assertions_test.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  /*
   * \note The values for the matrix used in tests below are generated with python/numpy
   */

  template <unsigned int N>
  void AssertIsUnitaryMatrix(const CORE::LINALG::Matrix<N, N>& M)
  {
    CORE::LINALG::Matrix<N, N> MHM(false);
    MHM.MultiplyTN(M, M);

    for (unsigned int i = 0; i < N; ++i)
      for (unsigned int j = 0; j < N; ++j) EXPECT_NEAR(MHM(i, j), i == j, 1e-9);
  }

  void AssertIsUnitaryMatrix(const CORE::LINALG::SerialDenseMatrix& M)
  {
    CORE::LINALG::SerialDenseMatrix MHM(M.numRows(), M.numCols(), false);

    CORE::LINALG::multiplyTN(MHM, M, M);

    for (int i = 0; i < M.numRows(); ++i)
      for (int j = 0; j < M.numCols(); ++j) EXPECT_NEAR(MHM(i, j), i == j, 1e-9);
  }

  template <unsigned int rows, unsigned int cols, size_t length>
  void AssertSVDResult(const CORE::LINALG::Matrix<rows, cols>& A,
      const CORE::LINALG::Matrix<rows, rows>& Q, const CORE::LINALG::Matrix<rows, cols>& S,
      const CORE::LINALG::Matrix<cols, cols>& VT, const std::array<double, length>& singularValues)
  {
    // check whether SVD fulfills: A = Q * S * VT
    CORE::LINALG::Matrix<rows, cols> QS(false);
    CORE::LINALG::Matrix<rows, cols> A_result(false);
    QS.MultiplyNN(Q, S);
    A_result.MultiplyNN(QS, VT);

    FOUR_C_EXPECT_NEAR(A, A_result, 1e-9);

    // check singular values
    for (unsigned int i = 0; i < rows; ++i)
      for (unsigned int j = 0; j < cols; ++j)
      {
        if (i == j && i < singularValues.size())
          EXPECT_NEAR(S(i, j), singularValues[i], 1e-9) << "for i=" << i << ", j=" << j;
        else
          EXPECT_NEAR(S(i, j), 0.0, 1e-9) << "for i=" << i << ", j=" << j;
      }

    // check Q and VT are unitary matrices
    AssertIsUnitaryMatrix(Q);
    AssertIsUnitaryMatrix(VT);
  }

  template <size_t length>
  void AssertSVDResult(const CORE::LINALG::SerialDenseMatrix& A,
      const CORE::LINALG::SerialDenseMatrix& Q, const CORE::LINALG::SerialDenseMatrix& S,
      const CORE::LINALG::SerialDenseMatrix& VT, const std::array<double, length>& singularValues)
  {
    int rows = A.numRows();
    int cols = A.numCols();
    // check whether SVD fulfills: A = Q * S * VT
    CORE::LINALG::SerialDenseMatrix QS(rows, cols, false);
    CORE::LINALG::SerialDenseMatrix A_result(rows, cols, false);
    CORE::LINALG::multiply(QS, Q, S);
    CORE::LINALG::multiply(A_result, QS, VT);

    FOUR_C_EXPECT_NEAR(A, A_result, 1e-9);

    // check singular values
    for (int i = 0; i < rows; ++i)
      for (int j = 0; j < cols; ++j)
      {
        if (i == j && static_cast<unsigned long>(i) < singularValues.size())
          EXPECT_NEAR(S(i, j), singularValues[i], 1e-9) << "for i=" << i << ", j=" << j;
        else
          EXPECT_NEAR(S(i, j), 0.0, 1e-9) << "for i=" << i << ", j=" << j;
      }

    // check Q and VT are unitary matrices
    AssertIsUnitaryMatrix(Q);
    AssertIsUnitaryMatrix(VT);
  }


  /*
   * \note The values for the matrix used in tests below are generated with python/numpy
   *
   */
  TEST(LinalgDenseMatrixSVDTest, SVD2x2Matrix)
  {
    constexpr int rows = 2;
    constexpr int cols = 2;
    CORE::LINALG::Matrix<rows, cols, double> A;
    A(0, 0) = 0.771320643266746;
    A(0, 1) = 0.0207519493594015;
    A(1, 0) = 0.6336482349262754;
    A(1, 1) = 0.7488038825386119;

    std::array singular_values{1.1469088916371344, 0.49212144085114373};

    CORE::LINALG::Matrix<rows, rows, double> Q;
    CORE::LINALG::Matrix<rows, cols, double> S;

    CORE::LINALG::Matrix<cols, cols, double> VT;

    CORE::LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  TEST(LinalgDenseMatrixSVDTest, SVD2x2SerialDenseMatrix)
  {
    constexpr int rows = 2;
    constexpr int cols = 2;

    CORE::LINALG::SerialDenseMatrix A(rows, cols, false);
    A(0, 0) = 0.771320643266746;
    A(0, 1) = 0.0207519493594015;
    A(1, 0) = 0.6336482349262754;
    A(1, 1) = 0.7488038825386119;

    std::array singular_values{1.1469088916371344, 0.49212144085114373};

    CORE::LINALG::SerialDenseMatrix Q(rows, rows, false);
    CORE::LINALG::SerialDenseMatrix S(rows, cols, false);

    CORE::LINALG::SerialDenseMatrix VT(cols, cols, false);

    CORE::LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  TEST(LinalgDenseMatrixSVDTest, SVD3x3Matrix)
  {
    constexpr int rows = 3;
    constexpr int cols = 3;

    CORE::LINALG::Matrix<rows, cols, double> A;
    A(0, 0) = 0.4985070123025904;
    A(0, 1) = 0.22479664553084766;
    A(0, 2) = 0.19806286475962398;
    A(1, 0) = 0.7605307121989587;
    A(1, 1) = 0.16911083656253545;
    A(1, 2) = 0.08833981417401027;
    A(2, 0) = 0.6853598183677972;
    A(2, 1) = 0.9533933461949365;
    A(2, 2) = 0.003948266327914451;

    std::array singular_values{1.4366496228962886, 0.5015270327633732, 0.12760122260790782};

    CORE::LINALG::Matrix<rows, rows, double> Q;
    CORE::LINALG::Matrix<rows, cols, double> S;

    CORE::LINALG::Matrix<cols, cols, double> VT;

    CORE::LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  TEST(LinalgDenseMatrixSVDTest, SVD3x3SerialDenseMatrix)
  {
    constexpr int rows = 3;
    constexpr int cols = 3;

    CORE::LINALG::SerialDenseMatrix A(rows, cols, false);
    A(0, 0) = 0.4985070123025904;
    A(0, 1) = 0.22479664553084766;
    A(0, 2) = 0.19806286475962398;
    A(1, 0) = 0.7605307121989587;
    A(1, 1) = 0.16911083656253545;
    A(1, 2) = 0.08833981417401027;
    A(2, 0) = 0.6853598183677972;
    A(2, 1) = 0.9533933461949365;
    A(2, 2) = 0.003948266327914451;

    std::array singular_values{1.4366496228962886, 0.5015270327633732, 0.12760122260790782};

    CORE::LINALG::SerialDenseMatrix Q(rows, rows, false);
    CORE::LINALG::SerialDenseMatrix S(rows, cols, false);

    CORE::LINALG::SerialDenseMatrix VT(cols, cols, false);

    CORE::LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  TEST(LinalgDenseMatrixSVDTest, SVD4x4Matrix)
  {
    constexpr int rows = 4;
    constexpr int cols = 4;

    CORE::LINALG::Matrix<rows, cols, double> A;
    A(0, 0) = 0.5121922633857766;
    A(0, 1) = 0.8126209616521135;
    A(0, 2) = 0.6125260668293881;
    A(0, 3) = 0.7217553174317995;
    A(1, 0) = 0.29187606817063316;
    A(1, 1) = 0.9177741225129434;
    A(1, 2) = 0.7145757833976906;
    A(1, 3) = 0.5425443680112613;
    A(2, 0) = 0.14217004760152696;
    A(2, 1) = 0.3733407600514692;
    A(2, 2) = 0.6741336150663453;
    A(2, 3) = 0.4418331744229961;
    A(3, 0) = 0.4340139933332937;
    A(3, 1) = 0.6177669784693172;
    A(3, 2) = 0.5131382425543909;
    A(3, 3) = 0.6503971819314672;

    std::array singular_values{
        2.3331219940832653, 0.33317539589309747, 0.2493797215886197, 0.015235517801683926};

    CORE::LINALG::Matrix<rows, rows, double> Q;
    CORE::LINALG::Matrix<rows, cols, double> S;

    CORE::LINALG::Matrix<cols, cols, double> VT;

    CORE::LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  TEST(LinalgDenseMatrixSVDTest, SVD4x4SerialDenseMatrix)
  {
    constexpr int rows = 4;
    constexpr int cols = 4;

    CORE::LINALG::SerialDenseMatrix A(rows, cols, false);
    A(0, 0) = 0.5121922633857766;
    A(0, 1) = 0.8126209616521135;
    A(0, 2) = 0.6125260668293881;
    A(0, 3) = 0.7217553174317995;
    A(1, 0) = 0.29187606817063316;
    A(1, 1) = 0.9177741225129434;
    A(1, 2) = 0.7145757833976906;
    A(1, 3) = 0.5425443680112613;
    A(2, 0) = 0.14217004760152696;
    A(2, 1) = 0.3733407600514692;
    A(2, 2) = 0.6741336150663453;
    A(2, 3) = 0.4418331744229961;
    A(3, 0) = 0.4340139933332937;
    A(3, 1) = 0.6177669784693172;
    A(3, 2) = 0.5131382425543909;
    A(3, 3) = 0.6503971819314672;

    std::array singular_values{
        2.3331219940832653, 0.33317539589309747, 0.2493797215886197, 0.015235517801683926};

    CORE::LINALG::SerialDenseMatrix Q(rows, rows, false);
    CORE::LINALG::SerialDenseMatrix S(rows, cols, false);

    CORE::LINALG::SerialDenseMatrix VT(cols, cols, false);

    CORE::LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  TEST(LinalgDenseMatrixSVDTest, SVD3x2Matrix)
  {
    constexpr int rows = 3;
    constexpr int cols = 2;

    CORE::LINALG::Matrix<rows, cols, double> A;
    A(0, 0) = 0.6010389534045444;
    A(0, 1) = 0.8052231968327465;
    A(1, 0) = 0.5216471523936341;
    A(1, 1) = 0.9086488808086682;
    A(2, 0) = 0.3192360889885453;
    A(2, 1) = 0.09045934927090737;

    std::array singular_values{1.4710166760448906, 0.23150653036895066};

    CORE::LINALG::Matrix<rows, rows, double> Q;
    CORE::LINALG::Matrix<rows, cols, double> S;

    CORE::LINALG::Matrix<cols, cols, double> VT;

    CORE::LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  TEST(LinalgDenseMatrixSVDTest, SVD3x2SerialDenseMatrix)
  {
    constexpr int rows = 3;
    constexpr int cols = 2;

    CORE::LINALG::SerialDenseMatrix A(rows, cols, false);
    A(0, 0) = 0.6010389534045444;
    A(0, 1) = 0.8052231968327465;
    A(1, 0) = 0.5216471523936341;
    A(1, 1) = 0.9086488808086682;
    A(2, 0) = 0.3192360889885453;
    A(2, 1) = 0.09045934927090737;

    std::array singular_values{1.4710166760448906, 0.23150653036895066};

    CORE::LINALG::SerialDenseMatrix Q(rows, rows, false);
    CORE::LINALG::SerialDenseMatrix S(rows, cols, false);

    CORE::LINALG::SerialDenseMatrix VT(cols, cols, false);

    CORE::LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  TEST(LinalgDenseMatrixSVDTest, SVD2x3Matrix)
  {
    constexpr int rows = 2;
    constexpr int cols = 3;

    CORE::LINALG::Matrix<rows, cols, double> A;
    A(0, 0) = 0.30070005663620336;
    A(0, 1) = 0.11398436186354977;
    A(0, 2) = 0.8286813263076767;
    A(1, 0) = 0.04689631938924976;
    A(1, 1) = 0.6262871483113925;
    A(1, 2) = 0.5475861559192435;

    std::array singular_values{1.1329579091315047, 0.44812669032826474};

    CORE::LINALG::Matrix<rows, rows, double> Q;
    CORE::LINALG::Matrix<rows, cols, double> S;

    CORE::LINALG::Matrix<cols, cols, double> VT;

    CORE::LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  TEST(LinalgDenseMatrixSVDTest, SVD2x3SerialDenseMatrix)
  {
    constexpr int rows = 2;
    constexpr int cols = 3;

    CORE::LINALG::SerialDenseMatrix A(rows, cols, false);
    A(0, 0) = 0.30070005663620336;
    A(0, 1) = 0.11398436186354977;
    A(0, 2) = 0.8286813263076767;
    A(1, 0) = 0.04689631938924976;
    A(1, 1) = 0.6262871483113925;
    A(1, 2) = 0.5475861559192435;

    std::array singular_values{1.1329579091315047, 0.44812669032826474};

    CORE::LINALG::SerialDenseMatrix Q(rows, rows, false);
    CORE::LINALG::SerialDenseMatrix S(rows, cols, false);

    CORE::LINALG::SerialDenseMatrix VT(cols, cols, false);

    CORE::LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  TEST(LinalgDenseMatrixSVDTest, SVD1x3Matrix)
  {
    constexpr int rows = 1;
    constexpr int cols = 3;

    CORE::LINALG::Matrix<rows, cols, double> A;
    A(0, 0) = 0.8192869956700687;
    A(0, 1) = 0.1989475396788123;
    A(0, 2) = 0.8568503024577332;

    std::array singular_values{1.202083085997074};

    CORE::LINALG::Matrix<rows, rows, double> Q;
    CORE::LINALG::Matrix<rows, cols, double> S;

    CORE::LINALG::Matrix<cols, cols, double> VT;

    CORE::LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  TEST(LinalgDenseMatrixSVDTest, SVD1x3SerialDenseMatrix)
  {
    constexpr int rows = 1;
    constexpr int cols = 3;

    CORE::LINALG::SerialDenseMatrix A(rows, cols, false);
    A(0, 0) = 0.8192869956700687;
    A(0, 1) = 0.1989475396788123;
    A(0, 2) = 0.8568503024577332;

    std::array singular_values{1.202083085997074};

    CORE::LINALG::SerialDenseMatrix Q(rows, rows, false);
    CORE::LINALG::SerialDenseMatrix S(rows, cols, false);

    CORE::LINALG::SerialDenseMatrix VT(cols, cols, false);

    CORE::LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  TEST(LinalgDenseMatrixSVDTest, SVD3x1Matrix)
  {
    constexpr int rows = 3;
    constexpr int cols = 1;

    CORE::LINALG::Matrix<rows, cols, double> A;
    A(0, 0) = 0.3516526394320879;
    A(1, 0) = 0.7546476915298572;
    A(2, 0) = 0.2959617068796787;

    std::array singular_values{0.88359835281084};

    CORE::LINALG::Matrix<rows, rows, double> Q;
    CORE::LINALG::Matrix<rows, cols, double> S;

    CORE::LINALG::Matrix<cols, cols, double> VT;

    CORE::LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  TEST(LinalgDenseMatrixSVDTest, SVD3x1SerialDenseMatrix)
  {
    constexpr int rows = 3;
    constexpr int cols = 1;

    CORE::LINALG::SerialDenseMatrix A(rows, cols, false);
    A(0, 0) = 0.3516526394320879;
    A(1, 0) = 0.7546476915298572;
    A(2, 0) = 0.2959617068796787;

    std::array singular_values{0.88359835281084};

    CORE::LINALG::SerialDenseMatrix Q(rows, rows, false);
    CORE::LINALG::SerialDenseMatrix S(rows, cols, false);

    CORE::LINALG::SerialDenseMatrix VT(cols, cols, false);

    CORE::LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }
}  // namespace

FOUR_C_NAMESPACE_CLOSE
