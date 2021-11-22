/*----------------------------------------------------------------------*/
/*! \file
\brief Unit tests for linalg SVD routine.

\level 0
*/


#ifndef UNIT_LINALG_UTILS_DENSEMATRIX_SVD_H
#define UNIT_LINALG_UTILS_DENSEMATRIX_SVD_H


#include "Unittests/src/common/unit_cxx_test_wrapper.H"

#include "src/linalg/linalg_utils_densematrix_svd.H"
#include "Unittests/src/common/special_assertions.H"


// Forward declaration.
namespace LINALG
{
  class Linalg_SVD_TestSuite;
}


/**
 * \brief Class to test linalg routines.
 *
 * \note The values for the matrix used in tests below are generated with python/numpy
 */
class LINALG::Linalg_SVD_TestSuite : public BACICxxTestWrapper
{
 private:
  template <unsigned int N>
  void AssertIsUnitaryMatrix(const LINALG::Matrix<N, N>& M)
  {
    LINALG::Matrix<N, N> MHM(false);
    MHM.MultiplyTN(M, M);

    for (unsigned int i = 0; i < N; ++i)
    {
      for (unsigned int j = 0; j < N; ++j)
      {
        TS_ASSERT_DELTA(MHM(i, j), i == j, 1e-9);
      }
    }
  }

  void AssertIsUnitaryMatrix(const LINALG::SerialDenseMatrix& M)
  {
    LINALG::SerialDenseMatrix MHM(M.RowDim(), M.ColDim(), false);

    MHM.Multiply('T', 'N', 1.0, M, M, 0.0);


    for (int i = 0; i < M.RowDim(); ++i)
    {
      for (int j = 0; j < M.ColDim(); ++j)
      {
        TS_ASSERT_DELTA(MHM(i, j), i == j, 1e-9);
      }
    }
  }

  template <unsigned int rows, unsigned int cols, unsigned long rank>
  void AssertSVDResult(const LINALG::Matrix<rows, cols>& A, const LINALG::Matrix<rows, rows>& Q,
      const LINALG::Matrix<rows, cols>& S, const LINALG::Matrix<cols, cols>& VT,
      const std::array<double, rank>& singularValues)
  {
    // check whether SVD fulfills: A = Q * S * VT
    LINALG::Matrix<rows, cols> QS(false);
    LINALG::Matrix<rows, cols> A_result(false);
    QS.MultiplyNN(Q, S);
    A_result.MultiplyNN(QS, VT);

    TESTING::AssertDelta(A, A_result, 1e-9);

    // check singular values
    for (unsigned int i = 0; i < rows; ++i)
    {
      for (unsigned int j = 0; j < cols; ++j)
      {
        if (i == j && i < rank)
        {
          TS_ASSERT_DELTA(S(i, j), singularValues[i], 1e-9);
        }
        else
        {
          TS_ASSERT_DELTA(S(i, j), 0.0, 1e-9);
        }
      }
    }

    // check Q and VT are unitary matrices
    AssertIsUnitaryMatrix(Q);
    AssertIsUnitaryMatrix(VT);
  }

  template <unsigned long rank>
  void AssertSVDResult(const LINALG::SerialDenseMatrix& A, const LINALG::SerialDenseMatrix& Q,
      const LINALG::SerialDenseMatrix& S, const LINALG::SerialDenseMatrix& VT,
      const std::array<double, rank>& singularValues)
  {
    int rows = A.RowDim();
    int cols = A.ColDim();
    // check whether SVD fulfills: A = Q * S * VT
    LINALG::SerialDenseMatrix QS(rows, cols, false);
    LINALG::SerialDenseMatrix A_result(rows, cols, false);
    QS.Multiply('N', 'N', 1.0, Q, S, 0.0);
    A_result.Multiply('N', 'N', 1.0, QS, VT, 0.0);

    TESTING::AssertDelta(A, A_result, 1e-9);

    // check singular values
    for (int i = 0; i < rows; ++i)
    {
      for (int j = 0; j < cols; ++j)
      {
        if (i == j && static_cast<unsigned long>(i) < rank)
        {
          TS_ASSERT_DELTA(S(i, j), singularValues[i], 1e-9);
        }
        else
        {
          TS_ASSERT_DELTA(S(i, j), 0.0, 1e-9);
        }
      }
    }

    // check Q and VT are unitary matrices
    AssertIsUnitaryMatrix(Q);
    AssertIsUnitaryMatrix(VT);
  }

 public:
  void TestSVD2x2Matrix()
  {
    constexpr int rows = 2;
    constexpr int cols = 2;

    LINALG::Matrix<rows, cols, double> A;
    A(0, 0) = 0.771320643266746;
    A(0, 1) = 0.0207519493594015;
    A(1, 0) = 0.6336482349262754;
    A(1, 1) = 0.7488038825386119;

    std::array<double, 2> singular_values = {1.1469088916371344, 0.49212144085114373};
    LINALG::Matrix<rows, rows, double> Q;
    LINALG::Matrix<rows, cols, double> S;

    LINALG::Matrix<cols, cols, double> VT;

    LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  void TestSVD2x2SerialDenseMatrix()
  {
    constexpr int rows = 2;
    constexpr int cols = 2;

    LINALG::SerialDenseMatrix A(rows, cols, false);
    A(0, 0) = 0.771320643266746;
    A(0, 1) = 0.0207519493594015;
    A(1, 0) = 0.6336482349262754;
    A(1, 1) = 0.7488038825386119;

    std::array<double, 2> singular_values = {1.1469088916371344, 0.49212144085114373};
    LINALG::SerialDenseMatrix Q(rows, rows, false);
    LINALG::SerialDenseMatrix S(rows, cols, false);

    LINALG::SerialDenseMatrix VT(cols, cols, false);

    LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  void TestSVD3x3Matrix()
  {
    constexpr int rows = 3;
    constexpr int cols = 3;

    LINALG::Matrix<rows, cols, double> A;
    A(0, 0) = 0.4985070123025904;
    A(0, 1) = 0.22479664553084766;
    A(0, 2) = 0.19806286475962398;
    A(1, 0) = 0.7605307121989587;
    A(1, 1) = 0.16911083656253545;
    A(1, 2) = 0.08833981417401027;
    A(2, 0) = 0.6853598183677972;
    A(2, 1) = 0.9533933461949365;
    A(2, 2) = 0.003948266327914451;

    std::array<double, 3> singular_values = {
        1.4366496228962886, 0.5015270327633732, 0.12760122260790782};
    LINALG::Matrix<rows, rows, double> Q;
    LINALG::Matrix<rows, cols, double> S;

    LINALG::Matrix<cols, cols, double> VT;

    LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  void TestSVD3x3SerialDenseMatrix()
  {
    constexpr int rows = 3;
    constexpr int cols = 3;

    LINALG::SerialDenseMatrix A(rows, cols, false);
    A(0, 0) = 0.4985070123025904;
    A(0, 1) = 0.22479664553084766;
    A(0, 2) = 0.19806286475962398;
    A(1, 0) = 0.7605307121989587;
    A(1, 1) = 0.16911083656253545;
    A(1, 2) = 0.08833981417401027;
    A(2, 0) = 0.6853598183677972;
    A(2, 1) = 0.9533933461949365;
    A(2, 2) = 0.003948266327914451;

    std::array<double, 3> singular_values = {
        1.4366496228962886, 0.5015270327633732, 0.12760122260790782};
    LINALG::SerialDenseMatrix Q(rows, rows, false);
    LINALG::SerialDenseMatrix S(rows, cols, false);

    LINALG::SerialDenseMatrix VT(cols, cols, false);

    LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  void TestSVD4x4Matrix()
  {
    constexpr int rows = 4;
    constexpr int cols = 4;

    LINALG::Matrix<rows, cols, double> A;
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

    std::array<double, 4> singular_values = {
        2.3331219940832653, 0.33317539589309747, 0.2493797215886197, 0.015235517801683926};
    LINALG::Matrix<rows, rows, double> Q;
    LINALG::Matrix<rows, cols, double> S;

    LINALG::Matrix<cols, cols, double> VT;

    LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  void TestSVD4x4SerialDenseMatrix()
  {
    constexpr int rows = 4;
    constexpr int cols = 4;

    LINALG::SerialDenseMatrix A(rows, cols, false);
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

    std::array<double, 4> singular_values = {
        2.3331219940832653, 0.33317539589309747, 0.2493797215886197, 0.015235517801683926};
    LINALG::SerialDenseMatrix Q(rows, rows, false);
    LINALG::SerialDenseMatrix S(rows, cols, false);

    LINALG::SerialDenseMatrix VT(cols, cols, false);

    LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  void TestSVD3x2Matrix()
  {
    constexpr int rows = 3;
    constexpr int cols = 2;

    LINALG::Matrix<rows, cols, double> A;
    A(0, 0) = 0.6010389534045444;
    A(0, 1) = 0.8052231968327465;
    A(1, 0) = 0.5216471523936341;
    A(1, 1) = 0.9086488808086682;
    A(2, 0) = 0.3192360889885453;
    A(2, 1) = 0.09045934927090737;

    std::array<double, 2> singular_values = {1.4710166760448906, 0.23150653036895066};
    LINALG::Matrix<rows, rows, double> Q;
    LINALG::Matrix<rows, cols, double> S;

    LINALG::Matrix<cols, cols, double> VT;

    LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  void TestSVD3x2SerialDenseMatrix()
  {
    constexpr int rows = 3;
    constexpr int cols = 2;

    LINALG::SerialDenseMatrix A(rows, cols, false);
    A(0, 0) = 0.6010389534045444;
    A(0, 1) = 0.8052231968327465;
    A(1, 0) = 0.5216471523936341;
    A(1, 1) = 0.9086488808086682;
    A(2, 0) = 0.3192360889885453;
    A(2, 1) = 0.09045934927090737;

    std::array<double, 2> singular_values = {1.4710166760448906, 0.23150653036895066};
    LINALG::SerialDenseMatrix Q(rows, rows, false);
    LINALG::SerialDenseMatrix S(rows, cols, false);

    LINALG::SerialDenseMatrix VT(cols, cols, false);

    LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  void TestSVD2x3Matrix()
  {
    constexpr int rows = 2;
    constexpr int cols = 3;

    LINALG::Matrix<rows, cols, double> A;
    A(0, 0) = 0.30070005663620336;
    A(0, 1) = 0.11398436186354977;
    A(0, 2) = 0.8286813263076767;
    A(1, 0) = 0.04689631938924976;
    A(1, 1) = 0.6262871483113925;
    A(1, 2) = 0.5475861559192435;

    std::array<double, 2> singular_values = {1.1329579091315047, 0.44812669032826474};
    LINALG::Matrix<rows, rows, double> Q;
    LINALG::Matrix<rows, cols, double> S;

    LINALG::Matrix<cols, cols, double> VT;

    LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  void TestSVD2x3SerialDenseMatrix()
  {
    constexpr int rows = 2;
    constexpr int cols = 3;

    LINALG::SerialDenseMatrix A(rows, cols, false);
    A(0, 0) = 0.30070005663620336;
    A(0, 1) = 0.11398436186354977;
    A(0, 2) = 0.8286813263076767;
    A(1, 0) = 0.04689631938924976;
    A(1, 1) = 0.6262871483113925;
    A(1, 2) = 0.5475861559192435;

    std::array<double, 2> singular_values = {1.1329579091315047, 0.44812669032826474};
    LINALG::SerialDenseMatrix Q(rows, rows, false);
    LINALG::SerialDenseMatrix S(rows, cols, false);

    LINALG::SerialDenseMatrix VT(cols, cols, false);

    LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  void TestSVD1x3Matrix()
  {
    constexpr int rows = 1;
    constexpr int cols = 3;

    LINALG::Matrix<rows, cols, double> A;
    A(0, 0) = 0.8192869956700687;
    A(0, 1) = 0.1989475396788123;
    A(0, 2) = 0.8568503024577332;

    std::array<double, 1> singular_values = {1.202083085997074};
    LINALG::Matrix<rows, rows, double> Q;
    LINALG::Matrix<rows, cols, double> S;

    LINALG::Matrix<cols, cols, double> VT;

    LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  void TestSVD1x3SerialDenseMatrix()
  {
    constexpr int rows = 1;
    constexpr int cols = 3;

    LINALG::SerialDenseMatrix A(rows, cols, false);
    A(0, 0) = 0.8192869956700687;
    A(0, 1) = 0.1989475396788123;
    A(0, 2) = 0.8568503024577332;

    std::array<double, 1> singular_values = {1.202083085997074};
    LINALG::SerialDenseMatrix Q(rows, rows, false);
    LINALG::SerialDenseMatrix S(rows, cols, false);

    LINALG::SerialDenseMatrix VT(cols, cols, false);

    LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  void TestSVD3x1Matrix()
  {
    constexpr int rows = 3;
    constexpr int cols = 1;

    LINALG::Matrix<rows, cols, double> A;
    A(0, 0) = 0.3516526394320879;
    A(1, 0) = 0.7546476915298572;
    A(2, 0) = 0.2959617068796787;

    std::array<double, 1> singular_values = {0.88359835281084};
    LINALG::Matrix<rows, rows, double> Q;
    LINALG::Matrix<rows, cols, double> S;

    LINALG::Matrix<cols, cols, double> VT;

    LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }

  void TestSVD3x1SerialDenseMatrix()
  {
    constexpr int rows = 3;
    constexpr int cols = 1;

    LINALG::SerialDenseMatrix A(rows, cols, false);
    A(0, 0) = 0.3516526394320879;
    A(1, 0) = 0.7546476915298572;
    A(2, 0) = 0.2959617068796787;

    std::array<double, 1> singular_values = {0.88359835281084};
    LINALG::SerialDenseMatrix Q(rows, rows, false);
    LINALG::SerialDenseMatrix S(rows, cols, false);

    LINALG::SerialDenseMatrix VT(cols, cols, false);

    LINALG::SVD(A, Q, S, VT);

    AssertSVDResult(A, Q, S, VT, singular_values);
  }
};

#endif
