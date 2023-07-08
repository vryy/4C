/*----------------------------------------------------------------------*/
/*! \file

\brief Unit tests for linalg dense eigen routines.

\level 0

*----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "linalg_utils_densematrix_eigen.H"

#include "unittest_utils_assertions.h"

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>

namespace
{
  /*
   * \note Helper functions for the dense eigen test routines.
   *
   */
  template <size_t length>
  void AssertEigenValues(
      const Epetra_SerialDenseVector& eigenvalues, const std::array<double, length>& eig_compare)
  {
    EXPECT_EQ(eigenvalues.Length(), length);
    for (unsigned i = 0; i < length; ++i) EXPECT_NEAR(eigenvalues[i], eig_compare[i], 1e-9);
  }

  template <unsigned int size, size_t length>
  void AssertEigenValues(
      const LINALG::Matrix<size, size>& eigenvalues, const std::array<double, length>& eig_compare)
  {
    for (unsigned i = 0; i < size; ++i)
    {
      for (unsigned j = 0; j < size; ++j)
      {
        if (i == j)
          EXPECT_NEAR(eigenvalues(i, j), eig_compare[i], 1e-9) << "for i=" << i << ", j=" << j;
        else
          EXPECT_NEAR(eigenvalues(i, j), 0.0, 1e-9) << "for i=" << i << ", j=" << j;
      }
    }
  }

  template <size_t length>
  void AssertEigenProblem(const Epetra_SerialDenseMatrix& A,
      const Epetra_SerialDenseVector& eigenvalues, const Epetra_SerialDenseMatrix& eigenvectors,
      const std::array<double, length>& eig_compare)
  {
    EXPECT_EQ(eigenvalues.Length(), length);
    EXPECT_EQ(A.RowDim(), length);
    EXPECT_EQ(A.ColDim(), length);
    EXPECT_EQ(eigenvectors.RowDim(), length);
    EXPECT_EQ(eigenvectors.ColDim(), length);

    AssertEigenValues(eigenvalues, eig_compare);

    Epetra_SerialDenseMatrix A_result(length, length, true);
    for (std::size_t i = 0; i < length; ++i)
    {
      Epetra_SerialDenseMatrix v(length, 1, false);
      for (std::size_t j = 0; j < length; ++j) v(j, 0) = eigenvectors(j, i);
      A_result.Multiply('N', 'T', eigenvalues(i), v, v, 1.0);
    }

    BACI_EXPECT_NEAR(A, A_result, 1e-9);
  }

  template <unsigned int size, size_t length>
  void AssertEigenProblem(const LINALG::Matrix<size, size>& A,
      const LINALG::Matrix<size, size>& eigenvalues, const LINALG::Matrix<size, size>& eigenvectors,
      const std::array<double, length>& eig_compare)
  {
    AssertEigenValues(eigenvalues, eig_compare);

    LINALG::Matrix<size, size> A_result(true);
    for (unsigned int i = 0; i < size; ++i)
    {
      LINALG::Matrix<size, 1> v(false);
      for (unsigned int j = 0; j < size; ++j) v(j, 0) = eigenvectors(j, i);
      A_result.MultiplyNT(eigenvalues(i, i), v, v, 1.0);
    }

    BACI_EXPECT_NEAR(A, A_result, 1e-9);
  }


  /*
   * \note The values for the matrix used in tests below are generated with python/numpy
   */
  TEST(LinalgDenseMatrixEigenTest, 2x2SymmetricEigenValues)
  {
    Epetra_SerialDenseMatrix A(2, 2, false);
    A(0, 0) = 0.9964456203546112;
    A(0, 1) = 0.490484665405466;
    A(1, 0) = 0.490484665405466;
    A(1, 1) = 0.5611378979071144;

    std::array eigenvalues{0.24218351254540577, 1.3154000057163198};

    Epetra_SerialDenseVector L(2);
    LINALG::SymmetricEigenValues(A, L, false);

    AssertEigenValues(L, eigenvalues);
  }

  TEST(LinalgDenseMatrixEigenTest, 2x2SymmetricEigenProblem)
  {
    Epetra_SerialDenseMatrix A(2, 2, false);
    A(0, 0) = 0.9964456203546112;
    A(0, 1) = 0.490484665405466;
    A(1, 0) = 0.490484665405466;
    A(1, 1) = 0.5611378979071144;

    std::array eigenvalues{0.24218351254540577, 1.3154000057163198};

    Epetra_SerialDenseMatrix eigenvectors(A);
    Epetra_SerialDenseVector L(2);
    LINALG::SymmetricEigenProblem(eigenvectors, L, false);

    AssertEigenProblem(A, L, eigenvectors, eigenvalues);
  }

  TEST(LinalgDenseMatrixEigenTest, 2x2SymmetricEigenNoVectors)
  {
    Epetra_SerialDenseMatrix A(2, 2, false);
    A(0, 0) = 0.9964456203546112;
    A(0, 1) = 0.490484665405466;
    A(1, 0) = 0.490484665405466;
    A(1, 1) = 0.5611378979071144;

    std::array eigenvalues{0.24218351254540577, 1.3154000057163198};

    Epetra_SerialDenseVector L(2);
    LINALG::SymmetricEigen(A, L, 'N', false);

    AssertEigenValues(L, eigenvalues);
  }

  TEST(LinalgDenseMatrixEigenTest, 2x2SymmetricEigenVectors)
  {
    Epetra_SerialDenseMatrix A(2, 2, false);
    A(0, 0) = 0.9964456203546112;
    A(0, 1) = 0.490484665405466;
    A(1, 0) = 0.490484665405466;
    A(1, 1) = 0.5611378979071144;

    std::array eigenvalues{0.24218351254540577, 1.3154000057163198};

    Epetra_SerialDenseMatrix eigenvectors(A);
    Epetra_SerialDenseVector L(2);
    LINALG::SymmetricEigen(eigenvectors, L, 'V', false);

    AssertEigenProblem(A, L, eigenvectors, eigenvalues);
  }

  TEST(LinalgDenseMatrixEigenTest, 2x2SYEV)
  {
    LINALG::Matrix<2, 2> A(false);
    A(0, 0) = 0.9964456203546112;
    A(0, 1) = 0.490484665405466;
    A(1, 0) = 0.490484665405466;
    A(1, 1) = 0.5611378979071144;

    std::array eigenvalues{0.24218351254540577, 1.3154000057163198};

    LINALG::Matrix<2, 2> V(false);
    LINALG::Matrix<2, 2> S(false);
    LINALG::SYEV(A, S, V);

    AssertEigenProblem(A, S, V, eigenvalues);
  }

  TEST(LinalgDenseMatrixEigenTest, 3x3SymmetricEigenValues)
  {
    Epetra_SerialDenseMatrix A(3, 3, false);
    A(0, 0) = 1.2966342861458506;
    A(0, 1) = 0.8940941796919223;
    A(0, 2) = 0.16862685184206302;
    A(1, 0) = 0.8940941796919223;
    A(1, 1) = 0.9880908794535803;
    A(1, 2) = 0.06322733832497837;
    A(2, 0) = 0.16862685184206302;
    A(2, 1) = 0.06322733832497837;
    A(2, 2) = 0.047048409972083906;

    std::array eigenvalues{0.01628207201103285, 0.2515293645924337, 2.0639621389680487};

    Epetra_SerialDenseVector L(3);
    LINALG::SymmetricEigenValues(A, L, false);

    AssertEigenValues(L, eigenvalues);
  }

  TEST(LinalgDenseMatrixEigenTest, 3x3SymmetricEigenProblem)
  {
    Epetra_SerialDenseMatrix A(3, 3, false);
    A(0, 0) = 1.2966342861458506;
    A(0, 1) = 0.8940941796919223;
    A(0, 2) = 0.16862685184206302;
    A(1, 0) = 0.8940941796919223;
    A(1, 1) = 0.9880908794535803;
    A(1, 2) = 0.06322733832497837;
    A(2, 0) = 0.16862685184206302;
    A(2, 1) = 0.06322733832497837;
    A(2, 2) = 0.047048409972083906;

    std::array eigenvalues{0.01628207201103285, 0.2515293645924337, 2.0639621389680487};

    Epetra_SerialDenseMatrix V(A);
    Epetra_SerialDenseVector L(3);
    LINALG::SymmetricEigenProblem(V, L, false);

    AssertEigenProblem(A, L, V, eigenvalues);
  }

  TEST(LinalgDenseMatrixEigenTest, 3x3SymmetricEigenNoVectors)
  {
    Epetra_SerialDenseMatrix A(3, 3, false);
    A(0, 0) = 1.2966342861458506;
    A(0, 1) = 0.8940941796919223;
    A(0, 2) = 0.16862685184206302;
    A(1, 0) = 0.8940941796919223;
    A(1, 1) = 0.9880908794535803;
    A(1, 2) = 0.06322733832497837;
    A(2, 0) = 0.16862685184206302;
    A(2, 1) = 0.06322733832497837;
    A(2, 2) = 0.047048409972083906;

    std::array eigenvalues{0.01628207201103285, 0.2515293645924337, 2.0639621389680487};

    Epetra_SerialDenseVector L(3);
    LINALG::SymmetricEigen(A, L, 'N', false);

    AssertEigenValues(L, eigenvalues);
  }

  TEST(LinalgDenseMatrixEigenTest, 3x3SymmetricEigenVectors)
  {
    Epetra_SerialDenseMatrix A(3, 3, false);
    A(0, 0) = 1.2966342861458506;
    A(0, 1) = 0.8940941796919223;
    A(0, 2) = 0.16862685184206302;
    A(1, 0) = 0.8940941796919223;
    A(1, 1) = 0.9880908794535803;
    A(1, 2) = 0.06322733832497837;
    A(2, 0) = 0.16862685184206302;
    A(2, 1) = 0.06322733832497837;
    A(2, 2) = 0.047048409972083906;

    std::array eigenvalues{0.01628207201103285, 0.2515293645924337, 2.0639621389680487};

    Epetra_SerialDenseMatrix V(A);
    Epetra_SerialDenseVector L(3);
    LINALG::SymmetricEigen(V, L, 'V', false);

    AssertEigenProblem(A, L, V, eigenvalues);
  }

  TEST(LinalgDenseMatrixEigenTest, 3x3SYEV)
  {
    LINALG::Matrix<3, 3> A(false);
    A(0, 0) = 1.2966342861458506;
    A(0, 1) = 0.8940941796919223;
    A(0, 2) = 0.16862685184206302;
    A(1, 0) = 0.8940941796919223;
    A(1, 1) = 0.9880908794535803;
    A(1, 2) = 0.06322733832497837;
    A(2, 0) = 0.16862685184206302;
    A(2, 1) = 0.06322733832497837;
    A(2, 2) = 0.047048409972083906;

    std::array eigenvalues{0.01628207201103285, 0.2515293645924337, 2.0639621389680487};

    LINALG::Matrix<3, 3> V(false);
    LINALG::Matrix<3, 3> S(false);
    LINALG::SYEV(A, S, V);

    AssertEigenProblem(A, S, V, eigenvalues);
  }

  TEST(LinalgDenseMatrixEigenTest, 4x4SymmetricEigenValues)
  {
    Epetra_SerialDenseMatrix A(4, 4, false);
    A(0, 0) = 0.5561130226871257;
    A(0, 1) = 1.0052918588741722;
    A(0, 2) = 0.8408494685470309;
    A(0, 3) = 0.8731301282118089;
    A(1, 0) = 1.0052918588741722;
    A(1, 1) = 2.023681530073728;
    A(1, 2) = 1.7222521019056944;
    A(1, 3) = 1.6511949164466262;
    A(2, 0) = 0.8408494685470309;
    A(2, 1) = 1.7222521019056944;
    A(2, 2) = 1.6035737196981317;
    A(2, 3) = 1.4613812746280035;
    A(3, 0) = 0.8731301282118089;
    A(3, 1) = 1.6511949164466262;
    A(3, 2) = 1.4613812746280035;
    A(3, 3) = 1.4335181777869124;

    std::array eigenvalues{
        0.00023212100268553735, 0.06219024553961773, 0.11100584442852221, 5.443458239275074};

    Epetra_SerialDenseVector L(4);
    LINALG::SymmetricEigenValues(A, L, false);

    AssertEigenValues(L, eigenvalues);
  }

  TEST(LinalgDenseMatrixEigenTest, 4x4SymmetricEigenProblem)
  {
    Epetra_SerialDenseMatrix A(4, 4, false);
    A(0, 0) = 0.5561130226871257;
    A(0, 1) = 1.0052918588741722;
    A(0, 2) = 0.8408494685470309;
    A(0, 3) = 0.8731301282118089;
    A(1, 0) = 1.0052918588741722;
    A(1, 1) = 2.023681530073728;
    A(1, 2) = 1.7222521019056944;
    A(1, 3) = 1.6511949164466262;
    A(2, 0) = 0.8408494685470309;
    A(2, 1) = 1.7222521019056944;
    A(2, 2) = 1.6035737196981317;
    A(2, 3) = 1.4613812746280035;
    A(3, 0) = 0.8731301282118089;
    A(3, 1) = 1.6511949164466262;
    A(3, 2) = 1.4613812746280035;
    A(3, 3) = 1.4335181777869124;

    std::array eigenvalues{
        0.00023212100268553735, 0.06219024553961773, 0.11100584442852221, 5.443458239275074};

    Epetra_SerialDenseMatrix V(A);
    Epetra_SerialDenseVector L(4);
    LINALG::SymmetricEigenProblem(V, L, false);

    AssertEigenProblem(A, L, V, eigenvalues);
  }

  TEST(LinalgDenseMatrixEigenTest, 4x4SymmetricEigenNoVectors)
  {
    Epetra_SerialDenseMatrix A(4, 4, false);
    A(0, 0) = 0.5561130226871257;
    A(0, 1) = 1.0052918588741722;
    A(0, 2) = 0.8408494685470309;
    A(0, 3) = 0.8731301282118089;
    A(1, 0) = 1.0052918588741722;
    A(1, 1) = 2.023681530073728;
    A(1, 2) = 1.7222521019056944;
    A(1, 3) = 1.6511949164466262;
    A(2, 0) = 0.8408494685470309;
    A(2, 1) = 1.7222521019056944;
    A(2, 2) = 1.6035737196981317;
    A(2, 3) = 1.4613812746280035;
    A(3, 0) = 0.8731301282118089;
    A(3, 1) = 1.6511949164466262;
    A(3, 2) = 1.4613812746280035;
    A(3, 3) = 1.4335181777869124;

    std::array eigenvalues{
        0.00023212100268553735, 0.06219024553961773, 0.11100584442852221, 5.443458239275074};

    Epetra_SerialDenseVector L(4);
    LINALG::SymmetricEigen(A, L, 'N', false);

    AssertEigenValues(L, eigenvalues);
  }

  TEST(LinalgDenseMatrixEigenTest, 4x4SymmetricEigenVectors)
  {
    Epetra_SerialDenseMatrix A(4, 4, false);
    A(0, 0) = 0.5561130226871257;
    A(0, 1) = 1.0052918588741722;
    A(0, 2) = 0.8408494685470309;
    A(0, 3) = 0.8731301282118089;
    A(1, 0) = 1.0052918588741722;
    A(1, 1) = 2.023681530073728;
    A(1, 2) = 1.7222521019056944;
    A(1, 3) = 1.6511949164466262;
    A(2, 0) = 0.8408494685470309;
    A(2, 1) = 1.7222521019056944;
    A(2, 2) = 1.6035737196981317;
    A(2, 3) = 1.4613812746280035;
    A(3, 0) = 0.8731301282118089;
    A(3, 1) = 1.6511949164466262;
    A(3, 2) = 1.4613812746280035;
    A(3, 3) = 1.4335181777869124;

    std::array eigenvalues{
        0.00023212100268553735, 0.06219024553961773, 0.11100584442852221, 5.443458239275074};

    Epetra_SerialDenseMatrix V(A);
    Epetra_SerialDenseVector L(4);
    LINALG::SymmetricEigen(V, L, 'V', false);

    AssertEigenProblem(A, L, V, eigenvalues);
  }

  TEST(LinalgDenseMatrixEigenTest, 4x4SYEV)
  {
    LINALG::Matrix<4, 4> A(false);
    A(0, 0) = 0.5561130226871257;
    A(0, 1) = 1.0052918588741722;
    A(0, 2) = 0.8408494685470309;
    A(0, 3) = 0.8731301282118089;
    A(1, 0) = 1.0052918588741722;
    A(1, 1) = 2.023681530073728;
    A(1, 2) = 1.7222521019056944;
    A(1, 3) = 1.6511949164466262;
    A(2, 0) = 0.8408494685470309;
    A(2, 1) = 1.7222521019056944;
    A(2, 2) = 1.6035737196981317;
    A(2, 3) = 1.4613812746280035;
    A(3, 0) = 0.8731301282118089;
    A(3, 1) = 1.6511949164466262;
    A(3, 2) = 1.4613812746280035;
    A(3, 3) = 1.4335181777869124;

    std::array eigenvalues{
        0.00023212100268553735, 0.06219024553961773, 0.11100584442852221, 5.443458239275074};

    LINALG::Matrix<4, 4> V(false);
    LINALG::Matrix<4, 4> S(false);
    LINALG::SYEV(A, S, V);

    AssertEigenProblem(A, S, V, eigenvalues);
  }
}  // namespace
