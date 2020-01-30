/*----------------------------------------------------------------------*/
/*! \file
\brief Unit tests for linalg routines.

\level 0
\maintainer Sebastian Fuchs
*/


#ifndef UNIT_LINALG_UTILS_DENSEMATRIX_DETERMINANT_H
#define UNIT_LINALG_UTILS_DENSEMATRIX_DETERMINANT_H


#include <cxxtest/TestSuite.h>

#include "src/linalg/linalg_utils_densematrix_determinant.H"


// Forward declaration.
namespace LINALG
{
  class Linalg_Determinant_TestSuite;
}


/**
 * \brief Class to test linalg routines.
 *
 * \note The values for the matrix used in tests below are generated with Mathematica:
 *       > SeedRandom[666];
 *       > A = Table[RandomReal[WorkingPrecision->50], {i, n}, {j, n}];
 *       where n needs do be replace by the dimension, e.g., n=2, n=3, or n=4
 */
class LINALG::Linalg_Determinant_TestSuite : public CxxTest::TestSuite
{
 public:
  void test_Determinant2x2()
  {
    LINALG::Matrix<2, 2, double> A;
    A(0, 0) = 0.72903241936703114203;
    A(1, 0) = 0.81862230026150939335;
    A(0, 1) = 0.32707405507901372465;
    A(1, 1) = 0.0052737129228371719370;

    TS_ASSERT_DELTA(LINALG::Determinant2x2(A), -0.2639054076334587, 1e-14);
  }

  void test_Determinant3x3()
  {
    LINALG::Matrix<3, 3, double> A;
    A(0, 0) = 0.72903241936703114203;
    A(1, 0) = 0.0052737129228371719370;
    A(2, 0) = 0.36847164343389089096;
    A(0, 1) = 0.32707405507901372465;
    A(1, 1) = 0.87570663114228933311;
    A(2, 1) = 0.76895132151127114661;
    A(0, 2) = 0.81862230026150939335;
    A(1, 2) = 0.64019842179333806573;
    A(2, 2) = 0.69378923027976465858;

    TS_ASSERT_DELTA(LINALG::Determinant3x3(A), -0.1008304741716571, 1e-14);
  }

  void test_Determinant4x4()
  {
    LINALG::Matrix<4, 4, double> A;
    A(0, 0) = 0.72903241936703114203;
    A(1, 0) = 0.87570663114228933311;
    A(2, 0) = 0.69378923027976465858;
    A(3, 0) = 0.019637190415090362652;
    A(0, 1) = 0.32707405507901372465;
    A(1, 1) = 0.64019842179333806573;
    A(2, 1) = 0.15928293569477706215;
    A(3, 1) = 0.13119201434024140151;
    A(0, 2) = 0.81862230026150939335;
    A(1, 2) = 0.36847164343389089096;
    A(2, 2) = 0.12278929762839221138;
    A(3, 2) = 0.12028240083390837511;
    A(0, 3) = 0.0052737129228371719370;
    A(1, 3) = 0.76895132151127114661;
    A(2, 3) = 0.024003735765356129168;
    A(3, 3) = 0.27465069811053651449;

    TS_ASSERT_DELTA(LINALG::Determinant4x4(A), -0.01620776397174742, 1e-14);
  }

  void test_DeterminantLU_2x2()
  {
    LINALG::SerialDenseMatrix A(2, 2, true);
    A(0, 0) = 0.72903241936703114203;
    A(1, 0) = 0.81862230026150939335;
    A(0, 1) = 0.32707405507901372465;
    A(1, 1) = 0.0052737129228371719370;

    TS_ASSERT_DELTA(LINALG::DeterminantLU(A), -0.2639054076334587, 1e-14);
  }

  void test_DeterminantLU_3x3()
  {
    LINALG::SerialDenseMatrix A(3, 3, true);
    A(0, 0) = 0.72903241936703114203;
    A(1, 0) = 0.0052737129228371719370;
    A(2, 0) = 0.36847164343389089096;
    A(0, 1) = 0.32707405507901372465;
    A(1, 1) = 0.87570663114228933311;
    A(2, 1) = 0.76895132151127114661;
    A(0, 2) = 0.81862230026150939335;
    A(1, 2) = 0.64019842179333806573;
    A(2, 2) = 0.69378923027976465858;

    TS_ASSERT_DELTA(LINALG::DeterminantLU(A), -0.1008304741716571, 1e-14);
  }

  void test_DeterminantLU_4x4()
  {
    LINALG::SerialDenseMatrix A(4, 4, true);
    A(0, 0) = 0.72903241936703114203;
    A(1, 0) = 0.87570663114228933311;
    A(2, 0) = 0.69378923027976465858;
    A(3, 0) = 0.019637190415090362652;
    A(0, 1) = 0.32707405507901372465;
    A(1, 1) = 0.64019842179333806573;
    A(2, 1) = 0.15928293569477706215;
    A(3, 1) = 0.13119201434024140151;
    A(0, 2) = 0.81862230026150939335;
    A(1, 2) = 0.36847164343389089096;
    A(2, 2) = 0.12278929762839221138;
    A(3, 2) = 0.12028240083390837511;
    A(0, 3) = 0.0052737129228371719370;
    A(1, 3) = 0.76895132151127114661;
    A(2, 3) = 0.024003735765356129168;
    A(3, 3) = 0.27465069811053651449;

    TS_ASSERT_DELTA(LINALG::DeterminantLU(A), -0.01620776397174742, 1e-14);
  }

  // fixme: fix bug in method LINALG::DeterminantSVD and remove leading "X" in method name
  void Xtest_DeterminantSVD_2x2()
  {
    LINALG::SerialDenseMatrix A(2, 2, true);
    A(0, 0) = 0.72903241936703114203;
    A(1, 0) = 0.81862230026150939335;
    A(0, 1) = 0.32707405507901372465;
    A(1, 1) = 0.0052737129228371719370;

    TS_ASSERT_DELTA(LINALG::DeterminantSVD(A), -0.2639054076334587, 1e-12);
  }

  // fixme: fix bug in method LINALG::DeterminantSVD and remove leading "X" in method name
  void Xtest_DeterminantSVD_3x3()
  {
    LINALG::SerialDenseMatrix A(3, 3, true);
    A(0, 0) = 0.72903241936703114203;
    A(1, 0) = 0.0052737129228371719370;
    A(2, 0) = 0.36847164343389089096;
    A(0, 1) = 0.32707405507901372465;
    A(1, 1) = 0.87570663114228933311;
    A(2, 1) = 0.76895132151127114661;
    A(0, 2) = 0.81862230026150939335;
    A(1, 2) = 0.64019842179333806573;
    A(2, 2) = 0.69378923027976465858;

    TS_ASSERT_DELTA(LINALG::DeterminantSVD(A), -0.1008304741716571, 1e-12);
  }

  // fixme: fix bug in method LINALG::DeterminantSVD and remove leading "X" in method name
  void Xtest_DeterminantSVD_4x4()
  {
    LINALG::SerialDenseMatrix A(4, 4, true);
    A(0, 0) = 0.72903241936703114203;
    A(1, 0) = 0.87570663114228933311;
    A(2, 0) = 0.69378923027976465858;
    A(3, 0) = 0.019637190415090362652;
    A(0, 1) = 0.32707405507901372465;
    A(1, 1) = 0.64019842179333806573;
    A(2, 1) = 0.15928293569477706215;
    A(3, 1) = 0.13119201434024140151;
    A(0, 2) = 0.81862230026150939335;
    A(1, 2) = 0.36847164343389089096;
    A(2, 2) = 0.12278929762839221138;
    A(3, 2) = 0.12028240083390837511;
    A(0, 3) = 0.0052737129228371719370;
    A(1, 3) = 0.76895132151127114661;
    A(2, 3) = 0.024003735765356129168;
    A(3, 3) = 0.27465069811053651449;

    TS_ASSERT_DELTA(LINALG::DeterminantSVD(A), -0.01620776397174742, 1e-12);
  }
};

#endif
