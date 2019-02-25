/*!
\file unit_linalg.H

\brief Unit tests for linalg routines.

\level 0
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef UNIT_LINALG_H_
#define UNIT_LINALG_H_


#include <cxxtest/TestSuite.h>

#include "src/linalg/linalg_utils.H"


// Forward declaration.
namespace LINALG
{
  class Linalg_TestSuite;
}


/**
 * \brief Class to test linalg routines.
 */
class LINALG::Linalg_TestSuite : public CxxTest::TestSuite
{
 public:
  /**
   * \brief Test inversion of 2x2 matrix.
   *
   * The values for this matrix are generated with Mathematica:
   * > SeedRandom[666];
   * > A = Table[RandomReal[WorkingPrecision->50], {i, 2}, {j, 2}];
   */
  void test_inverse_2x2()
  {
    // Variable definition.
    LINALG::TMatrix<double, 2, 2> A1;
    LINALG::TMatrix<double, 2, 2> A2;

    // Set the matrix values.
    A1(0, 0) = A2(0, 0) = 0.72903241936703114203;
    A1(1, 0) = A2(1, 0) = 0.81862230026150939335;
    A1(0, 1) = A2(0, 1) = 0.32707405507901372465;
    A1(1, 1) = A2(1, 1) = 0.0052737129228371719370;

    // Check the determinant function.
    TS_ASSERT_DELTA(LINALG::Determinant2x2(A1), -0.2639054076334587, 1e-14);

    // Invert the matrices.
    LINALG::Inverse2x2(A1);
    TS_ASSERT_EQUALS(true, LINALG::Inverse2x2DoNotThrowErrorOnZeroDeterminant(A2, 1e-12));

    // Check the inverse matrix.
    TS_ASSERT_DELTA(A1(0, 0), -0.019983345434747187407, 1e-14);
    TS_ASSERT_DELTA(A1(1, 0), 3.1019534900872646712, 1e-14);
    TS_ASSERT_DELTA(A1(0, 1), 1.2393609438018440619, 1e-14);
    TS_ASSERT_DELTA(A1(1, 1), -2.7624762444413145470, 1e-14);
    TS_ASSERT_DELTA(A2(0, 0), -0.019983345434747187407, 1e-14);
    TS_ASSERT_DELTA(A2(1, 0), 3.1019534900872646712, 1e-14);
    TS_ASSERT_DELTA(A2(0, 1), 1.2393609438018440619, 1e-14);
    TS_ASSERT_DELTA(A2(1, 1), -2.7624762444413145470, 1e-14);

    // Check that the default inverse function throws an error if the matrix is singular, and that
    // the DoNotThrowErrorOnZeroDeterminant does not.
    LINALG::TMatrix<double, 2, 2> A_singular(true);
    A_singular(0, 0) = 1.;
    TS_ASSERT_THROWS_ANYTHING(LINALG::Inverse2x2(A_singular));
    TS_ASSERT_EQUALS(false, LINALG::Inverse2x2DoNotThrowErrorOnZeroDeterminant(A_singular, 1e-12));
  }

  /**
   * \brief Test inversion of 3x3 matrix.
   *
   * The values for this matrix are generated with Mathematica:
   * > SeedRandom[666];
   * > A = Table[RandomReal[WorkingPrecision->50], {i, 3}, {j, 3}];
   */
  void test_inverse_3x3()
  {
    // Variable definition.
    LINALG::TMatrix<double, 3, 3> A1;
    LINALG::TMatrix<double, 3, 3> A2;

    // Set the matrix values.
    A1(0, 0) = A2(0, 0) = 0.72903241936703114203;
    A1(1, 0) = A2(1, 0) = 0.0052737129228371719370;
    A1(2, 0) = A2(2, 0) = 0.36847164343389089096;
    A1(0, 1) = A2(0, 1) = 0.32707405507901372465;
    A1(1, 1) = A2(1, 1) = 0.87570663114228933311;
    A1(2, 1) = A2(2, 1) = 0.76895132151127114661;
    A1(0, 2) = A2(0, 2) = 0.81862230026150939335;
    A1(1, 2) = A2(1, 2) = 0.64019842179333806573;
    A1(2, 2) = A2(2, 2) = 0.69378923027976465858;

    // Check the determinant function.
    TS_ASSERT_DELTA(LINALG::Determinant3x3(A1), -0.1008304741716571, 1e-14);

    // Invert the matrices.
    LINALG::Inverse3x3(A1);
    TS_ASSERT_EQUALS(true, LINALG::Inverse3x3DoNotThrowErrorOnZeroDeterminant(A2, 1e-12));

    // Check the inverse matrix.
    TS_ASSERT_DELTA(A1(0, 0), -1.1432496777455423383, 1e-14);
    TS_ASSERT_DELTA(A1(1, 0), -2.3032334349351217883, 1e-14);
    TS_ASSERT_DELTA(A1(2, 0), 3.1599358789014842560, 1e-14);
    TS_ASSERT_DELTA(A1(0, 1), -3.9924461924239508629, 1e-14);
    TS_ASSERT_DELTA(A1(1, 1), -2.0247424048120138037, 1e-14);
    TS_ASSERT_DELTA(A1(2, 1), 4.3644833698599101987, 1e-14);
    TS_ASSERT_DELTA(A1(0, 2), 5.0330089889776210806, 1e-14);
    TS_ASSERT_DELTA(A1(1, 2), 4.5859967347165455115, 1e-14);
    TS_ASSERT_DELTA(A1(2, 2), -6.3144960342296159564, 1e-14);
    TS_ASSERT_DELTA(A2(0, 0), -1.1432496777455423383, 1e-14);
    TS_ASSERT_DELTA(A2(1, 0), -2.3032334349351217883, 1e-14);
    TS_ASSERT_DELTA(A2(2, 0), 3.1599358789014842560, 1e-14);
    TS_ASSERT_DELTA(A2(0, 1), -3.9924461924239508629, 1e-14);
    TS_ASSERT_DELTA(A2(1, 1), -2.0247424048120138037, 1e-14);
    TS_ASSERT_DELTA(A2(2, 1), 4.3644833698599101987, 1e-14);
    TS_ASSERT_DELTA(A2(0, 2), 5.0330089889776210806, 1e-14);
    TS_ASSERT_DELTA(A2(1, 2), 4.5859967347165455115, 1e-14);
    TS_ASSERT_DELTA(A2(2, 2), -6.3144960342296159564, 1e-14);

    // Check that the default inverse function throws an error if the matrix is singular, and that
    // the DoNotThrowErrorOnZeroDeterminant does not.
    LINALG::TMatrix<double, 3, 3> A_singular(true);
    A_singular(0, 0) = 1.;
    A_singular(1, 1) = 1.;
    TS_ASSERT_THROWS_ANYTHING(LINALG::Inverse3x3(A_singular));
    TS_ASSERT_EQUALS(false, LINALG::Inverse3x3DoNotThrowErrorOnZeroDeterminant(A_singular, 1e-12));
  }

  /**
   * \brief Test inversion of 4x4 matrix.
   *
   * The values for this matrix are generated with Mathematica:
   * > SeedRandom[666];
   * > A = Table[RandomReal[WorkingPrecision->50], {i, 4}, {j, 4}];
   */
  void test_inverse_4x4()
  {
    // Variable definition.
    LINALG::TMatrix<double, 4, 4> A1;
    LINALG::TMatrix<double, 4, 4> A2;

    // Set the matrix values.
    A1(0, 0) = A2(0, 0) = 0.72903241936703114203;
    A1(1, 0) = A2(1, 0) = 0.87570663114228933311;
    A1(2, 0) = A2(2, 0) = 0.69378923027976465858;
    A1(3, 0) = A2(3, 0) = 0.019637190415090362652;
    A1(0, 1) = A2(0, 1) = 0.32707405507901372465;
    A1(1, 1) = A2(1, 1) = 0.64019842179333806573;
    A1(2, 1) = A2(2, 1) = 0.15928293569477706215;
    A1(3, 1) = A2(3, 1) = 0.13119201434024140151;
    A1(0, 2) = A2(0, 2) = 0.81862230026150939335;
    A1(1, 2) = A2(1, 2) = 0.36847164343389089096;
    A1(2, 2) = A2(2, 2) = 0.12278929762839221138;
    A1(3, 2) = A2(3, 2) = 0.12028240083390837511;
    A1(0, 3) = A2(0, 3) = 0.0052737129228371719370;
    A1(1, 3) = A2(1, 3) = 0.76895132151127114661;
    A1(2, 3) = A2(2, 3) = 0.024003735765356129168;
    A1(3, 3) = A2(3, 3) = 0.27465069811053651449;

    // Check the determinant function.
    TS_ASSERT_DELTA(LINALG::Determinant4x4(A1), -0.01620776397174742, 1e-14);

    // Invert the matrices.
    LINALG::Inverse4x4(A1);
    TS_ASSERT_EQUALS(true, LINALG::Inverse4x4DoNotThrowErrorOnZeroDeterminant(A2, 1e-12));

    // Check the inverse matrix.
    TS_ASSERT_DELTA(A1(0, 0), -0.43977637337572104289, 1e-14);
    TS_ASSERT_DELTA(A1(1, 0), 1.1896143886050492292, 1e-14);
    TS_ASSERT_DELTA(A1(2, 0), 1.1445983336121054974, 1e-14);
    TS_ASSERT_DELTA(A1(3, 0), -1.0380711684475831990, 1e-14);
    TS_ASSERT_DELTA(A1(0, 1), -1.4272443093098390369, 1e-14);
    TS_ASSERT_DELTA(A1(1, 1), 8.1870548809629493769, 1e-14);
    TS_ASSERT_DELTA(A1(2, 1), -1.9810792318962757373, 1e-14);
    TS_ASSERT_DELTA(A1(3, 1), -2.9410454529305114157, 1e-14);
    TS_ASSERT_DELTA(A1(0, 2), 3.6005221408100166266, 1e-14);
    TS_ASSERT_DELTA(A1(1, 2), -10.961503276990213139, 1e-14);
    TS_ASSERT_DELTA(A1(2, 2), 1.1442558862090945693, 1e-14);
    TS_ASSERT_DELTA(A1(3, 2), 4.4774097409218651830, 1e-14);
    TS_ASSERT_DELTA(A1(0, 3), 3.6896853967343974819, 1e-14);
    TS_ASSERT_DELTA(A1(1, 3), -21.986484534963294880, 1e-14);
    TS_ASSERT_DELTA(A1(2, 3), 5.4245294283636667402, 1e-14);
    TS_ASSERT_DELTA(A1(3, 3), 11.503778211354084434, 1e-14);
    TS_ASSERT_DELTA(A2(0, 0), -0.43977637337572104289, 1e-14);
    TS_ASSERT_DELTA(A2(1, 0), 1.1896143886050492292, 1e-14);
    TS_ASSERT_DELTA(A2(2, 0), 1.1445983336121054974, 1e-14);
    TS_ASSERT_DELTA(A2(3, 0), -1.0380711684475831990, 1e-14);
    TS_ASSERT_DELTA(A2(0, 1), -1.4272443093098390369, 1e-14);
    TS_ASSERT_DELTA(A2(1, 1), 8.1870548809629493769, 1e-14);
    TS_ASSERT_DELTA(A2(2, 1), -1.9810792318962757373, 1e-14);
    TS_ASSERT_DELTA(A2(3, 1), -2.9410454529305114157, 1e-14);
    TS_ASSERT_DELTA(A2(0, 2), 3.6005221408100166266, 1e-14);
    TS_ASSERT_DELTA(A2(1, 2), -10.961503276990213139, 1e-14);
    TS_ASSERT_DELTA(A2(2, 2), 1.1442558862090945693, 1e-14);
    TS_ASSERT_DELTA(A2(3, 2), 4.4774097409218651830, 1e-14);
    TS_ASSERT_DELTA(A2(0, 3), 3.6896853967343974819, 1e-14);
    TS_ASSERT_DELTA(A2(1, 3), -21.986484534963294880, 1e-14);
    TS_ASSERT_DELTA(A2(2, 3), 5.4245294283636667402, 1e-14);
    TS_ASSERT_DELTA(A2(3, 3), 11.503778211354084434, 1e-14);

    // Check that the default inverse function throws an error if the matrix is singular, and that
    // the DoNotThrowErrorOnZeroDeterminant does not.
    LINALG::TMatrix<double, 4, 4> A_singular(true);
    A_singular(0, 0) = 1.;
    A_singular(1, 1) = 1.;
    A_singular(2, 2) = 1.;
    TS_ASSERT_THROWS_ANYTHING(LINALG::Inverse4x4(A_singular));
    TS_ASSERT_EQUALS(false, LINALG::Inverse4x4DoNotThrowErrorOnZeroDeterminant(A_singular, 1e-12));
  }
};

#endif
