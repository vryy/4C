/*----------------------------------------------------------------------*/
/*! \file
\brief Unit tests for linalg routines.

\level 0
\maintainer Ivo Steinbrecher
*/


#ifndef UNIT_LINALG_UTILS_DENSEMATRIX_INVERSE_H
#define UNIT_LINALG_UTILS_DENSEMATRIX_INVERSE_H


#include <cxxtest/TestSuite.h>

#include "src/linalg/linalg_utils_densematrix_inverse.H"


// Forward declaration.
namespace LINALG
{
  class Linalg_Inverse_TestSuite;
}


/**
 * \brief Class to test linalg routines.
 *
 * \note The values for the matrix used in tests below are generated with Mathematica:
 *       > SeedRandom[666];
 *       > A = Table[RandomReal[WorkingPrecision->50], {i, n}, {j, n}];
 *       where n needs do be replace by the dimension, e.g., n=2, n=3, or n=4
 */
class LINALG::Linalg_Inverse_TestSuite : public CxxTest::TestSuite
{
 public:
  void test_Inverse2x2ReorderMatrixEntries()
  {
    LINALG::Matrix<2, 2, double> A;
    A(0, 0) = 0.72903241936703114203;
    A(1, 0) = 0.81862230026150939335;
    A(0, 1) = 0.32707405507901372465;
    A(1, 1) = 0.0052737129228371719370;

    LINALG::Matrix<2, 2, double> B(A);

    LINALG::Inverse2x2ReorderMatrixEntries(A);

    TS_ASSERT_DELTA(A(0, 0), B(1, 1), 1e-14);
    TS_ASSERT_DELTA(A(1, 1), B(0, 0), 1e-14);
    TS_ASSERT_DELTA(A(0, 1), -B(0, 1), 1e-14);
    TS_ASSERT_DELTA(A(1, 0), -B(1, 0), 1e-14);
  }

  void test_Inverse2x2()
  {
    LINALG::Matrix<2, 2, double> A;
    A(0, 0) = 0.72903241936703114203;
    A(1, 0) = 0.81862230026150939335;
    A(0, 1) = 0.32707405507901372465;
    A(1, 1) = 0.0052737129228371719370;

    LINALG::Inverse2x2(A);

    TS_ASSERT_DELTA(A(0, 0), -0.019983345434747187407, 1e-14);
    TS_ASSERT_DELTA(A(1, 0), 3.1019534900872646712, 1e-14);
    TS_ASSERT_DELTA(A(0, 1), 1.2393609438018440619, 1e-14);
    TS_ASSERT_DELTA(A(1, 1), -2.7624762444413145470, 1e-14);
  }

  void test_Inverse2x2_singular()
  {
    LINALG::Matrix<2, 2, double> A_singular(true);
    A_singular(0, 0) = 1.0;
    TS_ASSERT_THROWS_ANYTHING(LINALG::Inverse2x2(A_singular));
  }

  void test_Inverse2x2DoNotThrowErrorOnZeroDeterminant()
  {
    LINALG::Matrix<2, 2, double> A;
    A(0, 0) = 0.72903241936703114203;
    A(1, 0) = 0.81862230026150939335;
    A(0, 1) = 0.32707405507901372465;
    A(1, 1) = 0.0052737129228371719370;

    TS_ASSERT(LINALG::Inverse2x2DoNotThrowErrorOnZeroDeterminant(A, 1e-12));

    TS_ASSERT_DELTA(A(0, 0), -0.019983345434747187407, 1e-14);
    TS_ASSERT_DELTA(A(1, 0), 3.1019534900872646712, 1e-14);
    TS_ASSERT_DELTA(A(0, 1), 1.2393609438018440619, 1e-14);
    TS_ASSERT_DELTA(A(1, 1), -2.7624762444413145470, 1e-14);
  }

  void test_Inverse2x2DoNotThrowErrorOnZeroDeterminant_singular()
  {
    LINALG::Matrix<2, 2, double> A_singular(true);
    A_singular(0, 0) = 1.0;
    TS_ASSERT(not LINALG::Inverse2x2DoNotThrowErrorOnZeroDeterminant(A_singular, 1e-12));
  }

  void test_Inverse3x3ReorderMatrixEntries()
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

    LINALG::Matrix<3, 3, double> B(A);

    LINALG::Inverse3x3ReorderMatrixEntries(A);

    TS_ASSERT_DELTA(A(0, 0), B(1, 1) * B(2, 2) - B(2, 1) * B(1, 2), 1e-14);
    TS_ASSERT_DELTA(A(1, 0), -B(1, 0) * B(2, 2) + B(2, 0) * B(1, 2), 1e-14);
    TS_ASSERT_DELTA(A(2, 0), B(1, 0) * B(2, 1) - B(2, 0) * B(1, 1), 1e-14);
    TS_ASSERT_DELTA(A(0, 1), -B(0, 1) * B(2, 2) + B(2, 1) * B(0, 2), 1e-14);
    TS_ASSERT_DELTA(A(1, 1), B(0, 0) * B(2, 2) - B(2, 0) * B(0, 2), 1e-14);
    TS_ASSERT_DELTA(A(2, 1), -B(0, 0) * B(2, 1) + B(2, 0) * B(0, 1), 1e-14);
    TS_ASSERT_DELTA(A(0, 2), B(0, 1) * B(1, 2) - B(1, 1) * B(0, 2), 1e-14);
    TS_ASSERT_DELTA(A(1, 2), -B(0, 0) * B(1, 2) + B(1, 0) * B(0, 2), 1e-14);
    TS_ASSERT_DELTA(A(2, 2), B(0, 0) * B(1, 1) - B(1, 0) * B(0, 1), 1e-14);
  }

  void test_Inverse3x3()
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

    LINALG::Inverse3x3(A);

    TS_ASSERT_DELTA(A(0, 0), -1.1432496777455423383, 1e-14);
    TS_ASSERT_DELTA(A(1, 0), -2.3032334349351217883, 1e-14);
    TS_ASSERT_DELTA(A(2, 0), 3.1599358789014842560, 1e-14);
    TS_ASSERT_DELTA(A(0, 1), -3.9924461924239508629, 1e-14);
    TS_ASSERT_DELTA(A(1, 1), -2.0247424048120138037, 1e-14);
    TS_ASSERT_DELTA(A(2, 1), 4.3644833698599101987, 1e-14);
    TS_ASSERT_DELTA(A(0, 2), 5.0330089889776210806, 1e-14);
    TS_ASSERT_DELTA(A(1, 2), 4.5859967347165455115, 1e-14);
    TS_ASSERT_DELTA(A(2, 2), -6.3144960342296159564, 1e-14);
  }

  void test_Inverse3x3_singular()
  {
    LINALG::Matrix<3, 3, double> A_singular(true);
    A_singular(0, 0) = 1.0;
    A_singular(1, 1) = 1.0;
    TS_ASSERT_THROWS_ANYTHING(LINALG::Inverse3x3(A_singular));
  }

  void test_Inverse3x3DoNotThrowErrorOnZeroDeterminant()
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

    TS_ASSERT(LINALG::Inverse3x3DoNotThrowErrorOnZeroDeterminant(A, 1e-12));

    TS_ASSERT_DELTA(A(0, 0), -1.1432496777455423383, 1e-14);
    TS_ASSERT_DELTA(A(1, 0), -2.3032334349351217883, 1e-14);
    TS_ASSERT_DELTA(A(2, 0), 3.1599358789014842560, 1e-14);
    TS_ASSERT_DELTA(A(0, 1), -3.9924461924239508629, 1e-14);
    TS_ASSERT_DELTA(A(1, 1), -2.0247424048120138037, 1e-14);
    TS_ASSERT_DELTA(A(2, 1), 4.3644833698599101987, 1e-14);
    TS_ASSERT_DELTA(A(0, 2), 5.0330089889776210806, 1e-14);
    TS_ASSERT_DELTA(A(1, 2), 4.5859967347165455115, 1e-14);
    TS_ASSERT_DELTA(A(2, 2), -6.3144960342296159564, 1e-14);
  }

  void test_Inverse3x3DoNotThrowErrorOnZeroDeterminant_singular()
  {
    LINALG::Matrix<3, 3, double> A_singular(true);
    A_singular(0, 0) = 1.0;
    A_singular(1, 1) = 1.0;
    TS_ASSERT(not LINALG::Inverse3x3DoNotThrowErrorOnZeroDeterminant(A_singular, 1e-12));
  }

  void test_Inverse4x4()
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

    LINALG::Inverse4x4(A);

    TS_ASSERT_DELTA(A(0, 0), -0.43977637337572104289, 1e-14);
    TS_ASSERT_DELTA(A(1, 0), 1.1896143886050492292, 1e-14);
    TS_ASSERT_DELTA(A(2, 0), 1.1445983336121054974, 1e-14);
    TS_ASSERT_DELTA(A(3, 0), -1.0380711684475831990, 1e-14);
    TS_ASSERT_DELTA(A(0, 1), -1.4272443093098390369, 1e-14);
    TS_ASSERT_DELTA(A(1, 1), 8.1870548809629493769, 1e-14);
    TS_ASSERT_DELTA(A(2, 1), -1.9810792318962757373, 1e-14);
    TS_ASSERT_DELTA(A(3, 1), -2.9410454529305114157, 1e-14);
    TS_ASSERT_DELTA(A(0, 2), 3.6005221408100166266, 1e-14);
    TS_ASSERT_DELTA(A(1, 2), -10.961503276990213139, 1e-14);
    TS_ASSERT_DELTA(A(2, 2), 1.1442558862090945693, 1e-14);
    TS_ASSERT_DELTA(A(3, 2), 4.4774097409218651830, 1e-14);
    TS_ASSERT_DELTA(A(0, 3), 3.6896853967343974819, 1e-14);
    TS_ASSERT_DELTA(A(1, 3), -21.986484534963294880, 1e-14);
    TS_ASSERT_DELTA(A(2, 3), 5.4245294283636667402, 1e-14);
    TS_ASSERT_DELTA(A(3, 3), 11.503778211354084434, 1e-14);
  }

  void test_Inverse4x4_singular()
  {
    LINALG::Matrix<4, 4, double> A_singular(true);
    A_singular(0, 0) = 1.0;
    A_singular(1, 1) = 1.0;
    A_singular(2, 2) = 1.0;
    TS_ASSERT_THROWS_ANYTHING(LINALG::Inverse4x4(A_singular));
  }

  void test_Inverse4x4DoNotThrowErrorOnZeroDeterminant()
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

    TS_ASSERT(LINALG::Inverse4x4DoNotThrowErrorOnZeroDeterminant(A, 1e-12));

    TS_ASSERT_DELTA(A(0, 0), -0.43977637337572104289, 1e-14);
    TS_ASSERT_DELTA(A(1, 0), 1.1896143886050492292, 1e-14);
    TS_ASSERT_DELTA(A(2, 0), 1.1445983336121054974, 1e-14);
    TS_ASSERT_DELTA(A(3, 0), -1.0380711684475831990, 1e-14);
    TS_ASSERT_DELTA(A(0, 1), -1.4272443093098390369, 1e-14);
    TS_ASSERT_DELTA(A(1, 1), 8.1870548809629493769, 1e-14);
    TS_ASSERT_DELTA(A(2, 1), -1.9810792318962757373, 1e-14);
    TS_ASSERT_DELTA(A(3, 1), -2.9410454529305114157, 1e-14);
    TS_ASSERT_DELTA(A(0, 2), 3.6005221408100166266, 1e-14);
    TS_ASSERT_DELTA(A(1, 2), -10.961503276990213139, 1e-14);
    TS_ASSERT_DELTA(A(2, 2), 1.1442558862090945693, 1e-14);
    TS_ASSERT_DELTA(A(3, 2), 4.4774097409218651830, 1e-14);
    TS_ASSERT_DELTA(A(0, 3), 3.6896853967343974819, 1e-14);
    TS_ASSERT_DELTA(A(1, 3), -21.986484534963294880, 1e-14);
    TS_ASSERT_DELTA(A(2, 3), 5.4245294283636667402, 1e-14);
    TS_ASSERT_DELTA(A(3, 3), 11.503778211354084434, 1e-14);
  }

  void test_Inverse4x4DoNotThrowErrorOnZeroDeterminant_singular()
  {
    LINALG::Matrix<4, 4, double> A_singular(true);
    A_singular(0, 0) = 1.0;
    A_singular(1, 1) = 1.0;
    A_singular(2, 2) = 1.0;
    TS_ASSERT(not LINALG::Inverse4x4DoNotThrowErrorOnZeroDeterminant(A_singular, 1e-12));
  }
};

#endif
