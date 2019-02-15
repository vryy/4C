/*!---------------------------------------------------------------------*
\file unit_so_hex8_determinant_analysis.H

\brief Unittests for the determinant validity check

\level 3

\maintainer Anh-Tu Vuong
*-----------------------------------------------------------------------*/

#ifndef UNITTESTS_SRC_DRT_SO3_UNIT_SO_HEX8_DETERMINANT_ANALYSIS_H_
#define UNITTESTS_SRC_DRT_SO3_UNIT_SO_HEX8_DETERMINANT_ANALYSIS_H_

#include <cxxtest/TestSuite.h>
#include "src/drt_so3/so_hex8_determinant_analysis.H"

namespace DRT
{
  namespace ELEMENTS
  {
    class So_hex8_DetAnalysis_TestSuite;
  }
}  // namespace DRT

/**
 * Class to test determinant check of hex8 element
 */
class DRT::ELEMENTS::So_hex8_DetAnalysis_TestSuite : public CxxTest::TestSuite
{
  Teuchos::RCP<DRT::ELEMENTS::So_Hex8_Determinant_Analysis> analyser_;

 public:
  /**
   * Set up the testing environment.
   */
  void setUp() override { analyser_ = DRT::ELEMENTS::So_Hex8_Determinant_Analysis::create(); }

  /**
   * Delete pointers.
   */
  void tearDown() override { analyser_ = Teuchos::null; }

  /** Test coordinates of an undeformed element (sanity check)
   *  Expected result = VALID */
  void test_element_undeformed()
  {
#ifdef DEBUG_SO_HEX8_DET_ANALYSIS
    std::cout << "\n\n" << std::string(80, '*') << "\n" << __PRETTY_FUNCTION__;
    std::cout << "\n" << std::string(80, '*') << "\n\n";
#endif

    LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> x_test;
    x_test(0, 0) = 0.0;
    x_test(0, 1) = 1.0;
    x_test(0, 2) = 1.0;
    x_test(0, 3) = 0.0;
    x_test(1, 0) = 0.0;
    x_test(1, 1) = 0.0;
    x_test(1, 2) = 1.0;
    x_test(1, 3) = 1.0;
    x_test(2, 0) = 0.0;
    x_test(2, 1) = 0.0;
    x_test(2, 2) = 0.0;
    x_test(2, 3) = 0.0;

    x_test(0, 4) = 0.0;
    x_test(0, 5) = 1.0;
    x_test(0, 6) = 1.0;
    x_test(0, 7) = 0.0;
    x_test(1, 4) = 0.0;
    x_test(1, 5) = 0.0;
    x_test(1, 6) = 1.0;
    x_test(1, 7) = 1.0;
    x_test(2, 4) = 1.0;
    x_test(2, 5) = 1.0;
    x_test(2, 6) = 1.0;
    x_test(2, 7) = 1.0;

    unsigned rc = 0;
    TS_ASSERT_EQUALS(analyser_->isValid(x_test, &rc), true);
    TS_ASSERT_EQUALS(rc, 0);
  }

  /** Test coordinates following Fig. 5 in Johnen2017
   *  Expected result = INVALID */
  void test_element_fig5_johnen2017()
  {
#ifdef DEBUG_SO_HEX8_DET_ANALYSIS
    std::cout << "\n\n" << std::string(80, '*') << "\n" << __PRETTY_FUNCTION__;
    std::cout << "\n" << std::string(80, '*') << "\n\n";
#endif

    LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> x_test;
    x_test(0, 0) = 0.0;
    x_test(0, 1) = 1.0;
    x_test(0, 2) = 1.7615170641459;
    x_test(0, 3) = 0.438888411629833;
    x_test(1, 0) = 0.0;
    x_test(1, 1) = 0.0;
    x_test(1, 2) = 0.594764272968121;
    x_test(1, 3) = 1.53098020041072;
    x_test(2, 0) = 0.0;
    x_test(2, 1) = 0.0;
    x_test(2, 2) = 0.15552188663289;
    x_test(2, 3) = 0.185631029838277;

    x_test(0, 4) = 1.3859049651391;
    x_test(0, 5) = 1.22129676447071;
    x_test(0, 6) = 1.77365642274365;
    x_test(0, 7) = 0.0769922201302364;
    x_test(1, 4) = 0.0755794018509022;
    x_test(1, 5) = 0.271876165350328;
    x_test(1, 6) = 1.25103990471942;
    x_test(1, 7) = 0.940424880836765;
    x_test(2, 4) = 1.77483024073906;
    x_test(2, 5) = 0.630922158503566;
    x_test(2, 6) = 1.83300604452892;
    x_test(2, 7) = 1.45521546591891;

    unsigned rc = 0;
    TS_ASSERT_EQUALS(analyser_->isValid(x_test, &rc), false);
    TS_ASSERT_EQUALS(rc, 3);
  }

  /** Test coordinates following Fig. 6 in Johnen2017
   *  Expected result = VALID */
  void test_element_fig6_johnen2017()
  {
#ifdef DEBUG_SO_HEX8_DET_ANALYSIS
    std::cout << "\n\n" << std::string(80, '*') << "\n" << __PRETTY_FUNCTION__;
    std::cout << "\n" << std::string(80, '*') << "\n\n";
#endif

    LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> x_test;
    x_test(0, 0) = 0.0;
    x_test(0, 1) = 1.0;
    x_test(0, 2) = 1.539;
    x_test(0, 3) = 0.166589;
    x_test(1, 0) = 0.0;
    x_test(1, 1) = 0.0;
    x_test(1, 2) = 0.704696;
    x_test(1, 3) = 1.08208;
    x_test(2, 0) = 0.0;
    x_test(2, 1) = 0.0;
    x_test(2, 2) = 1.84011;
    x_test(2, 3) = 0.162539;

    x_test(0, 4) = 0.0501127;
    x_test(0, 5) = 0.422336;
    x_test(0, 6) = 0.509917;
    x_test(0, 7) = 0.40783;
    x_test(1, 4) = 1.96347;
    x_test(1, 5) = 0.00419138;
    x_test(1, 6) = 0.0214216;
    x_test(1, 7) = 1.73452;
    x_test(2, 4) = 1.56559;
    x_test(2, 5) = 1.43038;
    x_test(2, 6) = 1.55322;
    x_test(2, 7) = 1.93234;

    unsigned rc = 0;
    TS_ASSERT_EQUALS(analyser_->isValid(x_test, &rc), true);
    TS_ASSERT_EQUALS(rc, 1);
  }

  /** Test coordinates following Fig. 7 in Johnen2017
   *  Expected result = INVALID */
  void test_element_fig7_johnen2017()
  {
#ifdef DEBUG_SO_HEX8_DET_ANALYSIS
    std::cout << "\n\n" << std::string(80, '*') << "\n" << __PRETTY_FUNCTION__;
    std::cout << "\n" << std::string(80, '*') << "\n\n";
#endif

    LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> x_test;
    x_test(0, 0) = 0.464949491866817;
    x_test(0, 1) = 0.481795709097567;
    x_test(0, 2) = 0.482406079287087;
    x_test(0, 3) = 0.466719565416425;
    x_test(1, 0) = 0.358989226966155;
    x_test(1, 1) = 0.358745078890347;
    x_test(1, 2) = 0.351664784691916;
    x_test(1, 3) = 0.339945677053133;
    x_test(2, 0) = 0.0133365886410108;
    x_test(2, 1) = 0.0163884395886105;
    x_test(2, 2) = 0.0235297708059938;
    x_test(2, 3) = 0.0278023621326335;

    x_test(0, 4) = 0.465498825037385;
    x_test(0, 5) = 0.465987121189001;
    x_test(0, 6) = 0.501998962370677;
    x_test(0, 7) = 0.487166966765343;
    x_test(1, 4) = 0.320291756950591;
    x_test(1, 5) = 0.321085238196966;
    x_test(1, 6) = 0.322367015594958;
    x_test(1, 7) = 0.308816797387616;
    x_test(2, 4) = -0.00277718436231578;
    x_test(2, 5) = -0.0042420728171636;
    x_test(2, 6) = -0.0116275521103549;
    x_test(2, 7) = 0.0115054780724508;

    unsigned rc = 0;
    TS_ASSERT_EQUALS(analyser_->isValid(x_test, &rc), false);
    TS_ASSERT_EQUALS(rc, 0);
  }
};


#endif /* UNITTESTS_SRC_DRT_SO3_UNIT_SO_HEX8_DETERMINANT_ANALYSIS_H_ */
