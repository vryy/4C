/*----------------------------------------------------------------------*/
/*! \file
\brief Unit tests for mortar interface utitlies.

\level 0
*/
/*----------------------------------------------------------------------*/
#ifndef UNIT_MORTAR_INTERFACE_UTILS_H
#define UNIT_MORTAR_INTERFACE_UTILS_H

#include <cxxtest/TestSuite.h>

#include "src/drt_mortar/mortar_interface_utils.H"

// Declaration
namespace MORTAR
{
  class InterfaceUtils_TestSuite;
}

// Implementation
class MORTAR::InterfaceUtils_TestSuite : public CxxTest::TestSuite
{
 public:
  void test_ComputeParallelDistributionStatistics()
  {
    std::vector<int> quantityAcrossAllRanks = {2, 1, 4, 3};

    int minElement = 0;
    int maxElement = 0;
    double meanOfAllElements = 0;

    MORTAR::INTERFACEUTILS::ComputeParallelDistributionStatistics(
        quantityAcrossAllRanks, minElement, maxElement, meanOfAllElements);

    TS_ASSERT_EQUALS(minElement, 1);
    TS_ASSERT_EQUALS(maxElement, 4);
    TS_ASSERT_EQUALS(meanOfAllElements, 2.5);
  }
};

#endif /* UNIT_MORTAR_INTERFACE_UTILS_H */
