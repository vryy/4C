/*----------------------------------------------------------------------*/
/*! \file
\brief Unit tests for mortar interface utitlies.

\level 1
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "src/drt_mortar/mortar_interface_utils.H"

namespace
{
  TEST(InterfaceUtilsTestSuite, TestComputeParallelDistributionStatistics)
  {
    std::vector<int> quantityAcrossAllRanks = {2, 1, 4, 3};

    int minElement = 0;
    int maxElement = 0;
    double meanOfAllElements = 0;

    MORTAR::INTERFACEUTILS::ComputeParallelDistributionStatistics(
        quantityAcrossAllRanks, minElement, maxElement, meanOfAllElements);

    EXPECT_EQ(minElement, 1);
    EXPECT_EQ(maxElement, 4);
    EXPECT_EQ(meanOfAllElements, 2.5);
  }
}  // namespace
