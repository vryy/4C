/*----------------------------------------------------------------------*/
/*! \file
\brief Unit tests for mortar interface utitlies.

\level 1
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "baci_mortar_interface_utils.hpp"

namespace
{
  using namespace FourC;

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
