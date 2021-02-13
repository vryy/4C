/*-----------------------------------------------------------------------*/
/*! \file
\brief A set of utility functions for mortar interfaces

\level 1

*/
/*-----------------------------------------------------------------------*/
#include <algorithm>
#include <numeric>
#include <vector>

#include "mortar_interface_utils.H"

void MORTAR::INTERFACEUTILS::ComputeParallelDistributionStatistics(
    const std::vector<int>& quantityAcrossAllRanks, int& minOverAllRanks, int& maxOverAllRanks,
    double& meanOverAllRanks)
{
  const auto first = quantityAcrossAllRanks.begin();
  const auto last = quantityAcrossAllRanks.end();

  minOverAllRanks = *std::min_element(first, last);
  maxOverAllRanks = *std::max_element(first, last);
  meanOverAllRanks = static_cast<double>(std::accumulate(first, last, 0)) /
                     static_cast<double>(quantityAcrossAllRanks.size());
}
