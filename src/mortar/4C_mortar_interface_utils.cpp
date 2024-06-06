/*-----------------------------------------------------------------------*/
/*! \file
\brief A set of utility functions for mortar interfaces

\level 1

*/
/*-----------------------------------------------------------------------*/
#include "4C_mortar_interface_utils.hpp"

#include <algorithm>
#include <numeric>
#include <sstream>
#include <vector>

FOUR_C_NAMESPACE_OPEN

void Mortar::InterfaceUtils::ComputeParallelDistributionStatistics(
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

void Mortar::InterfaceUtils::ComputeAndPrintRowOfParallelDistributionStatisctics(
    const std::string& nameOfQuantity, const std::vector<int>& quantityAcrossAllRanks,
    const bool printOnThisRank)
{
  if (printOnThisRank)
  {
    int minOverAllRanks = -1;
    int maxOverAllRanks = -1;
    double meanOverAllRanks = -1.0;

    Mortar::InterfaceUtils::ComputeParallelDistributionStatistics(
        quantityAcrossAllRanks, minOverAllRanks, maxOverAllRanks, meanOverAllRanks);

    std::stringstream maxToMinRatio;
    if (minOverAllRanks > 0)
      maxToMinRatio << static_cast<double>(maxOverAllRanks) / static_cast<double>(minOverAllRanks);
    else
      maxToMinRatio << "---";

    printf("    | %20s | %14d | %14d | %15.1f | %16s |\n", nameOfQuantity.c_str(), minOverAllRanks,
        maxOverAllRanks, meanOverAllRanks, maxToMinRatio.str().c_str());
  }
}

FOUR_C_NAMESPACE_CLOSE
