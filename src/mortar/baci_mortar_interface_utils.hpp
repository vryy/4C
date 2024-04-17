/*-----------------------------------------------------------------------*/
/*! \file
\brief A set of utility functions for mortar interfaces

\level 1

*/
/*-----------------------------------------------------------------------*/
#ifndef FOUR_C_MORTAR_INTERFACE_UTILS_HPP
#define FOUR_C_MORTAR_INTERFACE_UTILS_HPP

#include "baci_config.hpp"

#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace MORTAR
{
  namespace INTERFACEUTILS
  {
    /*!
    \brief Compute statistics on parallel distribution

    \param[in] quantityAcrossAllRanks Vector with number of nodes/elements per rank
    \param[out] minOverAllRanks Smallest element in input vector
    \param[out] maxOverAllRanks Largest element in input vector
    \param[out] meanOverAllRanks Mean of all elements in input vector
    */
    void ComputeParallelDistributionStatistics(const std::vector<int>& quantityAcrossAllRanks,
        int& minOverAllRanks, int& maxOverAllRanks, double& meanOverAllRanks);

    /*!
    \brief Compute and print one row of statisitcs on parallel distribution table

    \param[in] nameOfQuantity Name describing the quantity to be printed
    \param[in] quantityAcrossAllRanks Vector with number of nodes/elements per rank
    \param[in] printOnThisRank Flag to en-/disable screen output on this rank
    */
    void ComputeAndPrintRowOfParallelDistributionStatisctics(const std::string& nameOfQuantity,
        const std::vector<int>& quantityAcrossAllRanks, const bool printOnThisRank);

  }  // namespace INTERFACEUTILS
}  // namespace MORTAR

FOUR_C_NAMESPACE_CLOSE

#endif
