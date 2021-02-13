/*-----------------------------------------------------------------------*/
/*! \file
\brief A set of utility functions for mortar interfaces

\level 1

*/
/*-----------------------------------------------------------------------*/
#ifndef MORTAR_INTERFACE_UTILS_H
#define MORTAR_INTERFACE_UTILS_H

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

  }  // namespace INTERFACEUTILS
}  // namespace MORTAR

#endif  // #ifndef MORTAR_INTERFACE_UTILS_H
