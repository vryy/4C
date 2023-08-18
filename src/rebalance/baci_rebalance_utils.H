/*-----------------------------------------------------------*/
/*! \file

\brief Utility functions for the partitioning/rebalancing

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_REBALANCE_UTILS_H
#define BACI_REBALANCE_UTILS_H

#include "baci_lib_discret.H"

namespace REBALANCE::UTILS
{
  /*! \brief Prints details on the distribution/partitioning of the distribution
   */
  void PrintParallelDistribution(const DRT::Discretization& dis);

}  // namespace REBALANCE::UTILS

#endif
