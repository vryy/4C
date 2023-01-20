/*-----------------------------------------------------------*/
/*! \file

\brief Utility functions for the partitioning/rebalancing

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef REBALANCE_UTILS_H
#define REBALANCE_UTILS_H

#include "lib_discret.H"

namespace REBALANCE::UTILS
{
  /*! \brief Prints details on the distribution/partitioning of the distribution
   */
  void PrintParallelDistribution(const DRT::Discretization& dis);

}  // namespace REBALANCE::UTILS

#endif
