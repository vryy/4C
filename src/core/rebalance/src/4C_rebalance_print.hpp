#ifndef FOUR_C_REBALANCE_PRINT_HPP
#define FOUR_C_REBALANCE_PRINT_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Rebalance::Utils
{
  /*! \brief Prints details on the distribution/partitioning of the distribution
   */
  void print_parallel_distribution(const Core::FE::Discretization& dis);

}  // namespace Core::Rebalance::Utils

FOUR_C_NAMESPACE_CLOSE

#endif
