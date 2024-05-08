/*! \file
\brief General rebalance functionality
\level 3
*/

#ifndef FOUR_C_REBALANCE_HPP
#define FOUR_C_REBALANCE_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CORE::REBALANCE
{

  enum class RebalanceType
  {
    none,                            //< no partitioning method
    hypergraph,                      //< hypergraph based partitioning
    recursive_coordinate_bisection,  //< recursive coordinate bisection, geometric based
                                     // partitioning
    monolithic  //< hypergraph based partitioning by using a global monolithic graph constructed
                // via a global collision search
  };
}

FOUR_C_NAMESPACE_CLOSE

#endif
