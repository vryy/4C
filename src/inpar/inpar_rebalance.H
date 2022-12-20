/*-----------------------------------------------------------*/
/*! \file

\brief input parameter for rebalancing the discretization

\level 2

*/
/*-----------------------------------------------------------*/

#ifndef INPAR_REBALANCE_H
#define INPAR_REBALANCE_H

#include "validparameters.H"

namespace INPAR::REBALANCE
{
  enum class RebalanceType
  {
    none,                           //< no partitioning method
    hypergraph,                     //< hypergraph based partitioning
    recursive_coordinate_bisection  //< recursive coordinate bisection, geometric based partitioning
  };

  //! set the parameters for the geometric search strategy
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);
}  // namespace INPAR::REBALANCE

#endif
