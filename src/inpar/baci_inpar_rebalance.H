/*-----------------------------------------------------------*/
/*! \file

\brief input parameter for rebalancing the discretization

\level 2

*/
/*-----------------------------------------------------------*/

#ifndef BACI_INPAR_REBALANCE_H
#define BACI_INPAR_REBALANCE_H

#include "baci_inpar_validparameters.H"

BACI_NAMESPACE_OPEN

namespace INPAR::REBALANCE
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

  //! set the parameters for the geometric search strategy
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);
}  // namespace INPAR::REBALANCE

BACI_NAMESPACE_CLOSE

#endif
