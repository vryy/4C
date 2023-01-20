/*-----------------------------------------------------------*/
/*! \file

\brief Data container class for rebalance parameters

\level 3

*/
/*-----------------------------------------------------------*/

#include <Teuchos_ParameterList.hpp>

#include "lib_globalproblem.H"

#include "rebalance_params.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
REBALANCE::RebalanceParams::RebalanceParams()
{
  Teuchos::ParameterList const& params_list = DRT::Problem::Instance()->MeshPartitioningParams();

  rebalance_method_ =
      Teuchos::getIntegralValue<INPAR::REBALANCE::RebalanceType>(params_list, "METHOD");

  imbalance_tolerance_ = params_list.get<double>("IMBALANCE_TOL");
}