/*-----------------------------------------------------------*/
/*! \file

\brief Data container class for rebalance parameters

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef REBALANCE_PARAMS_H
#define REBALANCE_PARAMS_H

#include "inpar_rebalance.H"

namespace REBALANCE
{
  //! Data container for rebalance parameters
  class RebalanceParams
  {
   public:
    RebalanceParams();

   private:
    INPAR::REBALANCE::RebalanceType rebalance_method_;
    double imbalance_tolerance_;
  };
}  // namespace REBALANCE

#endif
