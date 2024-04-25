/*-----------------------------------------------------------*/
/*! \file

\brief Input parameters for the constraint framework

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_INPAR_CONSTRAINT_FRAMEWORK_HPP
#define FOUR_C_INPAR_CONSTRAINT_FRAMEWORK_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace INPAR::CONSTRAINTS
{
  /// type of the submodel for constraintmodels
  enum class SubModelType
  {
    submodel_undefined,  ///< default
    submodel_pbc_rve     ///< apply periodic displacement bcs on
  };

  /// type of employed constraint enforcement strategy
  enum class Strategy
  {
    penalty_regularization  ///< penalty method
  };
}  // namespace INPAR::CONSTRAINTS

FOUR_C_NAMESPACE_CLOSE

#endif
