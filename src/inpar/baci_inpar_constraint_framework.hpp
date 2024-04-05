/*-----------------------------------------------------------*/
/*! \file

\brief Input parameters for the constraint framework

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_INPAR_CONSTRAINT_FRAMEWORK_HPP
#define FOUR_C_INPAR_CONSTRAINT_FRAMEWORK_HPP

#include "baci_config.hpp"

BACI_NAMESPACE_OPEN

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

BACI_NAMESPACE_CLOSE

#endif  // BACI_BACI_INPAR_CONSTRAINT_MODELS_H
