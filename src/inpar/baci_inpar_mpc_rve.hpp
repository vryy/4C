/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for multi point constraints used for periodic boundary conditions
 for representative volume elements (RVEs)
\level 3
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_INPAR_MPC_RVE_HPP
#define FOUR_C_INPAR_MPC_RVE_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_integration.hpp"
#include "baci_lib_conditiondefinition.hpp"

BACI_NAMESPACE_OPEN

namespace INPAR::RVE_MPC
{
  /// Definition Type of the MultiPoint Constraint
  enum multiPointConstraintType
  {
    none,
    coupled_equation,  ///< Manually enter the dofs as coupled equation
    periodic_rve       ///< Automatically apply the MPCs that describe periodicity

  };

  /// Definition Type of Coupled Equation
  /// (this enum represents the input file parameter COUPLED_DOF_EQUATIONS)
  enum ceType
  {
    ce_none,
    ce_linear
  };

  /// Strategy used to enforce the MPCs
  /// (this enum represents the input file parameter ENFORCEMENT)
  enum enforcementStrategy
  {
    penalty,            ///< Enforce the Multi-Point Constraint with the penalty method
    lagrangeMultiplier  ///< Enforce the Multi-Point Constraint with the Lagrange Multiplier
                        ///< Method
  };

  /// Methods used to determine the reference points for an periodic RVE
  /// (this enum represents the input file parameter RVE_REFERENCE_POINTS)
  enum rveReferenceDeformationDefinition
  {
    automatic,  ///< Automatically use the corner nodes for reference
    manual      ///< provide the reference nodes as condition

  };

  enum rveDimension
  {
    rve2d,
    rve3d
  };

  enum rveEdgeIdentifiers
  {
    Gamma_xm,
    Gamma_ym,

  };
  /// set the multi point constraint parameters
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  /// set multi point constraint specific conditions
  void SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist);
}  // namespace INPAR::RVE_MPC

BACI_NAMESPACE_CLOSE
#endif
