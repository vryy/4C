/*----------------------------------------------------------------------*/
/*! \file
\brief input parameters for solid-scatra-interaction

\level 2


 *------------------------------------------------------------------------------------------------*/


#ifndef FOUR_C_INPAR_SSI_HPP
#define FOUR_C_INPAR_SSI_HPP

#include "4C_config.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace Inpar
{
  namespace SSI
  {
    /// Type of coupling strategy for SSI problems
    enum class SolutionSchemeOverFields
    {
      ssi_OneWay_ScatraToSolid,
      ssi_OneWay_SolidToScatra,
      // ssi_SequStagg_ScatraToSolid,
      // ssi_SequStagg_SolidToScatra,
      ssi_IterStagg,
      ssi_IterStaggFixedRel_ScatraToSolid,
      ssi_IterStaggFixedRel_SolidToScatra,
      ssi_IterStaggAitken_ScatraToSolid,
      ssi_IterStaggAitken_SolidToScatra,
      // IterStaggAitkenIrons,
      ssi_Monolithic
    };

    /// Type of coupling strategy between the two fields of the SSI problems
    enum class FieldCoupling
    {
      volume_match,
      volume_nonmatch,
      boundary_nonmatch,
      volumeboundary_match
    };

    //! type of scalar transport time integration
    enum class ScaTraTimIntType
    {
      standard,
      cardiac_monodomain,
      elch
    };

    /// set the ssi parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific ssi conditions
    void SetValidConditions(
        std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist);

  }  // namespace SSI

}  // namespace Inpar

FOUR_C_NAMESPACE_CLOSE

#endif
