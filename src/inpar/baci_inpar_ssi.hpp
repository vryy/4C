/*----------------------------------------------------------------------*/
/*! \file
\brief input parameters for solid-scatra-interaction

\level 2


 *------------------------------------------------------------------------------------------------*/


#ifndef BACI_INPAR_SSI_HPP
#define BACI_INPAR_SSI_HPP

#include "baci_config.hpp"

#include "baci_inpar_parameterlist_utils.hpp"

BACI_NAMESPACE_OPEN

// forward declaration
namespace INPUT
{
  class ConditionDefinition;
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace INPAR
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
    void SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist);

  }  // namespace SSI

}  // namespace INPAR

BACI_NAMESPACE_CLOSE

#endif  // INPAR_SSI_H
