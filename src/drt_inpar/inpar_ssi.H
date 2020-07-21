/*----------------------------------------------------------------------*/
/*! \file
\brief input parameters for solid-scatra-interaction

\level 2


 *------------------------------------------------------------------------------------------------*/


#ifndef INPAR_SSI_H_
#define INPAR_SSI_H_

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace INPAR
{
  namespace SSI
  {
    /// Type of coupling strategy for SSI problems
    enum SolutionSchemeOverFields
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
    enum FieldCoupling
    {
      coupling_volume_match,
      coupling_volume_nonmatch,
      coupling_boundary_nonmatch,
      coupling_volumeboundary_match
    };

    //! type of scalar transport time integration
    enum ScaTraTimIntType
    {
      scatratiminttype_standard,
      scatratiminttype_cardiac_monodomain,
      scatratiminttype_elch
    };

    /// set the ssi parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific ssi conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

  }  // namespace SSI

}  // namespace INPAR

#endif /* INPAR_SSI_H_ */
