/*----------------------------------------------------------------------*/
/*! \file
\brief input quantities and globally accessible enumerations for scatra-scatra interface coupling

\level 2

\maintainer Christoph Schmidt

*/
/*----------------------------------------------------------------------*/
#ifndef INPAR_S2I_H
#define INPAR_S2I_H

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

namespace INPAR
{
  namespace S2I
  {
    //! type of interface side
    enum InterfaceSides
    {
      side_undefined,
      side_slave,
      side_master
    };

    //! type of mesh coupling
    enum CouplingType
    {
      coupling_undefined,
      coupling_matching_nodes,
      coupling_mortar_standard,
      coupling_mortar_saddlepoint_petrov,
      coupling_mortar_saddlepoint_bubnov,
      coupling_mortar_condensed_petrov,
      coupling_mortar_condensed_bubnov,
      coupling_nts_standard
    };

    //! type of interface layer growth evaluation
    enum GrowthEvaluation
    {
      growth_evaluation_none,
      growth_evaluation_monolithic,
      growth_evaluation_semi_implicit
    };

    //! models for interface layer growth kinetics
    enum GrowthKineticModels
    {
      growth_kinetics_butlervolmer
    };

    //! models for interface kinetics
    enum KineticModels
    {
      kinetics_constperm,
      kinetics_butlervolmer,
      kinetics_butlervolmerpeltier,
      kinetics_butlervolmerreduced,
      kinetics_butlervolmerresistance,
      kinetics_butlervolmerreducedwithresistance,
      kinetics_constantinterfaceresistance,
      kinetics_nointerfaceflux
    };

    //! actions for mortar cell evaluation
    enum EvaluationActions
    {
      evaluate_condition,
      evaluate_condition_nts,
      evaluate_condition_od,
      evaluate_mortar_matrices,
      evaluate_nodal_area_fractions
    };

    //! regularization types for plating reaction
    enum RegularizationType
    {
      regularization_undefined,
      regularization_none,
      regularization_hein,
      regularization_polynomial,
      regularization_trigonometrical
    };

    //! set valid parameters for scatra-scatra interface coupling
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    //! set valid conditions for scatra-scatra interface coupling
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);
  }  // namespace S2I
}  // namespace INPAR
#endif  // #ifndef INPAR_S2I_H
