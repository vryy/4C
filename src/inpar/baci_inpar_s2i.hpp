/*----------------------------------------------------------------------*/
/*! \file
\brief input quantities and globally accessible enumerations for scatra-scatra interface coupling

\level 2


*/
/*----------------------------------------------------------------------*/
#ifndef BACI_INPAR_S2I_HPP
#define BACI_INPAR_S2I_HPP

#include "baci_config.hpp"

#include "baci_inpar_parameterlist_utils.hpp"

BACI_NAMESPACE_OPEN

// forward declaration

namespace INPUT
{
  class ConditionDefinition;
}

namespace INPAR::S2I
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
    kinetics_linearperm,
    kinetics_butlervolmer,
    kinetics_butlervolmerlinearized,
    kinetics_butlervolmerpeltier,
    kinetics_butlervolmerreduced,
    kinetics_butlervolmerreducedcapacitance,
    kinetics_butlervolmerreducedlinearized,
    kinetics_butlervolmerresistance,
    kinetics_butlervolmerreducedresistance,
    kinetics_butlervolmerreducedthermoresistance,
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
  void SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist);
}  // namespace INPAR::S2I

BACI_NAMESPACE_CLOSE

#endif  // BACI_INPAR_S2I_H
