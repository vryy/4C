/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for constraint framework library

\level 2


*/
/*----------------------------------------------------------------------*/

#include "4C_inpar_constraint_framework.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
void Inpar::CONSTRAINTS::set_valid_parameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& embeddedmeshcoupling = list->sublist("EMBEDDED MESH COUPLING", false, "");
  {
    setStringToIntegralParameter<EmbeddedMeshCouplingStrategy>("COUPLING_STRATEGY", "none",
        "Strategy to couple background and overlapping mesh", tuple<std::string>("none", "mortar"),
        tuple<EmbeddedMeshCouplingStrategy>(
            EmbeddedMeshCouplingStrategy::none, EmbeddedMeshCouplingStrategy::mortar),
        &embeddedmeshcoupling);

    setStringToIntegralParameter<EmbeddedMeshConstraintEnforcement>("CONSTRAINT_ENFORCEMENT",
        "none", "Apply a constraint enforcement in the embedded mesh coupling strategy",
        tuple<std::string>("none", "penalty"),
        tuple<EmbeddedMeshConstraintEnforcement>(
            EmbeddedMeshConstraintEnforcement::none, EmbeddedMeshConstraintEnforcement::penalty),
        &embeddedmeshcoupling);

    Core::UTILS::double_parameter("CONSTRAINT_ENFORCEMENT_PENALTYPARAM", 0.0,
        "Penalty parameter for the constraint enforcement in embedded mesh coupling",
        &embeddedmeshcoupling);
  }
}

FOUR_C_NAMESPACE_CLOSE