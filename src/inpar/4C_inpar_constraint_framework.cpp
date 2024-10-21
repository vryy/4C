#include "4C_inpar_constraint_framework.hpp"

#include "4C_utils_parameter_list.hpp"


FOUR_C_NAMESPACE_OPEN

/**
 *
 */
void Inpar::CONSTRAINTS::set_valid_parameters(Teuchos::ParameterList& list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& embeddedmeshcoupling = list.sublist("EMBEDDED MESH COUPLING", false, "");
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

    Core::Utils::double_parameter("CONSTRAINT_ENFORCEMENT_PENALTYPARAM", 0.0,
        "Penalty parameter for the constraint enforcement in embedded mesh coupling",
        &embeddedmeshcoupling);
  }
}

FOUR_C_NAMESPACE_CLOSE