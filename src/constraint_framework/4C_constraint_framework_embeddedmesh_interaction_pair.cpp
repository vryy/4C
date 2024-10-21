#include "4C_constraint_framework_embeddedmesh_interaction_pair.hpp"

#include "4C_constraint_framework_embeddedmesh_solid_to_solid_utils.hpp"
#include "4C_cut_boundarycell.hpp"
#include "4C_cut_cutwizard.hpp"

#include <Teuchos_ENull.hpp>

FOUR_C_NAMESPACE_OPEN

CONSTRAINTS::EMBEDDEDMESH::SolidInteractionPair::SolidInteractionPair(
    Teuchos::RCP<Core::Elements::Element> element1, Core::Elements::Element* element2,
    CONSTRAINTS::EMBEDDEDMESH::EmbeddedMeshParams& params_ptr,
    Teuchos::RCP<Cut::CutWizard> cutwizard_ptr,
    std::vector<Teuchos::RCP<Cut::BoundaryCell>>& boundary_cells)
    : params_(params_ptr),
      element1_(element1),
      element2_(element2),
      cutwizard_ptr_(cutwizard_ptr),
      boundary_cells_(boundary_cells)
{
  // empty constructor
}

FOUR_C_NAMESPACE_CLOSE
