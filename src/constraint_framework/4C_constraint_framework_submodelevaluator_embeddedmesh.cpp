/*-----------------------------------------------------------*/
/*! \file

\brief  Manage embedded mesh methods

\level 3
 */
/*-----------------------------------------------------------*/

#include "4C_constraint_framework_submodelevaluator_embeddedmesh.hpp"

#include "4C_constraint_framework_embeddedmesh_solid_to_solid_mortar_manager.hpp"
#include "4C_constraint_framework_embeddedmesh_solid_to_solid_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_io_visualization_manager.hpp"

FOUR_C_NAMESPACE_OPEN

CONSTRAINTS::SUBMODELEVALUATOR::EmbeddedMeshConstraintManager::EmbeddedMeshConstraintManager(
    Teuchos::RCP<Core::FE::Discretization> discret_ptr, const Epetra_Vector& dispnp)
{
  // Get the parameter lists and get information from them
  auto embedded_mesh_parameter_list = Global::Problem::instance()->embedded_mesh_params();
  auto xfem_parameter_list = Global::Problem::instance()->xfem_general_params();
  auto cut_parameter_list = Global::Problem::instance()->cut_general_params();

  auto embedded_mesh_coupling_strategy =
      Teuchos::getIntegralValue<Inpar::CONSTRAINTS::EmbeddedMeshCouplingStrategy>(
          embedded_mesh_parameter_list, "COUPLING_STRATEGY");

  auto embedded_mesh_constraint_enforcement =
      Teuchos::getIntegralValue<Inpar::CONSTRAINTS::EmbeddedMeshConstraintEnforcement>(
          embedded_mesh_parameter_list, "CONSTRAINT_ENFORCEMENT");

  auto embedded_mesh_constraint_penalty_parameter =
      embedded_mesh_parameter_list.get<double>("CONSTRAINT_ENFORCEMENT_PENALTYPARAM");

  auto nodal_dofset_strategy = Core::UTILS::integral_value<Cut::NodalDofSetStrategy>(
      xfem_parameter_list, "NODAL_DOFSET_STRATEGY");
  auto volume_cell_gauss_point_by = Core::UTILS::integral_value<Cut::VCellGaussPts>(
      xfem_parameter_list, "VOLUME_GAUSS_POINTS_BY");
  auto bound_cell_gauss_point_by = Core::UTILS::integral_value<Cut::BCellGaussPts>(
      xfem_parameter_list, "BOUNDARY_GAUSS_POINTS_BY");

  bool gmsh_cut_out = (Core::UTILS::integral_value<int>(xfem_parameter_list, "GMSH_CUT_OUT"));
  bool cut_screen_output = (Core::UTILS::integral_value<int>(xfem_parameter_list, "PRINT_OUTPUT"));

  // Initialize embedded mesh coupling parameters
  CONSTRAINTS::EMBEDDEDMESH::EmbeddedMeshParams embedded_mesh_coupling_params = {
      .embedded_mesh_coupling_strategy_ = embedded_mesh_coupling_strategy,
      .embedded_mesh_constraint_enforcement_ = embedded_mesh_constraint_enforcement,
      .embedded_mesh_constraint_penalty_parameter_ = embedded_mesh_constraint_penalty_parameter,
      .embedded_mesh_mortar_shape_function_ =
          Inpar::CONSTRAINTS::SolidToSolidMortarShapefunctions::none,
      .xfem_nodal_dof_set_strategy_ = nodal_dofset_strategy,
      .xfem_volume_cell_gauss_point_by_ = volume_cell_gauss_point_by,
      .xfem_bcell_gauss_point_by_ = bound_cell_gauss_point_by,
      .gmsh_cut_out_ = gmsh_cut_out,
      .cut_screen_output_ = cut_screen_output,
      .cut_params_ = cut_parameter_list};

  // Initialize visualization manager
  auto visualization_manager = Teuchos::rcp(new Core::IO::VisualizationManager(
      Core::IO::visualization_parameters_factory(
          Global::Problem::instance()->io_params().sublist("RUNTIME VTK OUTPUT"),
          *Global::Problem::instance()->output_control_file(), 0.0),  // Fix time
      discret_ptr->get_comm(), "embedded_mesh"));

  mortar_manager_ = Teuchos::rcp<CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager>(
      new CONSTRAINTS::EMBEDDEDMESH::SolidToSolidMortarManager(discret_ptr, dispnp,
          embedded_mesh_coupling_params, visualization_manager,
          discret_ptr->dof_row_map()->MaxAllGID() + 1));
}

bool CONSTRAINTS::SUBMODELEVALUATOR::EmbeddedMeshConstraintManager::evaluate_force_stiff(
    const Epetra_Vector& displacement_vector,
    Teuchos::RCP<Solid::TimeInt::BaseDataGlobalState>& global_state_ptr,
    Teuchos::RCP<Core::LinAlg::SparseMatrix> me_stiff_ptr, Teuchos::RCP<Epetra_Vector> me_force_ptr)
{
  // Evaluate the global mortar matrices
  mortar_manager_->evaluate_global_coupling_contributions(displacement_vector);
  mortar_manager_->add_global_force_stiffness_penalty_contributions(
      global_state_ptr, me_stiff_ptr, me_force_ptr);

  return true;
}

void CONSTRAINTS::SUBMODELEVALUATOR::EmbeddedMeshConstraintManager::runtime_output_step_state(
    std::pair<double, int> output_time_and_step)
{
  // Write runtime output for the embedded mesh method
  mortar_manager_->write_output(output_time_and_step.first, output_time_and_step.second);
}

FOUR_C_NAMESPACE_CLOSE
