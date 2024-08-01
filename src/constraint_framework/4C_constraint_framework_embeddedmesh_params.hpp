/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief data container holding all solid to solid embedded mesh interaction
       input parameters

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_CONSTRAINT_FRAMEWORK_EMBEDDEDMESH_PARAMS_HPP
#define FOUR_C_CONSTRAINT_FRAMEWORK_EMBEDDEDMESH_PARAMS_HPP

#include "4C_config.hpp"

#include "4C_cut_enum.hpp"
#include "4C_fem_general_cell_type.hpp"
#include "4C_inpar_constraint_framework.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONSTRAINTS::EMBEDDEDMESH
{
  struct EmbeddedMeshParams
  {
    //! Strategy for coupling the embedded meshes
    enum Inpar::CONSTRAINTS::EmbeddedMeshCouplingStrategy embedded_mesh_coupling_strategy_ =
        Inpar::CONSTRAINTS::EmbeddedMeshCouplingStrategy::none;

    //! Constraint enforcement method
    enum Inpar::CONSTRAINTS::EmbeddedMeshConstraintEnforcement
        embedded_mesh_constraint_enforcement_ =
            Inpar::CONSTRAINTS::EmbeddedMeshConstraintEnforcement::none;

    //! Penalty parameter for coupling enforcement
    double embedded_mesh_constraint_penalty_parameter_ = 0.0;

    //! Shape function for the mortar Lagrange-multiplicators
    enum Inpar::CONSTRAINTS::SolidToSolidMortarShapefunctions embedded_mesh_mortar_shape_function_ =
        Inpar::CONSTRAINTS::SolidToSolidMortarShapefunctions::none;

    //! Nodal Dof set strategy for XFEM
    Core::Geo::Cut::NodalDofSetStrategy xfem_nodal_dof_set_strategy_ =
        Core::Geo::Cut::NodalDofSetStrategy::NDS_Strategy_OneDofset_PerNodeAndPosition;

    //! Integration method for volume cells for XFEM
    Core::Geo::Cut::VCellGaussPts xfem_volume_cell_gauss_point_by_ =
        Core::Geo::Cut::VCellGaussPts::VCellGaussPts_Tessellation;

    //! Integration method for boundary cells for XFEM
    Core::Geo::Cut::BCellGaussPts xfem_bcell_gauss_point_by_ =
        Core::Geo::Cut::BCellGaussPts::BCellGaussPts_Tessellation;

    //! Get gmsh output of cut
    bool gmsh_cut_out_ = false;

    //! Print coutput of cut on the screen
    bool cut_screen_output_ = false;

    //! Parameter list of cut
    Teuchos::ParameterList cut_params_;
  };
}  // namespace CONSTRAINTS::EMBEDDEDMESH

FOUR_C_NAMESPACE_CLOSE

#endif