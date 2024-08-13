/*! \file

\brief Utility functions for MULF element technology

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_LIB_MULF_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_LIB_MULF_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_fixedsizematrix_generators.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret::ELEMENTS
{
  /*!
   * @brief A container storing the history data for element implementations with MULF prestressing
   */
  template <Core::FE::CellType celltype>
  struct MulfHistoryData
  {
    Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>> inverse_jacobian =
        Core::LinAlg::IdentityMatrix<Core::FE::dim<celltype>>();
    Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>> deformation_gradient =
        Core::LinAlg::IdentityMatrix<Core::FE::dim<celltype>>();
    bool is_setup = false;
  };

  /*!
   * @brief Evalaute the update of the deformation gradient for MULF prestressing
   */
  template <Core::FE::CellType celltype>
  Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>>
  evaluate_mulf_deformation_gradient_update(
      const Discret::ELEMENTS::ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const Core::LinAlg::Matrix<Core::FE::num_nodes<celltype>, Core::FE::dim<celltype>>&
          nodal_displacements,
      const Discret::ELEMENTS::MulfHistoryData<celltype>& mulf_history_data)
  {
    Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::num_nodes<celltype>> N_xyz;

    N_xyz.multiply(mulf_history_data.inverse_jacobian, shape_functions.derivatives_);

    Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>> defgrd =
        Core::LinAlg::IdentityMatrix<Core::FE::dim<celltype>>();

    defgrd.multiply_tt(1.0, nodal_displacements, N_xyz, 1.0);

    return defgrd;
  }

  /*!
   * @brief Evaluate the spatial material mapping (deformation gradient) for MULF prestressing
   */
  template <Core::FE::CellType celltype>
  Discret::ELEMENTS::SpatialMaterialMapping<celltype> evaluate_mulf_spatial_material_mapping(
      const Discret::ELEMENTS::JacobianMapping<celltype>& jacobian_mapping,
      const Discret::ELEMENTS::ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const Core::LinAlg::Matrix<Core::FE::num_nodes<celltype>, Core::FE::dim<celltype>>&
          nodal_displacements,
      const Discret::ELEMENTS::MulfHistoryData<celltype>& mulf_history_data)
  {
    Discret::ELEMENTS::SpatialMaterialMapping<celltype> spatial_material_mapping;

    Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>> defgrd =
        evaluate_mulf_deformation_gradient_update(
            shape_functions, nodal_displacements, mulf_history_data);

    spatial_material_mapping.deformation_gradient_.multiply(
        defgrd, mulf_history_data.deformation_gradient);

    spatial_material_mapping.inverse_deformation_gradient_ =
        spatial_material_mapping.deformation_gradient_;
    spatial_material_mapping.determinant_deformation_gradient_ =
        spatial_material_mapping.inverse_deformation_gradient_.invert();

    return spatial_material_mapping;
  }
}  // namespace Discret::ELEMENTS

FOUR_C_NAMESPACE_CLOSE
#endif