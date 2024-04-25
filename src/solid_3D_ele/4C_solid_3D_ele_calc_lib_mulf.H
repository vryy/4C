/*! \file

\brief Utility functions for MULF element technology

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_LIB_MULF_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_LIB_MULF_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_cell_type.hpp"
#include "4C_discretization_fem_general_cell_type_traits.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_fixedsizematrix_generators.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{
  /*!
   * @brief A container storing the history data for element implementations with MULF prestressing
   */
  template <CORE::FE::CellType celltype>
  struct MulfHistoryData
  {
    CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> inverse_jacobian =
        CORE::LINALG::IdentityMatrix<CORE::FE::dim<celltype>>();
    CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> deformation_gradient =
        CORE::LINALG::IdentityMatrix<CORE::FE::dim<celltype>>();
    bool is_setup = false;
  };

  /*!
   * @brief Evalaute the update of the deformation gradient for MULF prestressing
   */
  template <CORE::FE::CellType celltype>
  CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>>
  EvaluateMulfDeformationGradientUpdate(
      const DRT::ELEMENTS::ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype>, CORE::FE::dim<celltype>>&
          nodal_displacements,
      const DRT::ELEMENTS::MulfHistoryData<celltype>& mulf_history_data)
  {
    CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::num_nodes<celltype>> N_xyz;

    N_xyz.Multiply(mulf_history_data.inverse_jacobian, shape_functions.derivatives_);

    CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> defgrd =
        CORE::LINALG::IdentityMatrix<CORE::FE::dim<celltype>>();

    defgrd.MultiplyTT(1.0, nodal_displacements, N_xyz, 1.0);

    return defgrd;
  }

  /*!
   * @brief Evaluate the spatial material mapping (deformation gradient) for MULF prestressing
   */
  template <CORE::FE::CellType celltype>
  DRT::ELEMENTS::SpatialMaterialMapping<celltype> EvaluateMulfSpatialMaterialMapping(
      const DRT::ELEMENTS::JacobianMapping<celltype>& jacobian_mapping,
      const DRT::ELEMENTS::ShapeFunctionsAndDerivatives<celltype>& shape_functions,
      const CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype>, CORE::FE::dim<celltype>>&
          nodal_displacements,
      const DRT::ELEMENTS::MulfHistoryData<celltype>& mulf_history_data)
  {
    DRT::ELEMENTS::SpatialMaterialMapping<celltype> spatial_material_mapping;

    CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> defgrd =
        EvaluateMulfDeformationGradientUpdate(
            shape_functions, nodal_displacements, mulf_history_data);

    spatial_material_mapping.deformation_gradient_.Multiply(
        defgrd, mulf_history_data.deformation_gradient);

    spatial_material_mapping.inverse_deformation_gradient_ =
        spatial_material_mapping.deformation_gradient_;
    spatial_material_mapping.determinant_deformation_gradient_ =
        spatial_material_mapping.inverse_deformation_gradient_.Invert();

    return spatial_material_mapping;
  }
}  // namespace DRT::ELEMENTS

FOUR_C_NAMESPACE_CLOSE
#endif