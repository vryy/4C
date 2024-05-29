/*! \file

\brief A displacement based solid element formulation with FBAR element technology

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_FBAR_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_FBAR_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_cell_type_traits.hpp"
#include "4C_discretization_fem_general_element.hpp"
#include "4C_solid_3D_ele_calc.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_solid_3D_ele_calc_lib_fbar.H"

FOUR_C_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{
  template <CORE::FE::CellType celltype>
  struct FBarPreparationData
  {
    /// jacobian mapping evaluated at element centroid
    JacobianMapping<celltype> jacobian_mapping_centroid;

    /// deformation gradient at element centroid
    SpatialMaterialMapping<celltype> spatial_material_mapping_centroid;
  };

  struct FBarHistoryData
  {
    // no history data needed
  };

  /*!
   * @brief A displacement based solid element formulation with FBAR element technology
   *
   * @tparam celltype
   */
  template <CORE::FE::CellType celltype>
  struct FBarFormulation
  {
    static constexpr bool has_gauss_point_history = false;
    static constexpr bool has_global_history = false;
    static constexpr bool has_preparation_data = true;

    using LinearizationContainer = FBarLinearizationContainer<celltype>;
    using PreparationData = FBarPreparationData<celltype>;


    static FBarPreparationData<celltype> Prepare(
        const CORE::Elements::Element& ele, const ElementNodes<celltype>& nodal_coordinates)
    {
      const JacobianMapping<celltype> jacobian_mapping_centroid =
          EvaluateJacobianMappingCentroid(nodal_coordinates);

      return {jacobian_mapping_centroid,
          EvaluateSpatialMaterialMapping(jacobian_mapping_centroid, nodal_coordinates)};
    }

    template <typename Evaluator>
    static auto Evaluate(const CORE::Elements::Element& ele,
        const ElementNodes<celltype>& nodal_coordinates,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const FBarPreparationData<celltype>& preparation_data, Evaluator evaluator)
    {
      const SpatialMaterialMapping<celltype> spatial_material_mapping =
          EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

      // factor (detF0/detF)^1/3
      const double fbar_factor = EvaluateFbarFactor(
          preparation_data.spatial_material_mapping_centroid.determinant_deformation_gradient_,
          spatial_material_mapping.determinant_deformation_gradient_);

      const FBarLinearizationContainer<celltype> linearization = std::invoke(
          [&]()
          {
            FBarLinearizationContainer<celltype> linearization{};
            linearization.Bop = EvaluateStrainGradient(jacobian_mapping, spatial_material_mapping);

            linearization.Hop = EvaluateFbarHOperator(jacobian_mapping.N_XYZ_,
                preparation_data.jacobian_mapping_centroid.N_XYZ_, spatial_material_mapping,
                preparation_data.spatial_material_mapping_centroid);

            linearization.fbar_factor = fbar_factor;

            linearization.cauchygreen = EvaluateCauchyGreen(spatial_material_mapping);

            return linearization;
          });

      // deformation gradient F_bar and resulting strains: F_bar = (detF_0/detF)^1/3 F
      const SpatialMaterialMapping<celltype> spatial_material_mapping_bar =
          EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates, fbar_factor);

      const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> cauchygreen_bar =
          EvaluateCauchyGreen(spatial_material_mapping_bar);

      CORE::LINALG::Matrix<DETAIL::num_str<celltype>, 1> gl_strain_bar =
          EvaluateGreenLagrangeStrain(cauchygreen_bar);

      return evaluator(
          spatial_material_mapping_bar.deformation_gradient_, gl_strain_bar, linearization);
    }

    static inline CORE::LINALG::Matrix<9, CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_displacements(const CORE::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
            deformation_gradient,
        const FBarPreparationData<celltype>& preparation_data)
    {
      FOUR_C_THROW(
          "This derivative of the deformation gradient w.r.t. the displacements is not "
          "implemented");
    }

    static inline CORE::LINALG::Matrix<9, CORE::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_xi(const CORE::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
            deformation_gradient,
        const FBarPreparationData<celltype>& preparation_data)
    {
      FOUR_C_THROW("This derivative of the deformation gradient w.r.t. xi is not implemented");
    }

    static inline CORE::LINALG::Matrix<9,
        CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype> * CORE::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_displacements_d_xi(const CORE::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
            deformation_gradient,
        const FBarPreparationData<celltype>& preparation_data)
    {
      FOUR_C_THROW(
          "This second derivative of the deformation gradient w.r.t. the displacements and xi is "
          "not implemented");
    }

    static CORE::LINALG::Matrix<DETAILS::num_str<celltype>,
        CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>
    GetLinearBOperator(const FBarLinearizationContainer<celltype>& linearization)
    {
      return linearization.Bop;
    }

    static void add_internal_force_vector(const FBarLinearizationContainer<celltype>& linearization,
        const Stress<celltype>& stress, const double integration_factor,
        const FBarPreparationData<celltype>& preparation_data,
        CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>, 1>&
            force_vector)
    {
      DRT::ELEMENTS::add_internal_force_vector(
          linearization.Bop, stress, integration_factor / linearization.fbar_factor, force_vector);
    }

    static void AddStiffnessMatrix(const FBarLinearizationContainer<celltype>& linearization,
        const JacobianMapping<celltype>& jacobian_mapping, const Stress<celltype>& stress,
        const double integration_factor, const FBarPreparationData<celltype>& preparation_data,
        CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>,
            CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>& stiffness_matrix)
    {
      DRT::ELEMENTS::AddElasticStiffnessMatrix(linearization.Bop, stress,
          integration_factor * linearization.fbar_factor, stiffness_matrix);
      DRT::ELEMENTS::AddGeometricStiffnessMatrix(jacobian_mapping.N_XYZ_, stress,
          integration_factor / linearization.fbar_factor, stiffness_matrix);

      // additional stiffness matrix needed for fbar method
      AddFbarStiffnessMatrix(linearization.Bop, linearization.Hop, linearization.fbar_factor,
          integration_factor, linearization.cauchygreen, stress, stiffness_matrix);
    }
  };

  template <CORE::FE::CellType celltype>
  using FBarSolidIntegrator = SolidEleCalc<celltype, FBarFormulation<celltype>>;


}  // namespace DRT::ELEMENTS

FOUR_C_NAMESPACE_CLOSE
#endif
