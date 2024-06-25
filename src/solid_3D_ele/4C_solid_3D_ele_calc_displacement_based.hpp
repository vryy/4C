/*! \file

\brief A displacement based solid element formulation

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_DISPLACEMENT_BASED_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_DISPLACEMENT_BASED_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_solid_3D_ele_calc.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret::ELEMENTS
{
  template <Core::FE::CellType celltype>
  struct DisplacementBasedLinearizationContainer
  {
    Core::LinAlg::Matrix<Details::num_str<celltype>,
        Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>
        Bop_{};
  };

  /*!
   * @brief A displacement based solid element formulation
   *
   * @tparam celltype
   */
  template <Core::FE::CellType celltype>
  struct DisplacementBasedFormulation
  {
    static constexpr bool has_gauss_point_history = false;
    static constexpr bool has_global_history = false;
    static constexpr bool has_preparation_data = false;

    using LinearizationContainer = DisplacementBasedLinearizationContainer<celltype>;

    template <typename Evaluator>
    static inline auto evaluate(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& nodal_coordinates,
        const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping, Evaluator evaluator)
    {
      const SpatialMaterialMapping<celltype> spatial_material_mapping =
          evaluate_spatial_material_mapping(jacobian_mapping, nodal_coordinates);

      const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>> cauchygreen =
          evaluate_cauchy_green(spatial_material_mapping);

      const Core::LinAlg::Matrix<DETAIL::num_str<celltype>, 1> gl_strain =
          evaluate_green_lagrange_strain(cauchygreen);

      const DisplacementBasedLinearizationContainer<celltype> linearization = std::invoke(
          [&]()
          {
            DisplacementBasedLinearizationContainer<celltype> linearization{};
            linearization.Bop_ =
                evaluate_strain_gradient(jacobian_mapping, spatial_material_mapping);

            return linearization;
          });

      return evaluator(spatial_material_mapping.deformation_gradient_, gl_strain, linearization);
    }

    static inline Core::LinAlg::Matrix<9, Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_displacements(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
            deformation_gradient)
    {
      Core::LinAlg::Matrix<9, Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>> d_F_dd{};

      // evaluate derivative w.r.t. displacements
      for (int i = 0; i < Core::FE::num_nodes<celltype>; ++i)
      {
        d_F_dd(0, Core::FE::dim<celltype> * i + 0) = jacobian_mapping.N_XYZ_(0, i);
        d_F_dd(1, Core::FE::dim<celltype> * i + 1) = jacobian_mapping.N_XYZ_(1, i);
        d_F_dd(2, Core::FE::dim<celltype> * i + 2) = jacobian_mapping.N_XYZ_(2, i);
        d_F_dd(3, Core::FE::dim<celltype> * i + 0) = jacobian_mapping.N_XYZ_(1, i);
        d_F_dd(4, Core::FE::dim<celltype> * i + 1) = jacobian_mapping.N_XYZ_(2, i);
        d_F_dd(5, Core::FE::dim<celltype> * i + 0) = jacobian_mapping.N_XYZ_(2, i);
        d_F_dd(6, Core::FE::dim<celltype> * i + 1) = jacobian_mapping.N_XYZ_(0, i);
        d_F_dd(7, Core::FE::dim<celltype> * i + 2) = jacobian_mapping.N_XYZ_(1, i);
        d_F_dd(8, Core::FE::dim<celltype> * i + 2) = jacobian_mapping.N_XYZ_(0, i);
      }

      return d_F_dd;
    }

    static inline Core::LinAlg::Matrix<9, Core::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_xi(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
            deformation_gradient)
    {
      Core::LinAlg::Matrix<9, Core::FE::dim<celltype>> d_F_dxi{};

      Core::LinAlg::Matrix<Core::FE::num_nodes<celltype>, Core::FE::dim<celltype>> xXF(true);
      Core::LinAlg::Matrix<Core::FE::dim<celltype>,
          Core::FE::DisTypeToNumDeriv2<celltype>::numderiv2>
          xXFsec(true);
      xXF.update(1.0, element_nodes.reference_coordinates_, 0.0);
      xXF.update(1.0, element_nodes.displacements_, 1.0);
      xXF.multiply_nt(-1.0, element_nodes.reference_coordinates_, deformation_gradient, 1.0);

      Core::LinAlg::Matrix<Core::FE::DisTypeToNumDeriv2<celltype>::numderiv2,
          Core::FE::num_nodes<celltype>>
          deriv2(true);
      Core::FE::shape_function_deriv2<celltype>(xi, deriv2);

      xXFsec.multiply_tt(1.0, xXF, deriv2, 0.0);

      for (int a = 0; a < Core::FE::dim<celltype>; ++a)
      {
        for (int b = 0; b < Core::FE::dim<celltype>; ++b)
        {
          using VoigtMapping = Core::LinAlg::Voigt::IndexMappings;
          d_F_dxi(VoigtMapping::non_symmetric_tensor_to_voigt9_index(a, b), 0) +=
              xXFsec(a, 0) * jacobian_mapping.inverse_jacobian_(b, 0) +
              xXFsec(a, 3) * jacobian_mapping.inverse_jacobian_(b, 1) +
              xXFsec(a, 4) * jacobian_mapping.inverse_jacobian_(b, 2);
          d_F_dxi(VoigtMapping::non_symmetric_tensor_to_voigt9_index(a, b), 1) +=
              xXFsec(a, 3) * jacobian_mapping.inverse_jacobian_(b, 0) +
              xXFsec(a, 1) * jacobian_mapping.inverse_jacobian_(b, 1) +
              xXFsec(a, 5) * jacobian_mapping.inverse_jacobian_(b, 2);
          d_F_dxi(VoigtMapping::non_symmetric_tensor_to_voigt9_index(a, b), 2) +=
              xXFsec(a, 4) * jacobian_mapping.inverse_jacobian_(b, 0) +
              xXFsec(a, 5) * jacobian_mapping.inverse_jacobian_(b, 1) +
              xXFsec(a, 2) * jacobian_mapping.inverse_jacobian_(b, 2);
        }
      }

      return d_F_dxi;
    }

    static inline Core::LinAlg::Matrix<9,
        Core::FE::num_nodes<celltype> * Core::FE::dim<celltype> * Core::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_displacements_d_xi(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
            deformation_gradient)
    {
      Core::LinAlg::Matrix<9,
          Core::FE::num_nodes<celltype> * Core::FE::dim<celltype> * Core::FE::dim<celltype>>
          d2_F_dxi_dd{};

      // evaluate derivative w.r.t. displacements
      Core::LinAlg::Matrix<Core::FE::DisTypeToNumDeriv2<celltype>::numderiv2,
          Core::FE::dim<celltype>>
          Xsec(true);
      Core::LinAlg::Matrix<Core::FE::num_nodes<celltype>,
          Core::FE::DisTypeToNumDeriv2<celltype>::numderiv2>
          N_XYZ_Xsec(true);

      Core::LinAlg::Matrix<Core::FE::DisTypeToNumDeriv2<celltype>::numderiv2,
          Core::FE::num_nodes<celltype>>
          deriv2(true);
      Core::FE::shape_function_deriv2<celltype>(xi, deriv2);
      Xsec.multiply(1.0, deriv2, element_nodes.reference_coordinates_, 0.0);
      N_XYZ_Xsec.multiply_tt(1.0, jacobian_mapping.N_XYZ_, Xsec, 0.0);

      for (int i = 0; i < Core::FE::dim<celltype>; ++i)
      {
        for (int j = 0; j < Core::FE::dim<celltype>; ++j)
        {
          for (int k = 0; k < Core::FE::num_nodes<celltype>; ++k)
          {
            using VoigtMapping = Core::LinAlg::Voigt::IndexMappings;
            d2_F_dxi_dd(VoigtMapping::non_symmetric_tensor_to_voigt9_index(i, j),
                Core::FE::dim<celltype> * (Core::FE::dim<celltype> * k + i) + 0) +=
                deriv2(0, k) * jacobian_mapping.inverse_jacobian_(j, 0) +
                deriv2(3, k) * jacobian_mapping.inverse_jacobian_(j, 1) +
                deriv2(4, k) * jacobian_mapping.inverse_jacobian_(j, 2) -
                N_XYZ_Xsec(k, 0) * jacobian_mapping.inverse_jacobian_(j, 0) -
                N_XYZ_Xsec(k, 3) * jacobian_mapping.inverse_jacobian_(j, 1) -
                N_XYZ_Xsec(k, 4) * jacobian_mapping.inverse_jacobian_(j, 2);

            d2_F_dxi_dd(VoigtMapping::non_symmetric_tensor_to_voigt9_index(i, j),
                Core::FE::dim<celltype> * (Core::FE::dim<celltype> * k + i) + 1) +=
                deriv2(3, k) * jacobian_mapping.inverse_jacobian_(j, 0) +
                deriv2(1, k) * jacobian_mapping.inverse_jacobian_(j, 1) +
                deriv2(5, k) * jacobian_mapping.inverse_jacobian_(j, 2) -
                N_XYZ_Xsec(k, 3) * jacobian_mapping.inverse_jacobian_(j, 0) -
                N_XYZ_Xsec(k, 1) * jacobian_mapping.inverse_jacobian_(j, 1) -
                N_XYZ_Xsec(k, 5) * jacobian_mapping.inverse_jacobian_(j, 2);

            d2_F_dxi_dd(VoigtMapping::non_symmetric_tensor_to_voigt9_index(i, j),
                Core::FE::dim<celltype> * (Core::FE::dim<celltype> * k + i) + 2) +=
                deriv2(4, k) * jacobian_mapping.inverse_jacobian_(j, 0) +
                deriv2(5, k) * jacobian_mapping.inverse_jacobian_(j, 1) +
                deriv2(2, k) * jacobian_mapping.inverse_jacobian_(j, 2) -
                N_XYZ_Xsec(k, 4) * jacobian_mapping.inverse_jacobian_(j, 0) -
                N_XYZ_Xsec(k, 5) * jacobian_mapping.inverse_jacobian_(j, 1) -
                N_XYZ_Xsec(k, 2) * jacobian_mapping.inverse_jacobian_(j, 2);
          }
        }
      }

      return d2_F_dxi_dd;
    }

    static Core::LinAlg::Matrix<Details::num_str<celltype>,
        Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>
    get_linear_b_operator(const DisplacementBasedLinearizationContainer<celltype>& linearization)
    {
      return linearization.Bop_;
    }

    static void add_internal_force_vector(
        const DisplacementBasedLinearizationContainer<celltype>& linearization,
        const Stress<celltype>& stress, const double integration_factor,
        Core::LinAlg::Matrix<Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>, 1>&
            force_vector)
    {
      Discret::ELEMENTS::add_internal_force_vector(
          linearization.Bop_, stress, integration_factor, force_vector);
    }

    static void add_stiffness_matrix(
        const DisplacementBasedLinearizationContainer<celltype>& linearization,
        const JacobianMapping<celltype>& jacobian_mapping, const Stress<celltype>& stress,
        const double integration_factor,
        Core::LinAlg::Matrix<Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>,
            Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>& stiffness_matrix)
    {
      Discret::ELEMENTS::add_elastic_stiffness_matrix(
          linearization.Bop_, stress, integration_factor, stiffness_matrix);
      Discret::ELEMENTS::add_geometric_stiffness_matrix(
          jacobian_mapping.N_XYZ_, stress, integration_factor, stiffness_matrix);
    }
  };

  template <Core::FE::CellType celltype>
  using DisplacementBasedSolidIntegrator =
      SolidEleCalc<celltype, DisplacementBasedFormulation<celltype>>;


}  // namespace Discret::ELEMENTS

FOUR_C_NAMESPACE_CLOSE
#endif
