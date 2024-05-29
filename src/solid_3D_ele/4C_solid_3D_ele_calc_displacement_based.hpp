/*! \file

\brief A displacement based solid element formulation

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_DISPLACEMENT_BASED_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_DISPLACEMENT_BASED_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_cell_type_traits.hpp"
#include "4C_discretization_fem_general_element.hpp"
#include "4C_solid_3D_ele_calc.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{
  template <CORE::FE::CellType celltype>
  struct DisplacementBasedLinearizationContainer
  {
    CORE::LINALG::Matrix<DETAILS::num_str<celltype>,
        CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>
        Bop_{};
  };

  /*!
   * @brief A displacement based solid element formulation
   *
   * @tparam celltype
   */
  template <CORE::FE::CellType celltype>
  struct DisplacementBasedFormulation
  {
    static constexpr bool has_gauss_point_history = false;
    static constexpr bool has_global_history = false;
    static constexpr bool has_preparation_data = false;

    using LinearizationContainer = DisplacementBasedLinearizationContainer<celltype>;

    template <typename Evaluator>
    static inline auto Evaluate(const CORE::Elements::Element& ele,
        const ElementNodes<celltype>& nodal_coordinates,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping, Evaluator evaluator)
    {
      const SpatialMaterialMapping<celltype> spatial_material_mapping =
          EvaluateSpatialMaterialMapping(jacobian_mapping, nodal_coordinates);

      const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>> cauchygreen =
          EvaluateCauchyGreen(spatial_material_mapping);

      const CORE::LINALG::Matrix<DETAIL::num_str<celltype>, 1> gl_strain =
          EvaluateGreenLagrangeStrain(cauchygreen);

      const DisplacementBasedLinearizationContainer<celltype> linearization = std::invoke(
          [&]()
          {
            DisplacementBasedLinearizationContainer<celltype> linearization{};
            linearization.Bop_ = EvaluateStrainGradient(jacobian_mapping, spatial_material_mapping);

            return linearization;
          });

      return evaluator(spatial_material_mapping.deformation_gradient_, gl_strain, linearization);
    }

    static inline CORE::LINALG::Matrix<9, CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_displacements(const CORE::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
            deformation_gradient)
    {
      CORE::LINALG::Matrix<9, CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>> d_F_dd{};

      // evaluate derivative w.r.t. displacements
      for (int i = 0; i < CORE::FE::num_nodes<celltype>; ++i)
      {
        d_F_dd(0, CORE::FE::dim<celltype> * i + 0) = jacobian_mapping.N_XYZ_(0, i);
        d_F_dd(1, CORE::FE::dim<celltype> * i + 1) = jacobian_mapping.N_XYZ_(1, i);
        d_F_dd(2, CORE::FE::dim<celltype> * i + 2) = jacobian_mapping.N_XYZ_(2, i);
        d_F_dd(3, CORE::FE::dim<celltype> * i + 0) = jacobian_mapping.N_XYZ_(1, i);
        d_F_dd(4, CORE::FE::dim<celltype> * i + 1) = jacobian_mapping.N_XYZ_(2, i);
        d_F_dd(5, CORE::FE::dim<celltype> * i + 0) = jacobian_mapping.N_XYZ_(2, i);
        d_F_dd(6, CORE::FE::dim<celltype> * i + 1) = jacobian_mapping.N_XYZ_(0, i);
        d_F_dd(7, CORE::FE::dim<celltype> * i + 2) = jacobian_mapping.N_XYZ_(1, i);
        d_F_dd(8, CORE::FE::dim<celltype> * i + 2) = jacobian_mapping.N_XYZ_(0, i);
      }

      return d_F_dd;
    }

    static inline CORE::LINALG::Matrix<9, CORE::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_xi(const CORE::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
            deformation_gradient)
    {
      CORE::LINALG::Matrix<9, CORE::FE::dim<celltype>> d_F_dxi{};

      CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype>, CORE::FE::dim<celltype>> xXF(true);
      CORE::LINALG::Matrix<CORE::FE::dim<celltype>,
          CORE::FE::DisTypeToNumDeriv2<celltype>::numderiv2>
          xXFsec(true);
      xXF.Update(1.0, element_nodes.reference_coordinates_, 0.0);
      xXF.Update(1.0, element_nodes.displacements_, 1.0);
      xXF.MultiplyNT(-1.0, element_nodes.reference_coordinates_, deformation_gradient, 1.0);

      CORE::LINALG::Matrix<CORE::FE::DisTypeToNumDeriv2<celltype>::numderiv2,
          CORE::FE::num_nodes<celltype>>
          deriv2(true);
      CORE::FE::shape_function_deriv2<celltype>(xi, deriv2);

      xXFsec.MultiplyTT(1.0, xXF, deriv2, 0.0);

      for (int a = 0; a < CORE::FE::dim<celltype>; ++a)
      {
        for (int b = 0; b < CORE::FE::dim<celltype>; ++b)
        {
          using VoigtMapping = CORE::LINALG::VOIGT::IndexMappings;
          d_F_dxi(VoigtMapping::NonSymToVoigt9(a, b), 0) +=
              xXFsec(a, 0) * jacobian_mapping.inverse_jacobian_(b, 0) +
              xXFsec(a, 3) * jacobian_mapping.inverse_jacobian_(b, 1) +
              xXFsec(a, 4) * jacobian_mapping.inverse_jacobian_(b, 2);
          d_F_dxi(VoigtMapping::NonSymToVoigt9(a, b), 1) +=
              xXFsec(a, 3) * jacobian_mapping.inverse_jacobian_(b, 0) +
              xXFsec(a, 1) * jacobian_mapping.inverse_jacobian_(b, 1) +
              xXFsec(a, 5) * jacobian_mapping.inverse_jacobian_(b, 2);
          d_F_dxi(VoigtMapping::NonSymToVoigt9(a, b), 2) +=
              xXFsec(a, 4) * jacobian_mapping.inverse_jacobian_(b, 0) +
              xXFsec(a, 5) * jacobian_mapping.inverse_jacobian_(b, 1) +
              xXFsec(a, 2) * jacobian_mapping.inverse_jacobian_(b, 2);
        }
      }

      return d_F_dxi;
    }

    static inline CORE::LINALG::Matrix<9,
        CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype> * CORE::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_displacements_d_xi(const CORE::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
            deformation_gradient)
    {
      CORE::LINALG::Matrix<9,
          CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype> * CORE::FE::dim<celltype>>
          d2_F_dxi_dd{};

      // evaluate derivative w.r.t. displacements
      CORE::LINALG::Matrix<CORE::FE::DisTypeToNumDeriv2<celltype>::numderiv2,
          CORE::FE::dim<celltype>>
          Xsec(true);
      CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype>,
          CORE::FE::DisTypeToNumDeriv2<celltype>::numderiv2>
          N_XYZ_Xsec(true);

      CORE::LINALG::Matrix<CORE::FE::DisTypeToNumDeriv2<celltype>::numderiv2,
          CORE::FE::num_nodes<celltype>>
          deriv2(true);
      CORE::FE::shape_function_deriv2<celltype>(xi, deriv2);
      Xsec.Multiply(1.0, deriv2, element_nodes.reference_coordinates_, 0.0);
      N_XYZ_Xsec.MultiplyTT(1.0, jacobian_mapping.N_XYZ_, Xsec, 0.0);

      for (int i = 0; i < CORE::FE::dim<celltype>; ++i)
      {
        for (int j = 0; j < CORE::FE::dim<celltype>; ++j)
        {
          for (int k = 0; k < CORE::FE::num_nodes<celltype>; ++k)
          {
            using VoigtMapping = CORE::LINALG::VOIGT::IndexMappings;
            d2_F_dxi_dd(VoigtMapping::NonSymToVoigt9(i, j),
                CORE::FE::dim<celltype> * (CORE::FE::dim<celltype> * k + i) + 0) +=
                deriv2(0, k) * jacobian_mapping.inverse_jacobian_(j, 0) +
                deriv2(3, k) * jacobian_mapping.inverse_jacobian_(j, 1) +
                deriv2(4, k) * jacobian_mapping.inverse_jacobian_(j, 2) -
                N_XYZ_Xsec(k, 0) * jacobian_mapping.inverse_jacobian_(j, 0) -
                N_XYZ_Xsec(k, 3) * jacobian_mapping.inverse_jacobian_(j, 1) -
                N_XYZ_Xsec(k, 4) * jacobian_mapping.inverse_jacobian_(j, 2);

            d2_F_dxi_dd(VoigtMapping::NonSymToVoigt9(i, j),
                CORE::FE::dim<celltype> * (CORE::FE::dim<celltype> * k + i) + 1) +=
                deriv2(3, k) * jacobian_mapping.inverse_jacobian_(j, 0) +
                deriv2(1, k) * jacobian_mapping.inverse_jacobian_(j, 1) +
                deriv2(5, k) * jacobian_mapping.inverse_jacobian_(j, 2) -
                N_XYZ_Xsec(k, 3) * jacobian_mapping.inverse_jacobian_(j, 0) -
                N_XYZ_Xsec(k, 1) * jacobian_mapping.inverse_jacobian_(j, 1) -
                N_XYZ_Xsec(k, 5) * jacobian_mapping.inverse_jacobian_(j, 2);

            d2_F_dxi_dd(VoigtMapping::NonSymToVoigt9(i, j),
                CORE::FE::dim<celltype> * (CORE::FE::dim<celltype> * k + i) + 2) +=
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

    static CORE::LINALG::Matrix<DETAILS::num_str<celltype>,
        CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>
    GetLinearBOperator(const DisplacementBasedLinearizationContainer<celltype>& linearization)
    {
      return linearization.Bop_;
    }

    static void add_internal_force_vector(
        const DisplacementBasedLinearizationContainer<celltype>& linearization,
        const Stress<celltype>& stress, const double integration_factor,
        CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>, 1>&
            force_vector)
    {
      DRT::ELEMENTS::add_internal_force_vector(
          linearization.Bop_, stress, integration_factor, force_vector);
    }

    static void AddStiffnessMatrix(
        const DisplacementBasedLinearizationContainer<celltype>& linearization,
        const JacobianMapping<celltype>& jacobian_mapping, const Stress<celltype>& stress,
        const double integration_factor,
        CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>,
            CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>& stiffness_matrix)
    {
      DRT::ELEMENTS::AddElasticStiffnessMatrix(
          linearization.Bop_, stress, integration_factor, stiffness_matrix);
      DRT::ELEMENTS::AddGeometricStiffnessMatrix(
          jacobian_mapping.N_XYZ_, stress, integration_factor, stiffness_matrix);
    }
  };

  template <CORE::FE::CellType celltype>
  using DisplacementBasedSolidIntegrator =
      SolidEleCalc<celltype, DisplacementBasedFormulation<celltype>>;


}  // namespace DRT::ELEMENTS

FOUR_C_NAMESPACE_CLOSE
#endif
