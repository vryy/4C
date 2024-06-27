/*! \file

\brief A displacement based solid element formulation with MULF prestressing

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_MULF_FBAR_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_MULF_FBAR_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_solid_3D_ele_calc.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_solid_3D_ele_calc_lib_fbar.H"
#include "4C_solid_3D_ele_calc_lib_io.hpp"
#include "4C_solid_3D_ele_calc_lib_mulf.H"
#include "4C_solid_3D_ele_formulation.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret::ELEMENTS
{
  template <Core::FE::CellType celltype>
  struct MulfFBarPreparationData
  {
    Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_nodes<celltype>> N_XYZ{};

    SpatialMaterialMapping<celltype> spatial_material_mapping{};
  };

  namespace Details
  {
    template <Core::FE::CellType celltype>
    SpatialMaterialMapping<celltype> evaluate_mulf_spatial_material_mapping_centroid(
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions_centeroid,
        const ElementNodes<celltype>& nodal_coordinates,
        const MulfHistoryData<celltype>& mulf_data_centeroid)
    {
      Core::LinAlg::Matrix<3, 3> delta_F = evaluate_mulf_deformation_gradient_update(
          shape_functions_centeroid, nodal_coordinates.displacements_, mulf_data_centeroid);

      SpatialMaterialMapping<celltype> spatial_material_mapping_centeroid{};

      spatial_material_mapping_centeroid.deformation_gradient_.multiply(
          delta_F, mulf_data_centeroid.deformation_gradient);

      spatial_material_mapping_centeroid.inverse_deformation_gradient_ =
          spatial_material_mapping_centeroid.deformation_gradient_;

      spatial_material_mapping_centeroid.determinant_deformation_gradient_ =
          spatial_material_mapping_centeroid.inverse_deformation_gradient_.invert();

      return spatial_material_mapping_centeroid;
    }

    template <Core::FE::CellType celltype>
    SpatialMaterialMapping<celltype> GetSpatialMaterialMappingBar(
        SpatialMaterialMapping<celltype> spatial_material_mapping, const double fbar_factor)
    {
      spatial_material_mapping.deformation_gradient_.scale(fbar_factor);
      spatial_material_mapping.determinant_deformation_gradient_ *=
          Core::FE::dim<celltype> * fbar_factor;
      spatial_material_mapping.inverse_deformation_gradient_.scale(1.0 / fbar_factor);

      return spatial_material_mapping;
    }

    /*!
     * @brief Do a MULF update step on the mulf data
     */
    template <Core::FE::CellType celltype>
    void UpdateMulfHistory(const ElementNodes<celltype>& element_nodes,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        MulfHistoryData<celltype>& mulf_data)
    {
      Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>> delta_defgrd =
          evaluate_mulf_deformation_gradient_update(
              shape_functions, element_nodes.displacements_, mulf_data);

      Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>> inv_delta_defgrd(true);
      inv_delta_defgrd.invert(delta_defgrd);

      Core::LinAlg::Matrix<3, 3> old_defgrd = mulf_data.deformation_gradient;
      mulf_data.deformation_gradient.multiply(delta_defgrd, old_defgrd);

      Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>> invJ_new;
      invJ_new.multiply_tn(inv_delta_defgrd, mulf_data.inverse_jacobian);
      mulf_data.inverse_jacobian = std::move(invJ_new);
    }
  }  // namespace Details



  /*!
   * @brief A solid element formulation with MULF prestressing and F-Bar
   *
   * @tparam celltype
   */
  template <Core::FE::CellType celltype>
  struct MulfFBarFormulation
  {
    static constexpr bool has_gauss_point_history = true;
    static constexpr bool has_global_history = true;
    static constexpr bool has_preparation_data = true;
    static constexpr bool is_prestress_updatable = true;

    using LinearizationContainer = FBarLinearizationContainer<celltype>;
    using GaussPointHistory = MulfHistoryData<celltype>;
    using GlobalHistory = MulfHistoryData<celltype>;
    using PreparationData = MulfFBarPreparationData<celltype>;

    static MulfFBarPreparationData<celltype> Prepare(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& nodal_coordinates, GlobalHistory& global_history)
    {
      if (!global_history.is_setup)
      {
        const JacobianMapping<celltype> jacobian_mapping =
            evaluate_jacobian_mapping_centroid(nodal_coordinates);

        global_history.inverse_jacobian = jacobian_mapping.inverse_jacobian_;
        global_history.is_setup = true;
      }

      // set coordinates in parameter space at centroid as zero -> xi = [0; 0; 0]
      Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> xi_centroid =
          EvaluateParameterCoordinateCentroid<celltype>();

      // shape functions and derivatives evaluated at element centroid
      const ShapeFunctionsAndDerivatives<celltype> shape_functions_centeroid =
          EvaluateShapeFunctionsAndDerivs<celltype>(xi_centroid, nodal_coordinates);

      JacobianMapping<celltype> jacobian_mapping =
          evaluate_jacobian_mapping_centroid(nodal_coordinates);

      Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::num_nodes<celltype>> N_XYZ_0;
      N_XYZ_0.multiply(jacobian_mapping.inverse_jacobian_, shape_functions_centeroid.derivatives_);

      return {N_XYZ_0, Details::evaluate_mulf_spatial_material_mapping_centroid(
                           shape_functions_centeroid, nodal_coordinates, global_history)};
    }

    template <typename Evaluator>
    static auto evaluate(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const MulfFBarPreparationData<celltype>& mapping_center,
        const MulfHistoryData<celltype>& mulf_data_centeroid,
        MulfHistoryData<celltype>& mulf_data_gp, Evaluator evaluator)
    {
      if (!mulf_data_gp.is_setup)
      {
        mulf_data_gp.inverse_jacobian = jacobian_mapping.inverse_jacobian_;
        mulf_data_gp.is_setup = true;
      }

      const SpatialMaterialMapping<celltype> spatial_material_mapping =
          evaluate_mulf_spatial_material_mapping(
              jacobian_mapping, shape_functions, element_nodes.displacements_, mulf_data_gp);

      const double fbar_factor = evaluate_fbar_factor(
          mapping_center.spatial_material_mapping.determinant_deformation_gradient_,
          spatial_material_mapping.determinant_deformation_gradient_);

      const FBarLinearizationContainer<celltype> linearization = std::invoke(
          [&]()
          {
            FBarLinearizationContainer<celltype> linearization{};
            linearization.Bop =
                evaluate_strain_gradient(jacobian_mapping, spatial_material_mapping);

            linearization.Hop =
                evaluate_fbar_h_operator(jacobian_mapping.N_XYZ_, mapping_center.N_XYZ,
                    spatial_material_mapping, mapping_center.spatial_material_mapping);

            linearization.fbar_factor = fbar_factor;

            linearization.cauchygreen = evaluate_cauchy_green(spatial_material_mapping);

            return linearization;
          });

      const SpatialMaterialMapping<celltype> spatial_material_mapping_bar =
          Details::GetSpatialMaterialMappingBar(spatial_material_mapping, fbar_factor);

      const Core::LinAlg::Matrix<Core::FE::dim<celltype>, Core::FE::dim<celltype>> cauchygreen_bar =
          evaluate_cauchy_green<celltype>(spatial_material_mapping_bar);

      const Core::LinAlg::Matrix<Details::num_str<celltype>, 1> gl_strain_bar =
          evaluate_green_lagrange_strain(cauchygreen_bar);

      return evaluator(
          spatial_material_mapping_bar.deformation_gradient_, gl_strain_bar, linearization);
    }


    static Core::LinAlg::Matrix<Details::num_str<celltype>,
        Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>
    get_linear_b_operator(const FBarLinearizationContainer<celltype>& linearization)
    {
      return linearization.Bop;
    }

    static void add_internal_force_vector(const FBarLinearizationContainer<celltype>& linearization,
        const Stress<celltype>& stress, const double integration_factor,
        const MulfFBarPreparationData<celltype>& mapping_center,
        MulfHistoryData<celltype>& mulf_data_centeroid, MulfHistoryData<celltype>& mulf_data_gp,
        Core::LinAlg::Matrix<Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>, 1>&
            force_vector)
    {
      Discret::ELEMENTS::add_internal_force_vector(
          linearization.Bop, stress, integration_factor / linearization.fbar_factor, force_vector);
    }

    static void add_stiffness_matrix(const FBarLinearizationContainer<celltype>& linearization,
        const JacobianMapping<celltype>& jacobian_mapping, const Stress<celltype>& stress,
        const double integration_factor, const MulfFBarPreparationData<celltype>& mapping_center,
        MulfHistoryData<celltype>& mulf_data_centeroid, MulfHistoryData<celltype>& mulf_data_gp,
        Core::LinAlg::Matrix<Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>,
            Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>& stiffness_matrix)
    {
      Discret::ELEMENTS::add_elastic_stiffness_matrix(linearization.Bop, stress,
          integration_factor * linearization.fbar_factor, stiffness_matrix);
      Discret::ELEMENTS::add_geometric_stiffness_matrix(jacobian_mapping.N_XYZ_, stress,
          integration_factor / linearization.fbar_factor, stiffness_matrix);

      // additional stiffness matrix needed for fbar method
      add_fbar_stiffness_matrix(linearization.Bop, linearization.Hop, linearization.fbar_factor,
          integration_factor, linearization.cauchygreen, stress, stiffness_matrix);
    }

    static void pack(
        const MulfHistoryData<celltype>& history_data, Core::Communication::PackBuffer& data)
    {
      Core::Communication::ParObject::add_to_pack(data, history_data.inverse_jacobian);
      Core::Communication::ParObject::add_to_pack(data, history_data.deformation_gradient);
      Core::Communication::ParObject::add_to_pack(data, history_data.is_setup);
    }

    static void unpack(std::vector<char>::size_type& position, const std::vector<char>& data,
        MulfHistoryData<celltype>& history_data)
    {
      Core::Communication::ParObject::extract_from_pack(
          position, data, history_data.inverse_jacobian);
      Core::Communication::ParObject::extract_from_pack(
          position, data, history_data.deformation_gradient);
      int is_setup_int;
      Core::Communication::ParObject::extract_from_pack(position, data, is_setup_int);
      history_data.is_setup = static_cast<bool>(is_setup_int);
    }

    static inline void update_prestress(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const MulfFBarPreparationData<celltype>& mapping_center,
        MulfHistoryData<celltype>& mulf_data_centeroid)
    {
      Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1> xi_centroid =
          EvaluateParameterCoordinateCentroid<celltype>();

      ShapeFunctionsAndDerivatives<celltype> shape_functions_centeroid =
          EvaluateShapeFunctionsAndDerivs<celltype>(xi_centroid, element_nodes);

      Details::UpdateMulfHistory(element_nodes, shape_functions_centeroid, mulf_data_centeroid);
    }

    static inline void update_prestress(const Core::Elements::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
            deformation_gradient,
        const MulfFBarPreparationData<celltype>& mapping_center,
        MulfHistoryData<celltype>& mulf_data_centeroid, MulfHistoryData<celltype>& mulf_data_gp)
    {
      Details::UpdateMulfHistory(element_nodes, shape_functions, mulf_data_gp);
    }
  };

  template <Core::FE::CellType celltype>
  using MulfFBarSolidIntegrator = SolidEleCalc<celltype, MulfFBarFormulation<celltype>>;



}  // namespace Discret::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif
