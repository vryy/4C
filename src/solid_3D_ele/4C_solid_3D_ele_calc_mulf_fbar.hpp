/*! \file

\brief A displacement based solid element formulation with MULF prestressing

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_MULF_FBAR_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_MULF_FBAR_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_cell_type.hpp"
#include "4C_discretization_fem_general_cell_type_traits.hpp"
#include "4C_lib_element.hpp"
#include "4C_solid_3D_ele_calc.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_solid_3D_ele_calc_lib_fbar.H"
#include "4C_solid_3D_ele_calc_lib_io.hpp"
#include "4C_solid_3D_ele_calc_lib_mulf.H"
#include "4C_solid_3D_ele_formulation.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{
  template <CORE::FE::CellType celltype>
  struct MulfFBarPreparationData
  {
    CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_nodes<celltype>> N_XYZ{};

    SpatialMaterialMapping<celltype> spatial_material_mapping{};
  };

  namespace DETAILS
  {
    template <CORE::FE::CellType celltype>
    SpatialMaterialMapping<celltype> EvaluateMulfSpatialMaterialMappingCentroid(
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions_centeroid,
        const ElementNodes<celltype>& nodal_coordinates,
        const MulfHistoryData<celltype>& mulf_data_centeroid)
    {
      CORE::LINALG::Matrix<3, 3> delta_F = EvaluateMulfDeformationGradientUpdate(
          shape_functions_centeroid, nodal_coordinates.displacements_, mulf_data_centeroid);

      SpatialMaterialMapping<celltype> spatial_material_mapping_centeroid{};

      spatial_material_mapping_centeroid.deformation_gradient_.Multiply(
          delta_F, mulf_data_centeroid.deformation_gradient);

      spatial_material_mapping_centeroid.inverse_deformation_gradient_ =
          spatial_material_mapping_centeroid.deformation_gradient_;

      spatial_material_mapping_centeroid.determinant_deformation_gradient_ =
          spatial_material_mapping_centeroid.inverse_deformation_gradient_.Invert();

      return spatial_material_mapping_centeroid;
    }

    template <CORE::FE::CellType celltype>
    SpatialMaterialMapping<celltype> GetSpatialMaterialMappingBar(
        SpatialMaterialMapping<celltype> spatial_material_mapping, const double fbar_factor)
    {
      spatial_material_mapping.deformation_gradient_.Scale(fbar_factor);
      spatial_material_mapping.determinant_deformation_gradient_ *=
          CORE::FE::dim<celltype> * fbar_factor;
      spatial_material_mapping.inverse_deformation_gradient_.Scale(1.0 / fbar_factor);

      return spatial_material_mapping;
    }

    /*!
     * @brief Do a MULF update step on the mulf data
     */
    template <CORE::FE::CellType celltype>
    void UpdateMulfHistory(const ElementNodes<celltype>& element_nodes,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        MulfHistoryData<celltype>& mulf_data)
    {
      CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> delta_defgrd =
          EvaluateMulfDeformationGradientUpdate(
              shape_functions, element_nodes.displacements_, mulf_data);

      CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> inv_delta_defgrd(true);
      inv_delta_defgrd.Invert(delta_defgrd);

      CORE::LINALG::Matrix<3, 3> old_defgrd = mulf_data.deformation_gradient;
      mulf_data.deformation_gradient.Multiply(delta_defgrd, old_defgrd);

      CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> invJ_new;
      invJ_new.MultiplyTN(inv_delta_defgrd, mulf_data.inverse_jacobian);
      mulf_data.inverse_jacobian = std::move(invJ_new);
    }
  }  // namespace DETAILS



  /*!
   * @brief A solid element formulation with MULF prestressing and F-Bar
   *
   * @tparam celltype
   */
  template <CORE::FE::CellType celltype>
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

    static MulfFBarPreparationData<celltype> Prepare(const DRT::Element& ele,
        const ElementNodes<celltype>& nodal_coordinates, GlobalHistory& global_history)
    {
      if (!global_history.is_setup)
      {
        const JacobianMapping<celltype> jacobian_mapping =
            EvaluateJacobianMappingCentroid(nodal_coordinates);

        global_history.inverse_jacobian = jacobian_mapping.inverse_jacobian_;
        global_history.is_setup = true;
      }

      // set coordinates in parameter space at centroid as zero -> xi = [0; 0; 0]
      CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> xi_centroid =
          EvaluateParameterCoordinateCentroid<celltype>();

      // shape functions and derivatives evaluated at element centroid
      const ShapeFunctionsAndDerivatives<celltype> shape_functions_centeroid =
          EvaluateShapeFunctionsAndDerivs<celltype>(xi_centroid, nodal_coordinates);

      JacobianMapping<celltype> jacobian_mapping =
          EvaluateJacobianMappingCentroid(nodal_coordinates);

      CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::num_nodes<celltype>> N_XYZ_0;
      N_XYZ_0.Multiply(jacobian_mapping.inverse_jacobian_, shape_functions_centeroid.derivatives_);

      return {N_XYZ_0, DETAILS::EvaluateMulfSpatialMaterialMappingCentroid(
                           shape_functions_centeroid, nodal_coordinates, global_history)};
    }

    template <typename Evaluator>
    static auto Evaluate(const DRT::Element& ele, const ElementNodes<celltype>& element_nodes,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
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
          EvaluateMulfSpatialMaterialMapping(
              jacobian_mapping, shape_functions, element_nodes.displacements_, mulf_data_gp);

      const double fbar_factor = EvaluateFbarFactor(
          mapping_center.spatial_material_mapping.determinant_deformation_gradient_,
          spatial_material_mapping.determinant_deformation_gradient_);

      const FBarLinearizationContainer<celltype> linearization = std::invoke(
          [&]()
          {
            FBarLinearizationContainer<celltype> linearization{};
            linearization.Bop = EvaluateStrainGradient(jacobian_mapping, spatial_material_mapping);

            linearization.Hop = EvaluateFbarHOperator(jacobian_mapping.N_XYZ_, mapping_center.N_XYZ,
                spatial_material_mapping, mapping_center.spatial_material_mapping);

            linearization.fbar_factor = fbar_factor;

            linearization.cauchygreen = EvaluateCauchyGreen(spatial_material_mapping);

            return linearization;
          });

      const SpatialMaterialMapping<celltype> spatial_material_mapping_bar =
          DETAILS::GetSpatialMaterialMappingBar(spatial_material_mapping, fbar_factor);

      const CORE::LINALG::Matrix<CORE::FE::dim<celltype>, CORE::FE::dim<celltype>> cauchygreen_bar =
          EvaluateCauchyGreen<celltype>(spatial_material_mapping_bar);

      const CORE::LINALG::Matrix<DETAILS::num_str<celltype>, 1> gl_strain_bar =
          EvaluateGreenLagrangeStrain(cauchygreen_bar);

      return evaluator(
          spatial_material_mapping_bar.deformation_gradient_, gl_strain_bar, linearization);
    }


    static CORE::LINALG::Matrix<DETAILS::num_str<celltype>,
        CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>
    GetLinearBOperator(const FBarLinearizationContainer<celltype>& linearization)
    {
      return linearization.Bop;
    }

    static void AddInternalForceVector(const FBarLinearizationContainer<celltype>& linearization,
        const Stress<celltype>& stress, const double integration_factor,
        const MulfFBarPreparationData<celltype>& mapping_center,
        MulfHistoryData<celltype>& mulf_data_centeroid, MulfHistoryData<celltype>& mulf_data_gp,
        CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>, 1>&
            force_vector)
    {
      DRT::ELEMENTS::AddInternalForceVector(
          linearization.Bop, stress, integration_factor / linearization.fbar_factor, force_vector);
    }

    static void AddStiffnessMatrix(const FBarLinearizationContainer<celltype>& linearization,
        const JacobianMapping<celltype>& jacobian_mapping, const Stress<celltype>& stress,
        const double integration_factor, const MulfFBarPreparationData<celltype>& mapping_center,
        MulfHistoryData<celltype>& mulf_data_centeroid, MulfHistoryData<celltype>& mulf_data_gp,
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

    static void Pack(const MulfHistoryData<celltype>& history_data, CORE::COMM::PackBuffer& data)
    {
      CORE::COMM::ParObject::AddtoPack(data, history_data.inverse_jacobian);
      CORE::COMM::ParObject::AddtoPack(data, history_data.deformation_gradient);
      CORE::COMM::ParObject::AddtoPack(data, history_data.is_setup);
    }

    static void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data,
        MulfHistoryData<celltype>& history_data)
    {
      CORE::COMM::ParObject::ExtractfromPack(position, data, history_data.inverse_jacobian);
      CORE::COMM::ParObject::ExtractfromPack(position, data, history_data.deformation_gradient);
      int is_setup_int;
      CORE::COMM::ParObject::ExtractfromPack(position, data, is_setup_int);
      history_data.is_setup = static_cast<bool>(is_setup_int);
    }

    static inline void UpdatePrestress(const DRT::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const MulfFBarPreparationData<celltype>& mapping_center,
        MulfHistoryData<celltype>& mulf_data_centeroid)
    {
      CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1> xi_centroid =
          EvaluateParameterCoordinateCentroid<celltype>();

      ShapeFunctionsAndDerivatives<celltype> shape_functions_centeroid =
          EvaluateShapeFunctionsAndDerivs<celltype>(xi_centroid, element_nodes);

      DETAILS::UpdateMulfHistory(element_nodes, shape_functions_centeroid, mulf_data_centeroid);
    }

    static inline void UpdatePrestress(const DRT::Element& ele,
        const ElementNodes<celltype>& element_nodes,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
            deformation_gradient,
        const MulfFBarPreparationData<celltype>& mapping_center,
        MulfHistoryData<celltype>& mulf_data_centeroid, MulfHistoryData<celltype>& mulf_data_gp)
    {
      DETAILS::UpdateMulfHistory(element_nodes, shape_functions, mulf_data_gp);
    }
  };

  template <CORE::FE::CellType celltype>
  using MulfFBarSolidIntegrator = SolidEleCalc<celltype, MulfFBarFormulation<celltype>>;



}  // namespace DRT::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif
