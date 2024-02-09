/*! \file

\brief A displacement based solid element formulation

\level 1
*/

#ifndef BACI_SOLID_ELE_CALC_DISPLACEMENT_BASED_HPP
#define BACI_SOLID_ELE_CALC_DISPLACEMENT_BASED_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_cell_type_traits.hpp"
#include "baci_lib_element.hpp"
#include "baci_solid_ele_calc.hpp"
#include "baci_solid_ele_calc_lib.hpp"

BACI_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{
  struct DisplacementBasedPreparationData
  {
    // no preparation data needed
  };

  struct DisplacementBasedHistoryData
  {
    // no history data needed
  };

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
    static DisplacementBasedPreparationData Prepare(const DRT::Element& ele,
        const ElementNodes<celltype>& nodal_coordinates, DisplacementBasedHistoryData& history_data)
    {
      // do nothing for simple displacement based evaluation
      return {};
    }

    template <typename Evaluator>
    static void Evaluate(const DRT::Element& ele, const ElementNodes<celltype>& nodal_coordinates,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const DisplacementBasedPreparationData& preparation_data,
        DisplacementBasedHistoryData& history_data, Evaluator evaluator)
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

      evaluator(spatial_material_mapping.deformation_gradient_, gl_strain, linearization);
    }

    static void AddInternalForceVector(
        const DisplacementBasedLinearizationContainer<celltype>& linearization,
        const Stress<celltype>& stress, const double integration_factor,
        const DisplacementBasedPreparationData& preparation_data,
        DisplacementBasedHistoryData& history_data,
        CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>, 1>&
            force_vector)
    {
      DRT::ELEMENTS::AddInternalForceVector(
          linearization.Bop_, stress, integration_factor, force_vector);
    }

    static void AddStiffnessMatrix(
        const DisplacementBasedLinearizationContainer<celltype>& linearization,
        const JacobianMapping<celltype>& jacobian_mapping, const Stress<celltype>& stress,
        const double integration_factor, const DisplacementBasedPreparationData& preparation_data,
        DisplacementBasedHistoryData& history_data,
        CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>,
            CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>& stiffness_matrix)
    {
      DRT::ELEMENTS::AddElasticStiffnessMatrix(
          linearization.Bop_, stress, integration_factor, stiffness_matrix);
      DRT::ELEMENTS::AddGeometricStiffnessMatrix(
          jacobian_mapping.N_XYZ_, stress, integration_factor, stiffness_matrix);
    }

    static void Pack(const DisplacementBasedHistoryData& history_data, CORE::COMM::PackBuffer& data)
    {
      // nothing to pack
    }

    static void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data,
        DisplacementBasedHistoryData& history_data)
    {
      // nothing to unpack
    }
  };

  template <CORE::FE::CellType celltype>
  using DisplacementBasedSolidIntegrator =
      SolidEleCalc<celltype, DisplacementBasedFormulation<celltype>,
          DisplacementBasedPreparationData, DisplacementBasedHistoryData>;


}  // namespace DRT::ELEMENTS

BACI_NAMESPACE_CLOSE
#endif  // BACI_SOLID_ELE_CALC_DISPLACEMENT_BASED_H
