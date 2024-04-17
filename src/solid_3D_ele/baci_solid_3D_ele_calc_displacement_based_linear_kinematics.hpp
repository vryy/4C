/*! \file

\brief A displacement based solid element formulation with linear kinematics

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_DISPLACEMENT_BASED_LINEAR_KINEMATICS_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_DISPLACEMENT_BASED_LINEAR_KINEMATICS_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_cell_type.hpp"
#include "baci_discretization_fem_general_cell_type_traits.hpp"
#include "baci_lib_element.hpp"
#include "baci_solid_3D_ele_calc.hpp"
#include "baci_solid_3D_ele_calc_lib.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{
  namespace DETAILS
  {
    template <CORE::FE::CellType celltype,
        std::enable_if_t<DETAIL::num_dim<celltype> == 3, int> = 0>
    CORE::LINALG::Matrix<DETAIL::num_str<celltype>,
        DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>
    EvaluateLinearStrainGradient(const JacobianMapping<celltype>& jacobian_mapping)
    {
      // B-operator
      CORE::LINALG::Matrix<DETAIL::num_str<celltype>,
          DETAIL::num_dim<celltype> * DETAIL::num_nodes<celltype>>
          Bop;
      for (int i = 0; i < DETAIL::num_nodes<celltype>; ++i)
      {
        for (int d = 0; d < DETAIL::num_dim<celltype>; ++d)
        {
          Bop(d, DETAIL::num_dim<celltype> * i + d) = jacobian_mapping.N_XYZ_(d, i);
        }

        Bop(3, DETAIL::num_dim<celltype> * i + 0) = jacobian_mapping.N_XYZ_(1, i);
        Bop(3, DETAIL::num_dim<celltype> * i + 1) = jacobian_mapping.N_XYZ_(0, i);
        Bop(3, DETAIL::num_dim<celltype> * i + 2) = 0;
        Bop(4, DETAIL::num_dim<celltype> * i + 0) = 0;
        Bop(4, DETAIL::num_dim<celltype> * i + 1) = jacobian_mapping.N_XYZ_(2, i);
        Bop(4, DETAIL::num_dim<celltype> * i + 2) = jacobian_mapping.N_XYZ_(1, i);
        Bop(5, DETAIL::num_dim<celltype> * i + 0) = jacobian_mapping.N_XYZ_(2, i);
        Bop(5, DETAIL::num_dim<celltype> * i + 1) = 0;
        Bop(5, DETAIL::num_dim<celltype> * i + 2) = jacobian_mapping.N_XYZ_(0, i);
      }

      return Bop;
    }
  }  // namespace DETAILS

  /*!
   * @brief A container holding the linearization of the displacement based linear kinematics solid
   * element formulation
   */
  template <CORE::FE::CellType celltype>
  struct DisplacementBasedLinearKinematicsLinearizationContainer
  {
    CORE::LINALG::Matrix<DETAILS::num_str<celltype>,
        CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>
        linear_b_operator_{};
  };

  /*!
   * @brief A displacement based solid element formulation with linear kinematics (i.e. small
   * displacements)
   *
   * @tparam celltype
   */
  template <CORE::FE::CellType celltype>
  struct DisplacementBasedLinearKinematicsFormulation
  {
    static constexpr bool has_gauss_point_history = false;
    static constexpr bool has_global_history = false;
    static constexpr bool has_preparation_data = false;


    using LinearizationContainer =
        DisplacementBasedLinearKinematicsLinearizationContainer<celltype>;

    template <typename Evaluator>
    static inline auto Evaluate(const DRT::Element& ele,
        const ElementNodes<celltype>& nodal_coordinates,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping, Evaluator evaluator)
    {
      const DisplacementBasedLinearKinematicsLinearizationContainer<celltype> linearization{
          DETAILS::EvaluateLinearStrainGradient(jacobian_mapping)};

      CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>, 1> nodal_displs(
          true);

      for (unsigned i = 0; i < CORE::FE::num_nodes<celltype>; ++i)
        for (unsigned j = 0; j < CORE::FE::dim<celltype>; ++j)
          nodal_displs(i * CORE::FE::dim<celltype> + j, 0) = nodal_coordinates.displacements_(i, j);

      CORE::LINALG::Matrix<DETAIL::num_str<celltype>, 1> gl_strain;
      gl_strain.Multiply(linearization.linear_b_operator_, nodal_displs);

      return evaluator(
          CORE::LINALG::IdentityMatrix<CORE::FE::dim<celltype>>(), gl_strain, linearization);
    }

    static inline SolidFormulationLinearization<celltype> EvaluateFullLinearization(
        const DRT::Element& ele, const ElementNodes<celltype>& nodal_coordinates,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
        const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const JacobianMapping<celltype>& jacobian_mapping,
        const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, DETAIL::num_dim<celltype>>&
            deformation_gradient)
    {
      dserror(
          "The full linearization is not yet implemented for the displacement based formulation "
          "with linear kinematics.");
    }

    static void AddInternalForceVector(
        const DisplacementBasedLinearKinematicsLinearizationContainer<celltype>& linearization,
        const Stress<celltype>& stress, const double integration_factor,
        CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>, 1>&
            force_vector)
    {
      DRT::ELEMENTS::AddInternalForceVector(
          linearization.linear_b_operator_, stress, integration_factor, force_vector);
    }

    static void AddStiffnessMatrix(
        const DisplacementBasedLinearKinematicsLinearizationContainer<celltype>& linearization,
        const JacobianMapping<celltype>& jacobian_mapping, const Stress<celltype>& stress,
        const double integration_factor,
        CORE::LINALG::Matrix<CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>,
            CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>& stiffness_matrix)
    {
      DRT::ELEMENTS::AddElasticStiffnessMatrix(
          linearization.linear_b_operator_, stress, integration_factor, stiffness_matrix);
    }
  };

  template <CORE::FE::CellType celltype>
  using DisplacementBasedLinearKinematicsSolidIntegrator =
      SolidEleCalc<celltype, DisplacementBasedLinearKinematicsFormulation<celltype>>;


}  // namespace DRT::ELEMENTS

FOUR_C_NAMESPACE_CLOSE
#endif