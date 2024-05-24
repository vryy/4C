/*! \file

\brief This file contains helper functions and type traits for the solid formulation

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_LIB_FORMULATION_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_LIB_FORMULATION_HPP

#include "4C_config.hpp"

#include "4C_lib_element.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_solid_3D_ele_formulation.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{

  template <CORE::FE::CellType celltype, typename SolidFormulation>
  class ElementFormulationDerivativeEvaluator
  {
   public:
    ElementFormulationDerivativeEvaluator(const DRT::Element& ele,
        const DRT::ELEMENTS::ElementNodes<celltype>& element_nodes,
        const CORE::LINALG::Matrix<3, 1>& xi,
        const DRT::ELEMENTS::ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const DRT::ELEMENTS::JacobianMapping<celltype>& jacobian_mapping,
        const CORE::LINALG::Matrix<3, 3>& deformation_gradient,
        const DRT::ELEMENTS::PreparationData<SolidFormulation>& preparation_data,
        DRT::ELEMENTS::SolidFormulationHistory<SolidFormulation>& history_data)
        : ele_(ele),
          element_nodes_(element_nodes),
          xi_(xi),
          shape_functions_(shape_functions),
          jacobian_mapping_(jacobian_mapping),
          deformation_gradient_(deformation_gradient),
          preparation_data_(preparation_data),
          history_data_(history_data)
    {
    }

    [[nodiscard]] CORE::LINALG::Matrix<9, CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_displacements() const
    {
      return DRT::ELEMENTS::evaluate_d_deformation_gradient_d_displacements(ele_, element_nodes_,
          xi_, shape_functions_, jacobian_mapping_, deformation_gradient_, preparation_data_,
          history_data_);
    }

    [[nodiscard]] CORE::LINALG::Matrix<9, CORE::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_xi() const
    {
      return DRT::ELEMENTS::evaluate_d_deformation_gradient_d_xi(ele_, element_nodes_, xi_,
          shape_functions_, jacobian_mapping_, deformation_gradient_, preparation_data_,
          history_data_);
    }

    [[nodiscard]] CORE::LINALG::Matrix<9,
        CORE::FE::num_nodes<celltype> * CORE::FE::dim<celltype> * CORE::FE::dim<celltype>>
    evaluate_d2_deformation_gradient_d_displacements_d_xi() const
    {
      return DRT::ELEMENTS::evaluate_d_deformation_gradient_d_displacements_d_xi(ele_,
          element_nodes_, xi_, shape_functions_, jacobian_mapping_, deformation_gradient_,
          preparation_data_, history_data_);
    }

   private:
    const DRT::Element& ele_;
    const DRT::ELEMENTS::ElementNodes<celltype>& element_nodes_;
    const CORE::LINALG::Matrix<3, 1>& xi_;
    const DRT::ELEMENTS::ShapeFunctionsAndDerivatives<celltype>& shape_functions_;
    const DRT::ELEMENTS::JacobianMapping<celltype>& jacobian_mapping_;
    const CORE::LINALG::Matrix<3, 3>& deformation_gradient_;
    const DRT::ELEMENTS::PreparationData<SolidFormulation>& preparation_data_;
    DRT::ELEMENTS::SolidFormulationHistory<SolidFormulation>& history_data_;
  };
}  // namespace DRT::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif