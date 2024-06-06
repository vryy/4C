/*! \file

\brief This file contains helper functions and type traits for the solid formulation

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_LIB_FORMULATION_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_LIB_FORMULATION_HPP

#include "4C_config.hpp"

#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_solid_3D_ele_formulation.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret::ELEMENTS
{

  template <Core::FE::CellType celltype, typename SolidFormulation>
  class ElementFormulationDerivativeEvaluator
  {
   public:
    ElementFormulationDerivativeEvaluator(const Core::Elements::Element& ele,
        const Discret::ELEMENTS::ElementNodes<celltype>& element_nodes,
        const Core::LinAlg::Matrix<3, 1>& xi,
        const Discret::ELEMENTS::ShapeFunctionsAndDerivatives<celltype>& shape_functions,
        const Discret::ELEMENTS::JacobianMapping<celltype>& jacobian_mapping,
        const Core::LinAlg::Matrix<3, 3>& deformation_gradient,
        const Discret::ELEMENTS::PreparationData<SolidFormulation>& preparation_data,
        Discret::ELEMENTS::SolidFormulationHistory<SolidFormulation>& history_data)
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

    [[nodiscard]] Core::LinAlg::Matrix<9, Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_displacements() const
    {
      return Discret::ELEMENTS::evaluate_d_deformation_gradient_d_displacements(ele_,
          element_nodes_, xi_, shape_functions_, jacobian_mapping_, deformation_gradient_,
          preparation_data_, history_data_);
    }

    [[nodiscard]] Core::LinAlg::Matrix<9, Core::FE::dim<celltype>>
    evaluate_d_deformation_gradient_d_xi() const
    {
      return Discret::ELEMENTS::evaluate_d_deformation_gradient_d_xi(ele_, element_nodes_, xi_,
          shape_functions_, jacobian_mapping_, deformation_gradient_, preparation_data_,
          history_data_);
    }

    [[nodiscard]] Core::LinAlg::Matrix<9,
        Core::FE::num_nodes<celltype> * Core::FE::dim<celltype> * Core::FE::dim<celltype>>
    evaluate_d2_deformation_gradient_d_displacements_d_xi() const
    {
      return Discret::ELEMENTS::evaluate_d_deformation_gradient_d_displacements_d_xi(ele_,
          element_nodes_, xi_, shape_functions_, jacobian_mapping_, deformation_gradient_,
          preparation_data_, history_data_);
    }

   private:
    const Core::Elements::Element& ele_;
    const Discret::ELEMENTS::ElementNodes<celltype>& element_nodes_;
    const Core::LinAlg::Matrix<3, 1>& xi_;
    const Discret::ELEMENTS::ShapeFunctionsAndDerivatives<celltype>& shape_functions_;
    const Discret::ELEMENTS::JacobianMapping<celltype>& jacobian_mapping_;
    const Core::LinAlg::Matrix<3, 3>& deformation_gradient_;
    const Discret::ELEMENTS::PreparationData<SolidFormulation>& preparation_data_;
    Discret::ELEMENTS::SolidFormulationHistory<SolidFormulation>& history_data_;
  };
}  // namespace Discret::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif