// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_3D_ELE_CALC_EAS_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_EAS_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_solid_3D_ele_calc_eas_helpers.hpp"
#include "4C_solid_3D_ele_calc_interface.hpp"
#include "4C_solid_3D_ele_interface_serializable.hpp"
#include "4C_solid_3D_ele_utils.hpp"

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  class So3Material;
}

namespace Solid::ModelEvaluator
{
  class GaussPointDataOutputManager;
}
namespace Discret
{

  namespace Elements
  {
    template <Core::FE::CellType celltype, Discret::Elements::EasType eastype,
        Inpar::Solid::KinemType>
    class SolidEleCalcEas
    {
     public:
      SolidEleCalcEas();

      void setup(
          Mat::So3Material& solid_material, const Core::IO::InputParameterContainer& container);

      void pack(Core::Communication::PackBuffer& data) const;

      void unpack(Core::Communication::UnpackBuffer& buffer);

      void material_post_setup(
          const Core::Elements::Element& ele, Mat::So3Material& solid_material);

      void evaluate_nonlinear_force_stiffness_mass(const Core::Elements::Element& ele,
          Mat::So3Material& solid_material, const Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params,
          Core::LinAlg::SerialDenseVector* force_vector,
          Core::LinAlg::SerialDenseMatrix* stiffness_matrix,
          Core::LinAlg::SerialDenseMatrix* mass_matrix);

      void recover(Core::Elements::Element& ele, const Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params);

      void calculate_stress(const Core::Elements::Element& ele, Mat::So3Material& solid_material,
          const StressIO& stressIO, const StrainIO& strainIO,
          const Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params);

      double calculate_internal_energy(const Core::Elements::Element& ele,
          Mat::So3Material& solid_material, const Core::FE::Discretization& discretization,
          const std::vector<int>& lm, Teuchos::ParameterList& params);

      void update(const Core::Elements::Element& ele, Mat::So3Material& solid_material,
          const Core::FE::Discretization& discretization, const std::vector<int>& lm,
          Teuchos::ParameterList& params);

      void initialize_gauss_point_data_output(const Core::Elements::Element& ele,
          const Mat::So3Material& solid_material,
          Solid::ModelEvaluator::GaussPointDataOutputManager& gp_data_output_manager) const;

      void evaluate_gauss_point_data_output(const Core::Elements::Element& ele,
          const Mat::So3Material& solid_material,
          Solid::ModelEvaluator::GaussPointDataOutputManager& gp_data_output_manager) const;

      void reset_to_last_converged(
          const Core::Elements::Element& ele, Mat::So3Material& solid_material);

      void for_each_gauss_point(const Core::Elements::Element& ele,
          Mat::So3Material& solid_material, const Core::FE::Discretization& discretization,
          const std::vector<int>& lm,
          const std::function<void(Mat::So3Material&, double integration_factor, int gp)>&
              integrator) const;

      void set_integration_rule(const Core::FE::GaussIntegration& integration_rule)
      {
        stiffness_matrix_integration_ = integration_rule;
        mass_matrix_integration_ = integration_rule;
      }

      const Core::FE::GaussIntegration& get_gauss_rule_stiffness_integration() const
      {
        return stiffness_matrix_integration_;
      };

     private:
      /// EAS matrices and vectors to be stored between iterations
      Discret::Elements::EasIterationData<celltype, eastype> eas_iteration_data_ = {};

      // line search parameter (old step length)
      double old_step_length_;

      /// static values for matrix sizes
      static constexpr int num_nodes_ = Core::FE::num_nodes<celltype>;
      static constexpr int num_dim_ = Core::FE::dim<celltype>;
      static constexpr int num_dof_per_ele_ = num_nodes_ * num_dim_;
      static constexpr int num_str_ = num_dim_ * (num_dim_ + 1) / 2;

      Core::FE::GaussIntegration stiffness_matrix_integration_;
      Core::FE::GaussIntegration mass_matrix_integration_;
    };  // class SolidEleCalcEas
  }     // namespace Elements
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
