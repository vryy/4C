/*! \file

\brief Declaration of routines for calculation of solid element with EAS element technology

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_EAS_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_EAS_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_solid_3D_ele_calc_interface.hpp"
#include "4C_solid_3D_ele_interface_serializable.hpp"
#include "4C_solid_3D_ele_utils.hpp"

#include <memory>
#include <string>
#include <unordered_map>

FOUR_C_NAMESPACE_OPEN

namespace Solid::ELEMENTS
{
  enum class EasType
  {
    soh8_easnone,
    eastype_h8_9,
    eastype_h8_21,
    eastype_sh8_7,
    eastype_sh18_9,
    eastype_undefined
  };

  template <Solid::ELEMENTS::EasType eastype>
  struct EasTypeToNumEas
  {
  };
  template <>
  struct EasTypeToNumEas<Solid::ELEMENTS::EasType::eastype_h8_9>
  {
    static constexpr int num_eas = 9;
  };
  template <>
  struct EasTypeToNumEas<Solid::ELEMENTS::EasType::eastype_h8_21>
  {
    static constexpr int num_eas = 21;
  };
  template <>
  struct EasTypeToNumEas<Solid::ELEMENTS::EasType::eastype_sh8_7>
  {
    static constexpr int num_eas = 7;
  };
  template <>
  struct EasTypeToNumEas<Solid::ELEMENTS::EasType::eastype_sh18_9>
  {
    static constexpr int num_eas = 9;
  };
  template <>
  struct EasTypeToNumEas<Solid::ELEMENTS::EasType::eastype_undefined>
  {
  };

}  // namespace Solid::ELEMENTS

namespace Mat
{
  class So3Material;
}

namespace Solid::MODELEVALUATOR
{
  class GaussPointDataOutputManager;
}
namespace Discret
{

  namespace ELEMENTS
  {
    /// struct for EAS matrices and vectors to be stored between iterations
    template <Core::FE::CellType celltype, Solid::ELEMENTS::EasType eastype>
    struct EasIterationData
    {
      constexpr static int num_eas = Solid::ELEMENTS::EasTypeToNumEas<eastype>::num_eas;

      /// inverse EAS matrix K_{alpha alpha}
      Core::LinAlg::Matrix<num_eas, num_eas> invKaa_{true};

      /// EAS matrix K_{d alpha}
      Core::LinAlg::Matrix<Core::FE::num_nodes<celltype> * Core::FE::dim<celltype>, num_eas> Kda_{
          true};

      /// EAS enhacement vector s
      Core::LinAlg::Matrix<num_eas, 1> s_{true};

      /// discrete enhanced strain scalars increment
      Core::LinAlg::Matrix<num_eas, 1> alpha_inc_{true};

      /// discrete enhanced strain scalars alpha
      Core::LinAlg::Matrix<num_eas, 1> alpha_{true};
    };

    template <Core::FE::CellType celltype, Solid::ELEMENTS::EasType eastype,
        Inpar::Solid::KinemType>
    class SolidEleCalcEas
    {
     public:
      SolidEleCalcEas();

      void setup(Mat::So3Material& solid_material, Input::LineDefinition* linedef);

      void pack(Core::Communication::PackBuffer& data) const;

      void unpack(std::vector<char>::size_type& position, const std::vector<char>& data);

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
          Solid::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const;

      void evaluate_gauss_point_data_output(const Core::Elements::Element& ele,
          const Mat::So3Material& solid_material,
          Solid::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const;

      void reset_to_last_converged(
          const Core::Elements::Element& ele, Mat::So3Material& solid_material);

      void for_each_gauss_point(const Core::Elements::Element& ele,
          Mat::So3Material& solid_material, const Core::FE::Discretization& discretization,
          const std::vector<int>& lm,
          const std::function<void(Mat::So3Material&, double integration_factor, int gp)>&
              integrator) const;

     private:
      /// EAS matrices and vectors to be stored between iterations
      Discret::ELEMENTS::EasIterationData<celltype, eastype> eas_iteration_data_ = {};

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
  }     // namespace ELEMENTS
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
