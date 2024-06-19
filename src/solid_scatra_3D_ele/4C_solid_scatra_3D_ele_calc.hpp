/*! \file

\brief Definition of a solid-scatra element

\level 1
*/

#ifndef FOUR_C_SOLID_SCATRA_3D_ELE_CALC_HPP
#define FOUR_C_SOLID_SCATRA_3D_ELE_CALC_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_solid_3D_ele_calc_interface.hpp"
#include "4C_solid_3D_ele_formulation.hpp"
#include "4C_solid_scatra_3D_ele_calc_lib_nitsche.hpp"
#include "4C_structure_new_gauss_point_data_output_manager.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  class So3Material;
}  // namespace Mat
namespace STR::ELEMENTS
{
  class ParamsInterface;
}  // namespace STR::ELEMENTS



namespace Discret::ELEMENTS
{
  template <Core::FE::CellType celltype, typename SolidFormulation>
  class SolidScatraEleCalc
  {
   public:
    SolidScatraEleCalc();

    void pack(Core::Communication::PackBuffer& data) const;

    void unpack(std::vector<char>::size_type& position, const std::vector<char>& data);

    void setup(Mat::So3Material& solid_material, Input::LineDefinition* linedef);

    void material_post_setup(const Core::Elements::Element& ele, Mat::So3Material& solid_material);

    void Recover(const Core::Elements::Element& ele, const Core::FE::Discretization& discretization,
        const Core::Elements::Element::LocationArray& la, Teuchos::ParameterList& params);

    void Update(const Core::Elements::Element& ele, Mat::So3Material& solid_material,
        const Core::FE::Discretization& discretization,
        const Core::Elements::Element::LocationArray& la, Teuchos::ParameterList& params);

    void calculate_stress(const Core::Elements::Element& ele, Mat::So3Material& solid_material,
        const StressIO& stressIO, const StrainIO& strainIO,
        const Core::FE::Discretization& discretization,
        const Core::Elements::Element::LocationArray& la, Teuchos::ParameterList& params);

    double calculate_internal_energy(const Core::Elements::Element& ele,
        Mat::So3Material& solid_material, const Core::FE::Discretization& discretization,
        const Core::Elements::Element::LocationArray& la, Teuchos::ParameterList& params);

    void initialize_gauss_point_data_output(const Core::Elements::Element& ele,
        const Mat::So3Material& solid_material,
        STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const;

    void evaluate_gauss_point_data_output(const Core::Elements::Element& ele,
        const Mat::So3Material& solid_material,
        STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const;

    void reset_to_last_converged(
        const Core::Elements::Element& ele, Mat::So3Material& solid_material);

    void evaluate_nonlinear_force_stiffness_mass(const Core::Elements::Element& ele,
        Mat::So3Material& solid_material, const Core::FE::Discretization& discretization,
        const Core::Elements::Element::LocationArray& la, Teuchos::ParameterList& params,
        Core::LinAlg::SerialDenseVector* force_vector,
        Core::LinAlg::SerialDenseMatrix* stiffness_matrix,
        Core::LinAlg::SerialDenseMatrix* mass_matrix);

    void evaluate_d_stress_d_scalar(const Core::Elements::Element& ele,
        Mat::So3Material& solid_material, const Core::FE::Discretization& discretization,
        const Core::Elements::Element::LocationArray& la, Teuchos::ParameterList& params,
        Core::LinAlg::SerialDenseMatrix& stiffness_matrix_dScalar);

    double get_normal_cauchy_stress_at_xi(const Core::Elements::Element& ele,
        Mat::So3Material& solid_material, const std::vector<double>& disp,
        const std::optional<std::vector<double>>& scalars, const Core::LinAlg::Matrix<3, 1>& xi,
        const Core::LinAlg::Matrix<3, 1>& n, const Core::LinAlg::Matrix<3, 1>& dir,
        SolidScatraCauchyNDirLinearizations<3>& linearizations);

   private:
    /// static values for matrix sizes
    static constexpr int num_nodes_ = Core::FE::num_nodes<celltype>;
    static constexpr int num_dim_ = Core::FE::dim<celltype>;
    static constexpr int num_dof_per_ele_ = num_nodes_ * num_dim_;
    static constexpr int num_str_ = num_dim_ * (num_dim_ + 1) / 2;

    Core::FE::GaussIntegration stiffness_matrix_integration_;
    Core::FE::GaussIntegration mass_matrix_integration_;

    SolidFormulationHistory<SolidFormulation> history_data_{};
  };
}  // namespace Discret::ELEMENTS


FOUR_C_NAMESPACE_CLOSE

#endif