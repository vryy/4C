/*! \file

\brief Definition of a solid-scatra element

\level 1
*/

#ifndef FOUR_C_SOLID_SCATRA_3D_ELE_CALC_HPP
#define FOUR_C_SOLID_SCATRA_3D_ELE_CALC_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_utils_gausspoints.hpp"
#include "4C_lib_element.hpp"
#include "4C_solid_3D_ele_calc_interface.hpp"
#include "4C_solid_3D_ele_formulation.hpp"
#include "4C_solid_scatra_3D_ele_calc_lib_nitsche.hpp"
#include "4C_structure_new_gauss_point_data_output_manager.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  class So3Material;
}  // namespace MAT
namespace STR::ELEMENTS
{
  class ParamsInterface;
}  // namespace STR::ELEMENTS



namespace DRT::ELEMENTS
{
  template <CORE::FE::CellType celltype, typename SolidFormulation>
  class SolidScatraEleCalc
  {
   public:
    SolidScatraEleCalc();

    void Pack(CORE::COMM::PackBuffer& data) const;

    void Unpack(std::vector<char>::size_type& position, const std::vector<char>& data);

    void Setup(MAT::So3Material& solid_material, INPUT::LineDefinition* linedef);

    void material_post_setup(const DRT::Element& ele, MAT::So3Material& solid_material);

    void Recover(const DRT::Element& ele, const DRT::Discretization& discretization,
        const DRT::Element::LocationArray& la, Teuchos::ParameterList& params);

    void Update(const DRT::Element& ele, MAT::So3Material& solid_material,
        const DRT::Discretization& discretization, const DRT::Element::LocationArray& la,
        Teuchos::ParameterList& params);

    void CalculateStress(const DRT::Element& ele, MAT::So3Material& solid_material,
        const StressIO& stressIO, const StrainIO& strainIO,
        const DRT::Discretization& discretization, const DRT::Element::LocationArray& la,
        Teuchos::ParameterList& params);

    double calculate_internal_energy(const DRT::Element& ele, MAT::So3Material& solid_material,
        const DRT::Discretization& discretization, const DRT::Element::LocationArray& la,
        Teuchos::ParameterList& params);

    void initialize_gauss_point_data_output(const DRT::Element& ele,
        const MAT::So3Material& solid_material,
        STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const;

    void evaluate_gauss_point_data_output(const DRT::Element& ele,
        const MAT::So3Material& solid_material,
        STR::MODELEVALUATOR::GaussPointDataOutputManager& gp_data_output_manager) const;

    void reset_to_last_converged(const DRT::Element& ele, MAT::So3Material& solid_material);

    void evaluate_nonlinear_force_stiffness_mass(const DRT::Element& ele,
        MAT::So3Material& solid_material, const DRT::Discretization& discretization,
        const DRT::Element::LocationArray& la, Teuchos::ParameterList& params,
        CORE::LINALG::SerialDenseVector* force_vector,
        CORE::LINALG::SerialDenseMatrix* stiffness_matrix,
        CORE::LINALG::SerialDenseMatrix* mass_matrix);

    void evaluate_d_stress_d_scalar(const DRT::Element& ele, MAT::So3Material& solid_material,
        const DRT::Discretization& discretization, const DRT::Element::LocationArray& la,
        Teuchos::ParameterList& params, CORE::LINALG::SerialDenseMatrix& stiffness_matrix_dScalar);

    double GetCauchyNDirAtXi(const DRT::Element& ele, MAT::So3Material& solid_material,
        const std::vector<double>& disp, const std::optional<std::vector<double>>& scalars,
        const CORE::LINALG::Matrix<3, 1>& xi, const CORE::LINALG::Matrix<3, 1>& n,
        const CORE::LINALG::Matrix<3, 1>& dir,
        SolidScatraCauchyNDirLinearizations<3>& linearizations);

   private:
    /// static values for matrix sizes
    static constexpr int num_nodes_ = CORE::FE::num_nodes<celltype>;
    static constexpr int num_dim_ = CORE::FE::dim<celltype>;
    static constexpr int num_dof_per_ele_ = num_nodes_ * num_dim_;
    static constexpr int num_str_ = num_dim_ * (num_dim_ + 1) / 2;

    CORE::FE::GaussIntegration stiffness_matrix_integration_;
    CORE::FE::GaussIntegration mass_matrix_integration_;

    SolidFormulationHistory<SolidFormulation> history_data_{};
  };
}  // namespace DRT::ELEMENTS


FOUR_C_NAMESPACE_CLOSE

#endif