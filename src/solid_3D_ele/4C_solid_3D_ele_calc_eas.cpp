// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_solid_3D_ele_calc_eas.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_fixedsizematrix_generators.hpp"
#include "4C_linalg_fixedsizematrix_solver.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_solid_3D_ele.hpp"
#include "4C_solid_3D_ele_calc_eas_helpers.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_solid_3D_ele_calc_lib_integration.hpp"
#include "4C_solid_3D_ele_calc_lib_io.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"
#include "4C_structure_new_gauss_point_data_output_manager.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_dyn_cast.hpp>
#include <Teuchos_ParameterList.hpp>

#include <memory>
#include <optional>

FOUR_C_NAMESPACE_OPEN

namespace
{
  /*!
   * @brief Solve for the inverse of a matrix and ignore any errors
   *
   * @tparam dim : matrix dimensions
   * @param matrix(in/out) : matrix to be inverted
   */
  template <unsigned int dim>
  void solve_for_inverse_ignoring_errors(Core::LinAlg::Matrix<dim, dim>& matrix)
  {
    Core::LinAlg::FixedSizeSerialDenseSolver<dim, dim, 1> solve_for_inverse;
    solve_for_inverse.set_matrix(matrix);

    solve_for_inverse.invert();
  }
}  // namespace

template <Core::FE::CellType celltype, Discret::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::SolidEleCalcEas()
    : stiffness_matrix_integration_(
          create_gauss_integration<celltype>(get_gauss_rule_stiffness_matrix<celltype>())),
      mass_matrix_integration_(
          create_gauss_integration<celltype>(get_gauss_rule_mass_matrix<celltype>()))
{
}

template <Core::FE::CellType celltype, Discret::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::pack(
    Core::Communication::PackBuffer& data) const
{
  add_to_pack(data, eas_iteration_data_.alpha_inc);
  add_to_pack(data, eas_iteration_data_.alpha);
  add_to_pack(data, eas_iteration_data_.s);
  add_to_pack(data, eas_iteration_data_.invKaa);
  add_to_pack(data, eas_iteration_data_.Kda);
  add_to_pack(data, old_step_length_);
};

template <Core::FE::CellType celltype, Discret::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::unpack(
    Core::Communication::UnpackBuffer& buffer)
{
  extract_from_pack(buffer, eas_iteration_data_.alpha_inc);
  extract_from_pack(buffer, eas_iteration_data_.alpha);
  extract_from_pack(buffer, eas_iteration_data_.s);
  extract_from_pack(buffer, eas_iteration_data_.invKaa);
  extract_from_pack(buffer, eas_iteration_data_.Kda);
  extract_from_pack(buffer, old_step_length_);
};

template <Core::FE::CellType celltype, Discret::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype,
    kinematic_type>::evaluate_nonlinear_force_stiffness_mass(const Core::Elements::Element& ele,
    Mat::So3Material& solid_material, const Core::FE::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params,
    Core::LinAlg::SerialDenseVector* force_vector,
    Core::LinAlg::SerialDenseMatrix* stiffness_matrix, Core::LinAlg::SerialDenseMatrix* mass_matrix)
{
  // Create views to SerialDenseMatrices
  std::optional<Core::LinAlg::Matrix<num_dof_per_ele_, num_dof_per_ele_>> stiff = {};
  std::optional<Core::LinAlg::Matrix<num_dof_per_ele_, num_dof_per_ele_>> mass = {};
  std::optional<Core::LinAlg::Matrix<num_dof_per_ele_, 1>> force = {};
  if (stiffness_matrix != nullptr) stiff.emplace(*stiffness_matrix, true);
  if (mass_matrix != nullptr) mass.emplace(*mass_matrix, true);
  if (force_vector != nullptr) force.emplace(*force_vector, true);

  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, lm);

  bool equal_integration_mass_stiffness =
      compare_gauss_integration(mass_matrix_integration_, stiffness_matrix_integration_);

  CentroidTransformation<celltype> centroid_transformation =
      evaluate_centroid_transformation<celltype>(nodal_coordinates);

  if (!ele.is_params_interface())
  {
    // Update alpha only in old time integration scheme
    update_alpha<celltype, eastype>(eas_iteration_data_, discretization, lm);
  }

  // clear for integration
  eas_iteration_data_.invKaa.clear();
  eas_iteration_data_.Kda.clear();
  eas_iteration_data_.s.clear();

  evaluate_centroid_coordinates_and_add_to_parameter_list<celltype>(nodal_coordinates, params);

  double element_mass = 0.0;
  double element_volume = 0.0;
  Discret::Elements::for_each_gauss_point<celltype>(nodal_coordinates,
      stiffness_matrix_integration_,
      [&](const Core::LinAlg::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        const EASKinematics<celltype, eastype> kinematic_quantitites =
            evaluate_eas_kinematics<celltype, eastype, kinematic_type>(nodal_coordinates,
                centroid_transformation, xi, jacobian_mapping, eas_iteration_data_);

        evaluate_gp_coordinates_and_add_to_parameter_list<celltype>(
            nodal_coordinates, shape_functions, params);

        const Stress<celltype> stress = evaluate_material_stress<celltype>(solid_material,
            kinematic_quantitites.enhanced_deformation_gradient, kinematic_quantitites.enhanced_gl,
            params, gp, ele.id());

        integrate_eas<celltype, eastype>(stress, kinematic_quantitites.m_tilde,
            kinematic_quantitites.b_op, integration_factor, eas_iteration_data_);

        if (force.has_value())
        {
          add_internal_force_vector(kinematic_quantitites.b_op, stress, integration_factor, *force);
        }

        if (stiff.has_value())
        {
          add_elastic_stiffness_matrix(
              kinematic_quantitites.b_op, stress, integration_factor, *stiff);
          add_geometric_stiffness_matrix(
              jacobian_mapping.N_XYZ_, stress, integration_factor, *stiff);
        }

        if (mass.has_value())
        {
          if (equal_integration_mass_stiffness)
          {
            add_mass_matrix(shape_functions, integration_factor, solid_material.density(gp), *mass);
          }
          else
          {
            element_mass += solid_material.density(gp) * integration_factor;
            element_volume += integration_factor;
          }
        }
      });

  // invert Kaa with solver. eas_iteration_data_.invKaa then is Kaa^{-1}
  solve_for_inverse_ignoring_errors(eas_iteration_data_.invKaa);

  // compute the product (- Kda Kaa^{-1}) which is later needed for force and stiffness update
  Core::LinAlg::Matrix<num_dof_per_ele_, Discret::Elements::EasTypeToNumEas<eastype>::num_eas>
      minusKdainvKaa(true);
  minusKdainvKaa.multiply_nn(-1.0, eas_iteration_data_.Kda, eas_iteration_data_.invKaa);

  if (force.has_value())
  {
    add_eas_internal_force<celltype, eastype>(minusKdainvKaa, eas_iteration_data_.s, *force);
  }

  if (stiff.has_value())
  {
    add_eas_stiffness_matrix<celltype, eastype>(minusKdainvKaa, eas_iteration_data_.Kda, *stiff);
  }

  if (mass.has_value() && !equal_integration_mass_stiffness)
  {
    // integrate mass matrix
    FOUR_C_ASSERT(element_mass > 0, "It looks like the element mass is 0.0");
    Discret::Elements::for_each_gauss_point<celltype>(nodal_coordinates, mass_matrix_integration_,
        [&](const Core::LinAlg::Matrix<num_dim_, 1>& xi,
            const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
            const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp) {
          add_mass_matrix(
              shape_functions, integration_factor, element_mass / element_volume, *mass);
        });
  }
}

template <Core::FE::CellType celltype, Discret::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::recover(
    Core::Elements::Element& ele, const Core::FE::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
  FourC::Solid::Elements::ParamsInterface& params_interface =
      *std::dynamic_pointer_cast<FourC::Solid::Elements::ParamsInterface>(
          ele.params_interface_ptr());

  const double step_length = params_interface.get_step_length();

  if (params_interface.is_default_step())
  {
    params_interface.sum_into_my_previous_sol_norm(NOX::Nln::StatusTest::quantity_eas,
        Discret::Elements::EasTypeToNumEas<eastype>::num_eas, &eas_iteration_data_.alpha(0, 0),
        ele.owner());

    // Update alpha
    update_alpha(eas_iteration_data_, discretization, lm, step_length);
  }
  else
  {
    correct_alpha(eas_iteration_data_, step_length, old_step_length_);
  }

  // store old step length
  old_step_length_ = step_length;

  params_interface.sum_into_my_update_norm(NOX::Nln::StatusTest::quantity_eas,
      Discret::Elements::EasTypeToNumEas<eastype>::num_eas, &eas_iteration_data_.alpha_inc(0, 0),
      &eas_iteration_data_.alpha(0, 0), step_length, ele.owner());
}

template <Core::FE::CellType celltype, Discret::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::update(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material,
    const Core::FE::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, lm);
  CentroidTransformation<celltype> centroid_transformation =
      evaluate_centroid_transformation<celltype>(nodal_coordinates);

  // No need to update alpha here. Update is called to copy states from t_{n+1} to t_{n} after the
  // time step and output. Hence, there are no more Newton iterations that would require an update
  // of alpha

  evaluate_centroid_coordinates_and_add_to_parameter_list<celltype>(nodal_coordinates, params);

  Discret::Elements::for_each_gauss_point<celltype>(nodal_coordinates,
      stiffness_matrix_integration_,
      [&](const Core::LinAlg::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        const EASKinematics<celltype, eastype> kinematic_quantitites =
            evaluate_eas_kinematics<celltype, eastype, kinematic_type>(nodal_coordinates,
                centroid_transformation, xi, jacobian_mapping, eas_iteration_data_);

        evaluate_gp_coordinates_and_add_to_parameter_list<celltype>(
            nodal_coordinates, shape_functions, params);

        solid_material.update(
            kinematic_quantitites.enhanced_deformation_gradient, gp, params, ele.id());
      });

  solid_material.update();
}


template <Core::FE::CellType celltype, Discret::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::calculate_stress(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material, const StressIO& stressIO,
    const StrainIO& strainIO, const Core::FE::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params)
{
  std::vector<char>& serialized_stress_data = stressIO.mutable_data;
  std::vector<char>& serialized_strain_data = strainIO.mutable_data;
  Core::LinAlg::SerialDenseMatrix stress_data(stiffness_matrix_integration_.num_points(), num_str_);
  Core::LinAlg::SerialDenseMatrix strain_data(stiffness_matrix_integration_.num_points(), num_str_);

  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, lm);

  CentroidTransformation<celltype> centroid_transformation =
      evaluate_centroid_transformation<celltype>(nodal_coordinates);

  evaluate_centroid_coordinates_and_add_to_parameter_list<celltype>(nodal_coordinates, params);

  Discret::Elements::for_each_gauss_point<celltype>(nodal_coordinates,
      stiffness_matrix_integration_,
      [&](const Core::LinAlg::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        const EASKinematics<celltype, eastype> kinematic_quantitites =
            evaluate_eas_kinematics<celltype, eastype, kinematic_type>(nodal_coordinates,
                centroid_transformation, xi, jacobian_mapping, eas_iteration_data_);

        evaluate_gp_coordinates_and_add_to_parameter_list<celltype>(
            nodal_coordinates, shape_functions, params);

        const Stress<celltype> stress = evaluate_material_stress<celltype>(solid_material,
            kinematic_quantitites.enhanced_deformation_gradient, kinematic_quantitites.enhanced_gl,
            params, gp, ele.id());

        assemble_strain_type_to_matrix_row<celltype>(kinematic_quantitites.enhanced_gl,
            kinematic_quantitites.enhanced_deformation_gradient, strainIO.type, strain_data, gp);
        assemble_stress_type_to_matrix_row(kinematic_quantitites.enhanced_deformation_gradient,
            stress, stressIO.type, stress_data, gp);
      });

  serialize(stress_data, serialized_stress_data);
  serialize(strain_data, serialized_strain_data);
}

template <Core::FE::CellType celltype, Discret::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
double
Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::calculate_internal_energy(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material,
    const Core::FE::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  double intenergy = 0.0;
  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, lm);

  CentroidTransformation<celltype> centroid_transformation =
      evaluate_centroid_transformation<celltype>(nodal_coordinates);

  Discret::Elements::for_each_gauss_point<celltype>(nodal_coordinates,
      stiffness_matrix_integration_,
      [&](const Core::LinAlg::Matrix<num_dim_, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        const EASKinematics<celltype, eastype> kinematic_quantitites =
            evaluate_eas_kinematics<celltype, eastype, kinematic_type>(nodal_coordinates,
                centroid_transformation, xi, jacobian_mapping, eas_iteration_data_);

        double psi = 0.0;
        solid_material.strain_energy(kinematic_quantitites.enhanced_gl, psi, gp, ele.id());

        intenergy += psi * integration_factor;
      });

  return intenergy;
}

template <Core::FE::CellType celltype, Discret::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::setup(
    Mat::So3Material& solid_material, const Core::IO::InputParameterContainer& container)
{
  solid_material.setup(stiffness_matrix_integration_.num_points(), container);
}

template <Core::FE::CellType celltype, Discret::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::material_post_setup(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material)
{
  Teuchos::ParameterList params{};

  // Check if element has fiber nodes, if so interpolate fibers to Gauss Points and add to params
  interpolate_fibers_to_gauss_points_and_add_to_parameter_list<celltype>(
      stiffness_matrix_integration_, ele, params);

  // Call post_setup of material
  solid_material.post_setup(params, ele.id());
}

template <Core::FE::CellType celltype, Discret::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype,
    kinematic_type>::initialize_gauss_point_data_output(const Core::Elements::Element& ele,
    const Mat::So3Material& solid_material,
    FourC::Solid::ModelEvaluator::GaussPointDataOutputManager& gp_data_output_manager) const
{
  FOUR_C_ASSERT(ele.is_params_interface(),
      "This action type should only be called from the new time integration framework!");

  ask_and_add_quantities_to_gauss_point_data_output(
      stiffness_matrix_integration_.num_points(), solid_material, gp_data_output_manager);
}

template <Core::FE::CellType celltype, Discret::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype,
    kinematic_type>::evaluate_gauss_point_data_output(const Core::Elements::Element& ele,
    const Mat::So3Material& solid_material,
    FourC::Solid::ModelEvaluator::GaussPointDataOutputManager& gp_data_output_manager) const
{
  FOUR_C_ASSERT(ele.is_params_interface(),
      "This action type should only be called from the new time integration framework!");

  collect_and_assemble_gauss_point_data_output<celltype>(
      stiffness_matrix_integration_, solid_material, ele, gp_data_output_manager);
}

template <Core::FE::CellType celltype, Discret::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::reset_to_last_converged(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material)
{
  solid_material.reset_step();
}

template <Core::FE::CellType celltype, Discret::Elements::EasType eastype,
    Inpar::Solid::KinemType kinematic_type>
void Discret::Elements::SolidEleCalcEas<celltype, eastype, kinematic_type>::for_each_gauss_point(
    const Core::Elements::Element& ele, Mat::So3Material& solid_material,
    const Core::FE::Discretization& discretization, const std::vector<int>& lm,
    const std::function<void(Mat::So3Material&, double, int)>& integrator) const
{
  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(ele, discretization, lm);

  Discret::Elements::for_each_gauss_point(nodal_coordinates, stiffness_matrix_integration_,
      [&](const Core::LinAlg::Matrix<Internal::num_dim<celltype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      { integrator(solid_material, integration_factor, gp); });
}

// template classes
template class Discret::Elements::SolidEleCalcEas<Core::FE::CellType::hex8,
    Discret::Elements::EasType::eastype_h8_9, Inpar::Solid::KinemType::nonlinearTotLag>;
template class Discret::Elements::SolidEleCalcEas<Core::FE::CellType::hex8,
    Discret::Elements::EasType::eastype_h8_21, Inpar::Solid::KinemType::nonlinearTotLag>;
template class Discret::Elements::SolidEleCalcEas<Core::FE::CellType::hex8,
    Discret::Elements::EasType::eastype_sh8_7, Inpar::Solid::KinemType::nonlinearTotLag>;
template class Discret::Elements::SolidEleCalcEas<Core::FE::CellType::hex8,
    Discret::Elements::EasType::eastype_h8_9, Inpar::Solid::KinemType::linear>;
template class Discret::Elements::SolidEleCalcEas<Core::FE::CellType::hex8,
    Discret::Elements::EasType::eastype_h8_21, Inpar::Solid::KinemType::linear>;

static_assert(
    Core::Communication::is_packable<Discret::Elements::SolidEleCalcEas<Core::FE::CellType::hex8,
        Discret::Elements::EasType::eastype_h8_9, Inpar::Solid::KinemType::nonlinearTotLag>>,
    "EAS needs to implement the method pack(Core::Communication::PackBuffer&) to be able to store "
    "history data!");
static_assert(
    Core::Communication::is_unpackable<Discret::Elements::SolidEleCalcEas<Core::FE::CellType::hex8,
        Discret::Elements::EasType::eastype_h8_9, Inpar::Solid::KinemType::nonlinearTotLag>>,
    "EAS needs to implement the method unpack(std::size_t, std::vector<char>&) to be able to store "
    "history data!");
static_assert(
    Core::Communication::is_packable<Discret::Elements::SolidEleCalcEas<Core::FE::CellType::hex8,
        Discret::Elements::EasType::eastype_h8_9, Inpar::Solid::KinemType::linear>>,
    "EAS needs to implement the method pack(Core::Communication::PackBuffer&) to be able to store "
    "history data!");
static_assert(
    Core::Communication::is_unpackable<Discret::Elements::SolidEleCalcEas<Core::FE::CellType::hex8,
        Discret::Elements::EasType::eastype_h8_9, Inpar::Solid::KinemType::linear>>,
    "EAS needs to implement the method unpack(std::size_t, std::vector<char>&) to be able to store "
    "history data!");
FOUR_C_NAMESPACE_CLOSE
