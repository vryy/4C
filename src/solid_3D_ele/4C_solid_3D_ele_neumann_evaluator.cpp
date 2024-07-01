/*! \file

\brief Evaluation of neumann loads

\level 1
*/

#include "4C_solid_3D_ele_neumann_evaluator.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_global_data.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_solid_3D_ele_calc_lib_integration.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN

void Discret::ELEMENTS::evaluate_neumann_by_element(Core::Elements::Element& element,
    const Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    const std::vector<int>& dof_index_array, Core::LinAlg::SerialDenseVector& element_force_vector,
    double total_time)
{
  switch (element.Shape())
  {
    case Core::FE::CellType::hex8:
      return evaluate_neumann<Core::FE::CellType::hex8>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case Core::FE::CellType::hex27:
      return evaluate_neumann<Core::FE::CellType::hex27>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case Core::FE::CellType::hex20:
      return evaluate_neumann<Core::FE::CellType::hex20>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case Core::FE::CellType::hex18:
      return evaluate_neumann<Core::FE::CellType::hex18>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case Core::FE::CellType::nurbs27:
      return evaluate_neumann<Core::FE::CellType::nurbs27>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case Core::FE::CellType::pyramid5:
      return evaluate_neumann<Core::FE::CellType::pyramid5>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case Core::FE::CellType::wedge6:
      return evaluate_neumann<Core::FE::CellType::wedge6>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case Core::FE::CellType::tet4:
      return evaluate_neumann<Core::FE::CellType::tet4>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case Core::FE::CellType::tet10:
      return evaluate_neumann<Core::FE::CellType::tet10>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    default:
      FOUR_C_THROW(
          "The cell type you are trying to evaluate the Neumann condition for is not yet "
          "implemented in EvaluateNeumannByElement.");
      break;
  }
}

template <Core::FE::CellType celltype>
void Discret::ELEMENTS::evaluate_neumann(Core::Elements::Element& element,
    const Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    const std::vector<int>& dof_index_array, Core::LinAlg::SerialDenseVector& element_force_vector,
    double total_time)
{
  constexpr auto numdim = Core::FE::dim<celltype>;
  constexpr auto numnod = Core::FE::num_nodes<celltype>;
  Core::FE::GaussIntegration gauss_integration = create_gauss_integration<celltype>(
      Discret::ELEMENTS::get_gauss_rule_stiffness_matrix<celltype>());

  // get values and switches from the condition
  const auto& onoff = condition.parameters().get<std::vector<int>>("onoff");
  const auto& value = condition.parameters().get<std::vector<double>>("val");

  // ensure that at least as many curves/functs as dofs are available
  if (onoff.size() < numdim)
    FOUR_C_THROW("Fewer functions or curves defined than the element's dimension.");

  for (std::size_t checkdof = numdim; checkdof < onoff.size(); ++checkdof)
  {
    if (onoff[checkdof] != 0)
    {
      FOUR_C_THROW(
          "You have activated more than %d dofs in your Neumann boundary condition. This is higher "
          "than the dimension of the element.",
          numdim);
    }
  }

  // get ids of functions of space and time
  const auto& function_ids = condition.parameters().get<std::vector<int>>("funct");

  const ElementNodes<celltype> nodal_coordinates =
      evaluate_element_nodes<celltype>(element, discretization, dof_index_array);

  ForEachGaussPoint<celltype>(nodal_coordinates, gauss_integration,
      [&](const Core::LinAlg::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<celltype>& shape_functions,
          const JacobianMapping<celltype>& jacobian_mapping, double integration_factor, int gp)
      {
        if (jacobian_mapping.determinant_ == 0.0)
          FOUR_C_THROW(
              "The determinant of the jacobian is zero for element with id %i", element.Id());
        else if (jacobian_mapping.determinant_ < 0.0)
          FOUR_C_THROW("The determinant of the jacobian is negative (%d) for element with id %i",
              jacobian_mapping.determinant_, element.Id());

        // material/reference co-ordinates of Gauss point
        Core::LinAlg::Matrix<numdim, 1> gauss_point_reference_coordinates;
        gauss_point_reference_coordinates.multiply_tn(
            nodal_coordinates.reference_coordinates, shape_functions.shapefunctions_);

        for (auto dim = 0; dim < numdim; dim++)
        {
          if (onoff[dim])
          {
            // function evaluation
            const int function_number = function_ids[dim];
            const double function_scale_factor =
                (function_number > 0)
                    ? Global::Problem::Instance()
                          ->FunctionById<Core::UTILS::FunctionOfSpaceTime>(function_number - 1)
                          .evaluate(gauss_point_reference_coordinates.data(), total_time, dim)
                    : 1.0;

            const double value_times_integration_factor =
                value[dim] * function_scale_factor * integration_factor;

            for (auto nodeid = 0; nodeid < numnod; ++nodeid)
            {
              // Evaluates the Neumann boundary condition: f_{x,y,z}^i=\sum_j N^i(xi^j) * value(t) *
              // integration_factor_j
              // assembles the element force vector [f_x^1, f_y^1, f_z^1, ..., f_x^n, f_y^n, f_z^n]
              element_force_vector[nodeid * numdim + dim] +=
                  shape_functions.shapefunctions_(nodeid) * value_times_integration_factor;
            }
          }
        }
      });
}
FOUR_C_NAMESPACE_CLOSE
