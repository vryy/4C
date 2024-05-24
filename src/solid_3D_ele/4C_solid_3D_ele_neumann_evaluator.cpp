/*! \file

\brief Evaluation of neumann loads

\level 1
*/

#include "4C_solid_3D_ele_neumann_evaluator.hpp"

#include "4C_discretization_fem_general_utils_gausspoints.hpp"
#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_element.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_solid_3D_ele_calc_lib_integration.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN

void DRT::ELEMENTS::EvaluateNeumannByElement(DRT::Element& element,
    const DRT::Discretization& discretization, CORE::Conditions::Condition& condition,
    const std::vector<int>& dof_index_array, CORE::LINALG::SerialDenseVector& element_force_vector,
    double total_time)
{
  switch (element.Shape())
  {
    case CORE::FE::CellType::hex8:
      return evaluate_neumann<CORE::FE::CellType::hex8>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case CORE::FE::CellType::hex27:
      return evaluate_neumann<CORE::FE::CellType::hex27>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case CORE::FE::CellType::hex20:
      return evaluate_neumann<CORE::FE::CellType::hex20>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case CORE::FE::CellType::hex18:
      return evaluate_neumann<CORE::FE::CellType::hex18>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case CORE::FE::CellType::nurbs27:
      return evaluate_neumann<CORE::FE::CellType::nurbs27>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case CORE::FE::CellType::pyramid5:
      return evaluate_neumann<CORE::FE::CellType::pyramid5>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case CORE::FE::CellType::wedge6:
      return evaluate_neumann<CORE::FE::CellType::wedge6>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case CORE::FE::CellType::tet4:
      return evaluate_neumann<CORE::FE::CellType::tet4>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case CORE::FE::CellType::tet10:
      return evaluate_neumann<CORE::FE::CellType::tet10>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    default:
      FOUR_C_THROW(
          "The cell type you are trying to evaluate the Neumann condition for is not yet "
          "implemented in EvaluateNeumannByElement.");
      break;
  }
}

template <CORE::FE::CellType celltype>
void DRT::ELEMENTS::evaluate_neumann(DRT::Element& element,
    const DRT::Discretization& discretization, CORE::Conditions::Condition& condition,
    const std::vector<int>& dof_index_array, CORE::LINALG::SerialDenseVector& element_force_vector,
    double total_time)
{
  constexpr auto numdim = CORE::FE::dim<celltype>;
  constexpr auto numnod = CORE::FE::num_nodes<celltype>;
  CORE::FE::GaussIntegration gauss_integration =
      CreateGaussIntegration<celltype>(DRT::ELEMENTS::GetGaussRuleStiffnessMatrix<celltype>());

  // get values and switches from the condition
  const auto& onoff = condition.parameters().Get<std::vector<int>>("onoff");
  const auto& value = condition.parameters().Get<std::vector<double>>("val");

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
  const auto& function_ids = condition.parameters().Get<std::vector<int>>("funct");

  const ElementNodes<celltype> nodal_coordinates =
      EvaluateElementNodes<celltype>(element, discretization, dof_index_array);

  ForEachGaussPoint<celltype>(nodal_coordinates, gauss_integration,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<celltype>, 1>& xi,
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
        CORE::LINALG::Matrix<numdim, 1> gauss_point_reference_coordinates;
        gauss_point_reference_coordinates.MultiplyTN(
            nodal_coordinates.reference_coordinates_, shape_functions.shapefunctions_);

        for (auto dim = 0; dim < numdim; dim++)
        {
          if (onoff[dim])
          {
            // function evaluation
            const int function_number = function_ids[dim];
            const double function_scale_factor =
                (function_number > 0)
                    ? GLOBAL::Problem::Instance()
                          ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(function_number - 1)
                          .Evaluate(gauss_point_reference_coordinates.A(), total_time, dim)
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
