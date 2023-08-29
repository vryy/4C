/*! \file

\brief Evaluation of neumann loads

\level 1
*/

#include "baci_solid_ele_neumann_evaluator.H"

#include "baci_discretization_fem_general_utils_gausspoints.H"
#include "baci_discretization_fem_general_utils_local_connectivity_matrices.H"
#include "baci_lib_element.H"
#include "baci_lib_function.H"
#include "baci_lib_globalproblem.H"
#include "baci_solid_ele_calc_lib.H"

void DRT::ELEMENTS::EvaluateNeumannByElement(DRT::Element& element,
    const DRT::Discretization& discretization, DRT::Condition& condition,
    const std::vector<int>& dof_index_array, CORE::LINALG::SerialDenseVector& element_force_vector,
    double total_time)
{
  switch (element.Shape())
  {
    case DRT::Element::hex8:
      return EvaluateNeumann<DRT::Element::hex8>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case DRT::Element::hex27:
      return EvaluateNeumann<DRT::Element::hex27>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case DRT::Element::hex20:
      return EvaluateNeumann<DRT::Element::hex20>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case DRT::Element::hex18:
      return EvaluateNeumann<DRT::Element::hex18>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case DRT::Element::pyramid5:
      return EvaluateNeumann<DRT::Element::pyramid5>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case DRT::Element::wedge6:
      return EvaluateNeumann<DRT::Element::wedge6>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case DRT::Element::tet4:
      return EvaluateNeumann<DRT::Element::tet4>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    case DRT::Element::tet10:
      return EvaluateNeumann<DRT::Element::tet10>(
          element, discretization, condition, dof_index_array, element_force_vector, total_time);
      break;
    default:
      dserror(
          "The discretization type you are trying to evaluate the Neumann condition for is not yet "
          "implemented in EvaluateNeumannByElement.");
      break;
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::EvaluateNeumann(DRT::Element& element,
    const DRT::Discretization& discretization, DRT::Condition& condition,
    const std::vector<int>& dof_index_array, CORE::LINALG::SerialDenseVector& element_force_vector,
    double total_time)
{
  constexpr auto numdim = CORE::DRT::UTILS::DisTypeToDim<distype>::dim;
  constexpr auto numnod = CORE::DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
  CORE::DRT::UTILS::GaussIntegration gauss_integration =
      CreateGaussIntegration<distype>(DRT::ELEMENTS::GetGaussRuleStiffnessMatrix<distype>());

  // get values and switches from the condition
  const auto& onoff = *condition.Get<std::vector<int>>("onoff");
  const auto& value = *condition.Get<std::vector<double>>("val");

  // ensure that at least as many curves/functs as dofs are available
  if (onoff.size() < numdim)
    dserror("Fewer functions or curves defined than the element's dimension.");

  for (std::size_t checkdof = numdim; checkdof < onoff.size(); ++checkdof)
  {
    if (onoff[checkdof] != 0)
    {
      dserror(
          "You have activated more than %d dofs in your Neumann boundary condition. This is higher "
          "than the dimension of the element.",
          numdim);
    }
  }

  // get ids of functions of space and time
  const auto* function_ids = condition.Get<std::vector<int>>("funct");

  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(element, discretization, dof_index_array);

  ForEachGaussPoint<distype>(nodal_coordinates, gauss_integration,
      [&](const CORE::LINALG::Matrix<DETAIL::num_dim<distype>, 1>& xi,
          const ShapeFunctionsAndDerivatives<distype>& shape_functions,
          const JacobianMapping<distype>& jacobian_mapping, double integration_factor, int gp)
      {
        if (jacobian_mapping.determinant_ == 0.0)
          dserror("The determinant of the jacobian is zero for element with id %i", element.Id());
        else if (jacobian_mapping.determinant_ < 0.0)
          dserror("The determinant of the jacobian is negative (%d) for element with id %i",
              jacobian_mapping.determinant_, element.Id());

        // material/reference co-ordinates of Gauss point
        CORE::LINALG::Matrix<numdim, 1> gauss_point_reference_coordinates;
        gauss_point_reference_coordinates.MultiplyTN(
            nodal_coordinates.reference_, shape_functions.shapefunctions_);

        for (auto dim = 0; dim < numdim; dim++)
        {
          if (onoff[dim])
          {
            // function evaluation
            const int function_number = (function_ids != nullptr) ? (*function_ids)[dim] : -1;
            const double function_scale_factor =
                (function_number > 0)
                    ? DRT::Problem::Instance()
                          ->FunctionById<DRT::UTILS::FunctionOfSpaceTime>(function_number - 1)
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