/*! \file

\brief Implementation of the solid element

This file contains the element-specific evaluation routines such as
evaluate(...), evaluate_neumann(...), etc.

\level 1
*/

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_elements_paramsinterface.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_solid_3D_ele.hpp"
#include "4C_solid_3D_ele_calc_interface.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_solid_3D_ele_calc_lib_io.hpp"
#include "4C_solid_3D_ele_calc_lib_nitsche.hpp"
#include "4C_solid_3D_ele_calc_mulf.hpp"
#include "4C_solid_3D_ele_neumann_evaluator.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  std::vector<double> get_acceleration_vector(
      const Core::FE::Discretization& discretization, const std::vector<int>& lm)
  {
    const Epetra_Vector& acceleration = *discretization.GetState("acceleration");
    std::vector<double> my_acceleration(lm.size());
    Core::FE::ExtractMyValues(acceleration, my_acceleration, lm);

    return my_acceleration;
  }

  void evaluate_inertia_force(const Core::LinAlg::SerialDenseMatrix& mass_matrix,
      const Core::LinAlg::SerialDenseVector& acceleration,
      Core::LinAlg::SerialDenseVector& inertia_force)
  {
    inertia_force.putScalar();
    inertia_force.multiply(Teuchos::ETransp::NO_TRANS, Teuchos::ETransp::NO_TRANS, 1.0, mass_matrix,
        acceleration, 0.0);
  }

  void evaluate_inertia_force(const Core::FE::Discretization& discretization,
      const std::vector<int>& lm, const Core::LinAlg::SerialDenseMatrix& mass_matrix,
      Core::LinAlg::SerialDenseVector& inertia_force)
  {
    std::vector<double> my_acceleration = get_acceleration_vector(discretization, lm);

    Core::LinAlg::SerialDenseVector acceleration_vector(
        Teuchos::DataAccess::View, my_acceleration.data(), static_cast<int>(lm.size()));

    evaluate_inertia_force(mass_matrix, acceleration_vector, inertia_force);
  }
}  // namespace

int Discret::ELEMENTS::Solid::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  if (!material_post_setup_)
  {
    std::visit([&](auto& interface) { interface->material_post_setup(*this, *SolidMaterial()); },
        solid_calc_variant_);
    material_post_setup_ = true;
  }

  // get ptr to interface to time integration
  set_params_interface_ptr(params);

  const Core::Elements::ActionType action = std::invoke(
      [&]()
      {
        if (IsParamsInterface())
          return params_interface().get_action_type();
        else
          return Core::Elements::String2ActionType(params.get<std::string>("action", "none"));
      });

  switch (action)
  {
    case Core::Elements::struct_calc_nlnstiff:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->evaluate_nonlinear_force_stiffness_mass(
                *this, *SolidMaterial(), discretization, lm, params, &elevec1, &elemat1, nullptr);
          },
          solid_calc_variant_);
      return 0;
    }
    case Core::Elements::struct_calc_internalforce:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->evaluate_nonlinear_force_stiffness_mass(
                *this, *SolidMaterial(), discretization, lm, params, &elevec1, nullptr, nullptr);
          },
          solid_calc_variant_);

      return 0;
    }
    case Core::Elements::struct_calc_nlnstiffmass:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->evaluate_nonlinear_force_stiffness_mass(
                *this, *SolidMaterial(), discretization, lm, params, &elevec1, &elemat1, &elemat2);
          },
          solid_calc_variant_);

      evaluate_inertia_force(discretization, lm, elemat2, elevec2);
      return 0;
    }
    case Core::Elements::struct_calc_nlnstifflmass:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->evaluate_nonlinear_force_stiffness_mass(
                *this, *SolidMaterial(), discretization, lm, params, &elevec1, &elemat1, &elemat2);
          },
          solid_calc_variant_);

      LumpMatrix(elemat2);

      evaluate_inertia_force(discretization, lm, elemat2, elevec2);
      return 0;
    }
    case Core::Elements::struct_calc_internalinertiaforce:
    {
      const int num_dof_per_ele = static_cast<int>(lm.size());
      Core::LinAlg::SerialDenseMatrix mass_matrix(num_dof_per_ele, num_dof_per_ele);

      std::visit(
          [&](auto& interface)
          {
            interface->evaluate_nonlinear_force_stiffness_mass(*this, *SolidMaterial(),
                discretization, lm, params, &elevec1, nullptr, &mass_matrix);
          },
          solid_calc_variant_);

      evaluate_inertia_force(discretization, lm, mass_matrix, elevec2);
      return 0;
    }
    case Core::Elements::struct_calc_update_istep:
    {
      std::visit([&](auto& interface)
          { interface->Update(*this, *SolidMaterial(), discretization, lm, params); },
          solid_calc_variant_);

      return 0;
    }
    case Core::Elements::struct_update_prestress:
    {
      update_prestress(solid_calc_variant_, *this, *SolidMaterial(), discretization, lm, params);
      return 0;
    }
    case Core::Elements::struct_calc_recover:
    {
      std::visit([&](auto& interface) { interface->Recover(*this, discretization, lm, params); },
          solid_calc_variant_);

      return 0;
    }
    case Core::Elements::struct_calc_stress:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->calculate_stress(*this, *SolidMaterial(),
                StressIO{get_io_stress_type(*this, params), get_stress_data(*this, params)},
                StrainIO{get_io_strain_type(*this, params), get_strain_data(*this, params)},
                discretization, lm, params);
          },
          solid_calc_variant_);

      return 0;
    }
    case Core::Elements::struct_calc_energy:
    {
      double int_energy = std::visit(
          [&](auto& interface)
          {
            return interface->calculate_internal_energy(
                *this, *SolidMaterial(), discretization, lm, params);
          },
          solid_calc_variant_);


      if (IsParamsInterface())
      {
        // new structural time integration
        params_interface().add_contribution_to_energy_type(int_energy, STR::internal_energy);
      }
      else
      {
        // old structural time integration
        // check length of elevec1
        if (elevec1.length() < 1) FOUR_C_THROW("The given result vector is too short.");

        elevec1(0) = int_energy;
      }
      return 0;
    }
    case Core::Elements::struct_init_gauss_point_data_output:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->initialize_gauss_point_data_output(
                *this, *SolidMaterial(), *params_interface().gauss_point_data_output_manager_ptr());
          },
          solid_calc_variant_);

      return 0;
    }
    case Core::Elements::struct_gauss_point_data_output:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->evaluate_gauss_point_data_output(
                *this, *SolidMaterial(), *params_interface().gauss_point_data_output_manager_ptr());
          },
          solid_calc_variant_);

      return 0;
    }
    case Core::Elements::struct_calc_reset_istep:
    {
      std::visit([&](auto& interface)
          { interface->reset_to_last_converged(*this, *SolidMaterial()); },
          solid_calc_variant_);

      return 0;
    }
    case Core::Elements::struct_calc_predict:
    {
      // do nothing for now
      return 0;
    }
    default:
      FOUR_C_THROW("The element action %s is not yet implemented for the new solid elements",
          ActionType2String(action).c_str());
  }

  return 0;
}
int Discret::ELEMENTS::Solid::evaluate_neumann(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  set_params_interface_ptr(params);

  const double time = std::invoke(
      [&]()
      {
        if (IsParamsInterface())
          return params_interface().get_total_time();
        else
          return params.get("total time", -1.0);
      });

  Discret::ELEMENTS::evaluate_neumann_by_element(
      *this, discretization, condition, lm, elevec1, time);
  return 0;
}

template <int dim>
double Discret::ELEMENTS::Solid::get_normal_cauchy_stress_at_xi(const std::vector<double>& disp,
    const Core::LinAlg::Matrix<dim, 1>& xi, const Core::LinAlg::Matrix<dim, 1>& n,
    const Core::LinAlg::Matrix<dim, 1>& dir, CauchyNDirLinearizations<dim>& linearizations)
{
  return Discret::ELEMENTS::get_normal_cauchy_stress_at_xi<dim>(
      solid_calc_variant_, *this, *SolidMaterial(), disp, xi, n, dir, linearizations);
}

template double Discret::ELEMENTS::Solid::get_normal_cauchy_stress_at_xi<3>(
    const std::vector<double>&, const Core::LinAlg::Matrix<3, 1>&,
    const Core::LinAlg::Matrix<3, 1>&, const Core::LinAlg::Matrix<3, 1>&,
    CauchyNDirLinearizations<3>&);
template double Discret::ELEMENTS::Solid::get_normal_cauchy_stress_at_xi<2>(
    const std::vector<double>&, const Core::LinAlg::Matrix<2, 1>&,
    const Core::LinAlg::Matrix<2, 1>&, const Core::LinAlg::Matrix<2, 1>&,
    CauchyNDirLinearizations<2>&);

FOUR_C_NAMESPACE_CLOSE
