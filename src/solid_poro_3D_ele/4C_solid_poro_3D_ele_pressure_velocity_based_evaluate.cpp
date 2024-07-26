/*! \file

\brief Evaluation routines for the solid-poro pressure-velocity-based element

This file contains the element-specific evaluation routines such as
evaluate(...), evaluate_neumann(...), etc.

\level 1
*/

#include "4C_mat_structporo.hpp"
#include "4C_solid_3D_ele_calc_lib.hpp"
#include "4C_solid_poro_3D_ele_pressure_velocity_based.hpp"

FOUR_C_NAMESPACE_OPEN


int Discret::ELEMENTS::SolidPoroPressureVelocityBased::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Elements::Element::LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  if (!material_post_setup_)
  {
    std::visit([&](auto& interface)
        { interface->material_post_setup(*this, struct_poro_material()); },
        solid_calc_variant_);
    material_post_setup_ = true;
  }

  // get ptr to interface to time integration
  set_params_interface_ptr(params);

  const Core::Elements::ActionType action = std::invoke(
      [&]()
      {
        if (is_params_interface())
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
            interface->evaluate_nonlinear_force_stiffness_mass(*this, this->struct_poro_material(),
                discretization, la[0].lm_, params, &elevec1, &elemat1, nullptr);
          },
          solid_calc_variant_);

      if (la.size() > 1)
      {
        if (discretization.has_state(1, "fluidvel"))
        {
          std::visit(
              [&](auto& interface)
              {
                interface->evaluate_nonlinear_force_stiffness(*this, this->struct_poro_material(),
                    this->fluid_poro_material(), this->get_anisotropic_permeability_property(),
                    this->kinematic_type(), discretization, la, params, &elevec1, &elemat1,
                    &elemat2);
              },
              solidporo_press_vel_based_calc_variant_);
        }
      }
      return 0;
    }
    case Core::Elements::struct_calc_internalforce:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->evaluate_nonlinear_force_stiffness_mass(*this, this->struct_poro_material(),
                discretization, la[0].lm_, params, &elevec1, nullptr, nullptr);
          },
          solid_calc_variant_);

      if (la.size() > 1)
      {
        if (discretization.has_state(1, "fluidvel"))
        {
          std::visit(
              [&](auto& interface)
              {
                interface->evaluate_nonlinear_force_stiffness(*this, this->struct_poro_material(),
                    this->fluid_poro_material(), this->get_anisotropic_permeability_property(),
                    this->kinematic_type(), discretization, la, params, &elevec1, nullptr, nullptr);
              },
              solidporo_press_vel_based_calc_variant_);
        }
      }
      return 0;
    }
    case Core::Elements::struct_calc_nlnstiffmass:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->evaluate_nonlinear_force_stiffness_mass(*this, this->struct_poro_material(),
                discretization, la[0].lm_, params, &elevec1, &elemat1, &elemat2);
          },
          solid_calc_variant_);

      if (la.size() > 1 and this->num_material() > 1)
      {
        if (discretization.has_state(1, "fluidvel"))
        {
          std::visit(
              [&](auto& interface)
              {
                interface->evaluate_nonlinear_force_stiffness(*this, this->struct_poro_material(),
                    this->fluid_poro_material(), this->get_anisotropic_permeability_property(),
                    this->kinematic_type(), discretization, la, params, &elevec1, &elemat1,
                    nullptr);
              },
              solidporo_press_vel_based_calc_variant_);
        }
      }
      return 0;
    }
    case Core::Elements::struct_calc_nlnstifflmass:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->evaluate_nonlinear_force_stiffness_mass(*this, this->struct_poro_material(),
                discretization, la[0].lm_, params, &elevec1, &elemat1, &elemat2);
          },
          solid_calc_variant_);
      Discret::ELEMENTS::LumpMatrix(elemat2);
      return 0;
    }
    case Core::Elements::struct_poro_calc_scatracoupling:
    {
      // no coupling -> return
      return 0;
    }
    case Core::Elements::struct_poro_calc_fluidcoupling:
    {
      if (discretization.has_state(1, "fluidvel"))
      {
        std::visit(
            [&](auto& interface)
            {
              interface->evaluate_nonlinear_force_stiffness_od(*this, this->struct_poro_material(),
                  this->fluid_poro_material(), this->get_anisotropic_permeability_property(),
                  this->kinematic_type(), discretization, la, params, &elemat1);
            },
            solidporo_press_vel_based_calc_variant_);
      }
      return 0;
    }
    case Core::Elements::struct_calc_update_istep:
    {
      std::visit([&](auto& interface)
          { interface->update(*this, solid_poro_material(), discretization, la[0].lm_, params); },
          solid_calc_variant_);
      return 0;
    }
    case Core::Elements::struct_calc_recover:
    {
      std::visit([&](auto& interface)
          { interface->recover(*this, discretization, la[0].lm_, params); },
          solid_calc_variant_);
      return 0;
    }
    case Core::Elements::struct_calc_stress:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->calculate_stress(*this, this->struct_poro_material(),
                StressIO{get_io_stress_type(*this, params), get_stress_data(*this, params)},
                StrainIO{get_io_strain_type(*this, params), get_strain_data(*this, params)},
                discretization, la[0].lm_, params);
          },
          solid_calc_variant_);


      if (discretization.has_state(1, "fluidvel"))
      {
        std::visit(
            [&](auto& interface)
            {
              interface->coupling_stress_poroelast(*this, this->struct_poro_material(),
                  this->kinematic_type(),
                  CouplStressIO{
                      get_io_couplstress_type(*this, params), get_couplstress_data(*this, params)},
                  discretization, la, params);
            },
            solidporo_press_vel_based_calc_variant_);
      }



      return 0;
    }
    case Core::Elements::struct_init_gauss_point_data_output:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->initialize_gauss_point_data_output(*this, solid_poro_material(),
                *params_interface().gauss_point_data_output_manager_ptr());
          },
          solid_calc_variant_);
      return 0;
    }
    case Core::Elements::struct_gauss_point_data_output:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->evaluate_gauss_point_data_output(*this, solid_poro_material(),
                *params_interface().gauss_point_data_output_manager_ptr());
          },
          solid_calc_variant_);
      return 0;
    }
    case Core::Elements::struct_calc_predict:
    case Core::Elements::none:
    {
      // do nothing for now
      return 0;
    }
    default:
      FOUR_C_THROW("The element action %s is not yet implemented for the new solid elements",
          ActionType2String(action).c_str());
      // do nothing (no error because there are some actions the poro element is supposed to ignore)
      return 0;
  }
}

int Discret::ELEMENTS::SolidPoroPressureVelocityBased::evaluate_neumann(
    Teuchos::ParameterList& params, Core::FE::Discretization& discretization,
    Core::Conditions::Condition& condition, std::vector<int>& lm,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseMatrix* elemat1)
{
  FOUR_C_THROW(
      "Cannot yet evaluate volume neumann forces within the pressure-velocity based poro "
      "implementation.");
}
FOUR_C_NAMESPACE_CLOSE
