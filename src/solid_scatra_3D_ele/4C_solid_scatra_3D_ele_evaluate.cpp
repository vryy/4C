/*! \file

\brief Evaluation routines for the solid-scatra element

This file contains the element-specific evaluation routines such as
Evaluate(...), evaluate_neumann(...), etc.

\level 1
*/

#include "4C_discretization_fem_general_elements_paramsinterface.hpp"
#include "4C_solid_3D_ele_neumann_evaluator.hpp"
#include "4C_solid_scatra_3D_ele.hpp"
#include "4C_solid_scatra_3D_ele_calc_lib_nitsche.hpp"

FOUR_C_NAMESPACE_OPEN

int Discret::ELEMENTS::SolidScatra::Evaluate(Teuchos::ParameterList& params,
    Discret::Discretization& discretization, Core::Elements::Element::LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  if (!material_post_setup_)
  {
    std::visit([&](auto& interface) { interface->material_post_setup(*this, SolidMaterial()); },
        solid_scatra_calc_variant_);
    material_post_setup_ = true;
  }

  // get ptr to interface to time integration
  set_params_interface_ptr(params);

  const Core::Elements::ActionType action = std::invoke(
      [&]()
      {
        if (IsParamsInterface())
          return params_interface().GetActionType();
        else
          return Core::Elements::String2ActionType(params.get<std::string>("action", "none"));
      });

  switch (action)
  {
    case Core::Elements::ActionType::calc_struct_stiffscalar:
      // Evaluate off-diagonal block for SSI
      std::visit(
          [&](auto& interface)
          {
            interface->evaluate_d_stress_d_scalar(
                *this, SolidMaterial(), discretization, la, params, elemat1);
          },
          solid_scatra_calc_variant_);
      return 0;
    case Core::Elements::struct_calc_nlnstiff:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->evaluate_nonlinear_force_stiffness_mass(
                *this, SolidMaterial(), discretization, la, params, &elevec1, &elemat1, nullptr);
          },
          solid_scatra_calc_variant_);
      return 0;
    }
    case Core::Elements::struct_calc_nlnstiffmass:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->evaluate_nonlinear_force_stiffness_mass(
                *this, SolidMaterial(), discretization, la, params, &elevec1, &elemat1, &elemat2);
          },
          solid_scatra_calc_variant_);

      return 0;
    }
    case Core::Elements::struct_calc_internalforce:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->evaluate_nonlinear_force_stiffness_mass(
                *this, SolidMaterial(), discretization, la, params, &elevec1, nullptr, nullptr);
          },
          solid_scatra_calc_variant_);

      return 0;
    }
    case Core::Elements::struct_calc_update_istep:
    {
      std::visit([&](auto& interface)
          { interface->Update(*this, SolidMaterial(), discretization, la, params); },
          solid_scatra_calc_variant_);

      return 0;
    }
    case Core::Elements::struct_calc_recover:
    {
      std::visit([&](auto& interface) { interface->Recover(*this, discretization, la, params); },
          solid_scatra_calc_variant_);

      return 0;
    }
    case Core::Elements::struct_calc_stress:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->CalculateStress(*this, SolidMaterial(),
                StressIO{GetIOStressType(*this, params), GetStressData(*this, params)},
                StrainIO{GetIOStrainType(*this, params), GetStrainData(*this, params)},
                discretization, la, params);
          },
          solid_scatra_calc_variant_);

      return 0;
    }
    case Core::Elements::struct_calc_predict:
      // do nothing for now
      return 0;
    default:
      FOUR_C_THROW("The element action %s is not yet implemented for the new solid-scatra elements",
          ActionType2String(action).c_str());
  }

  return 0;
}

int Discret::ELEMENTS::SolidScatra::evaluate_neumann(Teuchos::ParameterList& params,
    Discret::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  set_params_interface_ptr(params);

  const double time = std::invoke(
      [&]()
      {
        if (IsParamsInterface())
          return params_interface().GetTotalTime();
        else
          return params.get("total time", -1.0);
      });

  Discret::ELEMENTS::EvaluateNeumannByElement(*this, discretization, condition, lm, elevec1, time);
  return 0;
}

double Discret::ELEMENTS::SolidScatra::GetCauchyNDirAtXi(const std::vector<double>& disp,
    const std::optional<std::vector<double>>& scalars, const Core::LinAlg::Matrix<3, 1>& xi,
    const Core::LinAlg::Matrix<3, 1>& n, const Core::LinAlg::Matrix<3, 1>& dir,
    Discret::ELEMENTS::SolidScatraCauchyNDirLinearizations<3>& linearizations)
{
  return Discret::ELEMENTS::GetCauchyNDirAtXi(solid_scatra_calc_variant_, *this, SolidMaterial(),
      disp, scalars, xi, n, dir, linearizations);
}

FOUR_C_NAMESPACE_CLOSE