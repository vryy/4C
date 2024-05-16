/*! \file

\brief Evaluation routines for the solid-scatra element

This file contains the element-specific evaluation routines such as
Evaluate(...), EvaluateNeumann(...), etc.

\level 1
*/

#include "4C_lib_elements_paramsinterface.hpp"
#include "4C_solid_3D_ele_neumann_evaluator.hpp"
#include "4C_solid_scatra_3D_ele.hpp"

FOUR_C_NAMESPACE_OPEN

int DRT::ELEMENTS::SolidScatra::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  if (!material_post_setup_)
  {
    std::visit([&](auto& interface) { interface->MaterialPostSetup(*this, SolidMaterial()); },
        solid_scatra_calc_variant_);
    material_post_setup_ = true;
  }

  // get ptr to interface to time integration
  SetParamsInterfacePtr(params);

  const ELEMENTS::ActionType action = std::invoke(
      [&]()
      {
        if (IsParamsInterface())
          return ParamsInterface().GetActionType();
        else
          return String2ActionType(params.get<std::string>("action", "none"));
      });

  switch (action)
  {
    case ELEMENTS::ActionType::calc_struct_stiffscalar:
      // Evaluate off-diagonal block for SSI
      std::visit(
          [&](auto& interface) {
            interface->EvaluateDStressDScalar(
                *this, SolidMaterial(), discretization, la, params, elemat1);
          },
          solid_scatra_calc_variant_);
      return 0;
    case DRT::ELEMENTS::struct_calc_nlnstiff:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->EvaluateNonlinearForceStiffnessMass(
                *this, SolidMaterial(), discretization, la, params, &elevec1, &elemat1, nullptr);
          },
          solid_scatra_calc_variant_);
      return 0;
    }
    case struct_calc_nlnstiffmass:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->EvaluateNonlinearForceStiffnessMass(
                *this, SolidMaterial(), discretization, la, params, &elevec1, &elemat1, &elemat2);
          },
          solid_scatra_calc_variant_);

      return 0;
    }
    case struct_calc_internalforce:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->EvaluateNonlinearForceStiffnessMass(
                *this, SolidMaterial(), discretization, la, params, &elevec1, nullptr, nullptr);
          },
          solid_scatra_calc_variant_);

      return 0;
    }
    case DRT::ELEMENTS::struct_calc_update_istep:
    {
      std::visit([&](auto& interface)
          { interface->Update(*this, SolidMaterial(), discretization, la, params); },
          solid_scatra_calc_variant_);

      return 0;
    }
    case DRT::ELEMENTS::struct_calc_recover:
    {
      std::visit([&](auto& interface) { interface->Recover(*this, discretization, la, params); },
          solid_scatra_calc_variant_);

      return 0;
    }
    case struct_calc_stress:
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
    case DRT::ELEMENTS::struct_calc_predict:
      // do nothing for now
      return 0;
    default:
      FOUR_C_THROW("The element action %s is not yet implemented for the new solid-scatra elements",
          ActionType2String(action).c_str());
  }

  return 0;
}

int DRT::ELEMENTS::SolidScatra::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, CORE::Conditions::Condition& condition,
    std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseMatrix* elemat1)
{
  SetParamsInterfacePtr(params);

  const double time = std::invoke(
      [&]()
      {
        if (IsParamsInterface())
          return ParamsInterface().GetTotalTime();
        else
          return params.get("total time", -1.0);
      });

  DRT::ELEMENTS::EvaluateNeumannByElement(*this, discretization, condition, lm, elevec1, time);
  return 0;
}
FOUR_C_NAMESPACE_CLOSE