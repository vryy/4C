/*! \file

\brief Implementation of the solid element

This file contains the element-specific evaluation routines such as
Evaluate(...), EvaluateNeumann(...), etc.

\level 1
*/

#include "baci_lib_elements_paramsinterface.H"
#include "baci_solid_ele.H"
#include "baci_solid_ele_calc_interface.H"
#include "baci_solid_ele_calc_lib.H"
#include "baci_solid_ele_calc_lib_io.H"
#include "baci_solid_ele_calc_mulf.H"
#include "baci_solid_ele_neumann_evaluator.H"
#include "baci_structure_new_elements_paramsinterface.H"
#include "baci_utils_exceptions.H"

BACI_NAMESPACE_OPEN

int DRT::ELEMENTS::Solid::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  if (!material_post_setup_)
  {
    std::visit([&](auto& interface) { interface->MaterialPostSetup(*this, *SolidMaterial()); },
        solid_calc_variant_);
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
    case DRT::ELEMENTS::struct_calc_nlnstiff:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->EvaluateNonlinearForceStiffnessMass(
                *this, *SolidMaterial(), discretization, lm, params, &elevec1, &elemat1, nullptr);
          },
          solid_calc_variant_);
      return 0;
    }
    case DRT::ELEMENTS::struct_calc_nlnstiff_gemm:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->EvaluateNonlinearForceStiffnessMassGEMM(
                *this, *SolidMaterial(), discretization, lm, params, &elevec1, &elemat1, nullptr);
          },
          solid_calc_variant_);

      return 0;
    }
    case struct_calc_internalforce:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->EvaluateNonlinearForceStiffnessMass(
                *this, *SolidMaterial(), discretization, lm, params, &elevec1, nullptr, nullptr);
          },
          solid_calc_variant_);

      return 0;
    }
    case struct_calc_nlnstiffmass:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->EvaluateNonlinearForceStiffnessMass(
                *this, *SolidMaterial(), discretization, lm, params, &elevec1, &elemat1, &elemat2);
          },
          solid_calc_variant_);

      return 0;
    }
    case struct_calc_nlnstifflmass:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->EvaluateNonlinearForceStiffnessMass(
                *this, *SolidMaterial(), discretization, lm, params, &elevec1, &elemat1, &elemat2);
          },
          solid_calc_variant_);

      LumpMatrix(elemat2);
      return 0;
    }
    case DRT::ELEMENTS::struct_calc_update_istep:
    {
      std::visit([&](auto& interface)
          { interface->Update(*this, *SolidMaterial(), discretization, lm, params); },
          solid_calc_variant_);

      return 0;
    }
    case DRT::ELEMENTS::struct_update_prestress:
    {
      UpdatePrestress(solid_calc_variant_, *this, *SolidMaterial(), discretization, lm, params);
      return 0;
    }
    case DRT::ELEMENTS::struct_calc_recover:
    {
      std::visit([&](auto& interface) { interface->Recover(*this, discretization, lm, params); },
          solid_calc_variant_);

      return 0;
    }
    case struct_calc_stress:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->CalculateStress(*this, *SolidMaterial(),
                StressIO{GetIOStressType(*this, params), GetStressData(*this, params)},
                StrainIO{GetIOStrainType(*this, params), GetStrainData(*this, params)},
                discretization, lm, params);
          },
          solid_calc_variant_);

      return 0;
    }
    case DRT::ELEMENTS::struct_calc_energy:
    {
      double int_energy = std::visit(
          [&](auto& interface) {
            return interface->CalculateInternalEnergy(
                *this, *SolidMaterial(), discretization, lm, params);
          },
          solid_calc_variant_);


      if (IsParamsInterface())
      {
        // new structural time integration
        ParamsInterface().AddContributionToEnergyType(int_energy, STR::internal_energy);
      }
      else
      {
        // old structural time integration
        // check length of elevec1
        if (elevec1.length() < 1) dserror("The given result vector is too short.");

        elevec1(0) = int_energy;
      }
      return 0;
    }
    case struct_init_gauss_point_data_output:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->InitializeGaussPointDataOutput(
                *this, *SolidMaterial(), *ParamsInterface().GaussPointDataOutputManagerPtr());
          },
          solid_calc_variant_);

      return 0;
    }
    case struct_gauss_point_data_output:
    {
      std::visit(
          [&](auto& interface)
          {
            interface->EvaluateGaussPointDataOutput(
                *this, *SolidMaterial(), *ParamsInterface().GaussPointDataOutputManagerPtr());
          },
          solid_calc_variant_);

      return 0;
    }
    case ELEMENTS::struct_calc_reset_istep:
    {
      std::visit([&](auto& interface) { interface->ResetToLastConverged(*this, *SolidMaterial()); },
          solid_calc_variant_);

      return 0;
    }
    case DRT::ELEMENTS::struct_calc_predict:
      // do nothing for now
      return 0;
    default:
      dserror("The element action %s is not yet implemented for the new solid elements",
          ActionType2String(action).c_str());
  }

  return 0;
}
int DRT::ELEMENTS::Solid::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)
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

BACI_NAMESPACE_CLOSE
