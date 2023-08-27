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
#include "baci_solid_ele_factory.H"
#include "baci_solid_ele_neumann_evaluator.H"
#include "baci_structure_new_elements_paramsinterface.H"
#include "baci_utils_exceptions.H"


int DRT::ELEMENTS::Solid::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  if (!material_post_setup_)
  {
    solid_interface_->MaterialPostSetup(*this, *SolidMaterial());
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
      solid_interface_->EvaluateNonlinearForceStiffnessMass(
          *this, *SolidMaterial(), discretization, lm, params, &elevec1, &elemat1, nullptr);
      return 0;
    }
    case DRT::ELEMENTS::struct_calc_nlnstiff_gemm:
    {
      solid_interface_->EvaluateNonlinearForceStiffnessMassGEMM(
          *this, *SolidMaterial(), discretization, lm, params, &elevec1, &elemat1, nullptr);
      return 0;
    }
    case struct_calc_internalforce:
    {
      solid_interface_->EvaluateNonlinearForceStiffnessMass(
          *this, *SolidMaterial(), discretization, lm, params, &elevec1, nullptr, nullptr);
      return 0;
    }
    case struct_calc_nlnstiffmass:
    {
      solid_interface_->EvaluateNonlinearForceStiffnessMass(
          *this, *SolidMaterial(), discretization, lm, params, &elevec1, &elemat1, &elemat2);
      return 0;
    }
    case struct_calc_nlnstifflmass:
    {
      solid_interface_->EvaluateNonlinearForceStiffnessMass(
          *this, *SolidMaterial(), discretization, lm, params, &elevec1, &elemat1, &elemat2);
      LumpMatrix(elemat2);
      return 0;
    }
    case DRT::ELEMENTS::struct_calc_update_istep:
    {
      solid_interface_->Update(*this, *SolidMaterial(), discretization, lm, params);
      return 0;
    }
    case DRT::ELEMENTS::struct_calc_recover:
    {
      solid_interface_->Recover(*this, discretization, lm, params);
      return 0;
    }
    case struct_calc_stress:
    {
      solid_interface_->CalculateStress(*this, *SolidMaterial(),
          StressIO{GetIOStressType(*this, params), GetStressData(*this, params)},
          StrainIO{GetIOStrainType(*this, params), GetStrainData(*this, params)}, discretization,
          lm, params);
      return 0;
    }
    case DRT::ELEMENTS::struct_calc_energy:
    {
      double int_energy = solid_interface_->CalculateInternalEnergy(
          *this, *SolidMaterial(), discretization, lm, params);

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
      solid_interface_->InitializeGaussPointDataOutput(
          *this, *SolidMaterial(), *ParamsInterface().GaussPointDataOutputManagerPtr());
      return 0;
    }
    case struct_gauss_point_data_output:
    {
      solid_interface_->EvaluateGaussPointDataOutput(
          *this, *SolidMaterial(), *ParamsInterface().GaussPointDataOutputManagerPtr());
      return 0;
    }
    case ELEMENTS::struct_calc_reset_all:
    {
      solid_interface_->ResetAll(*this, *SolidMaterial());
      return 0;
    }
    case ELEMENTS::struct_calc_reset_istep:
    {
      solid_interface_->ResetToLastConverged(*this, *SolidMaterial());
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
