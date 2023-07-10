/*----------------------------------------------------------------------*/
/*! \file

\brief A C++ wrapper for the solid element

This file contains the element-specific evaluation routines such as
Evaluate(...), EvaluateNeumann(...), etc.

\level 1
*/

#include "lib_elements_paramsinterface.H"
#include "structure_new_elements_paramsinterface.H"
#include "lib_elements_paramsinterface.H"
#include "solid_ele_neumann_evaluator.H"
#include "structure_new_elements_paramsinterface.H"
#include "utils_exceptions.H"
#include "solid_ele.H"
#include "solid_ele_factory.H"
#include "solid_ele_calc_interface.H"
#include "solid_ele_calc_lib.H"


int DRT::ELEMENTS::Solid::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  if (!material_post_setup_)
  {
    DRT::ELEMENTS::SolidFactory::ProvideImpl(
        this, this->GetEleTech(), this->GetKinemType(), this->GetEAStype())
        ->MaterialPostSetup(*this, *SolidMaterial());
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
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, this->GetEleTech(), this->GetKinemType(), this->GetEAStype())
          ->EvaluateNonlinearForceStiffnessMass(
              *this, *SolidMaterial(), discretization, lm, params, &elevec1, &elemat1, nullptr);
      return 0;
    }
    case struct_calc_internalforce:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, this->GetEleTech(), this->GetKinemType(), this->GetEAStype())
          ->EvaluateNonlinearForceStiffnessMass(
              *this, *SolidMaterial(), discretization, lm, params, &elevec1, nullptr, nullptr);
      return 0;
    }
    case struct_calc_nlnstiffmass:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, this->GetEleTech(), this->GetKinemType(), this->GetEAStype())
          ->EvaluateNonlinearForceStiffnessMass(
              *this, *SolidMaterial(), discretization, lm, params, &elevec1, &elemat1, &elemat2);
      return 0;
    }
    case struct_calc_nlnstifflmass:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, this->GetEleTech(), this->GetKinemType(), this->GetEAStype())
          ->EvaluateNonlinearForceStiffnessMass(
              *this, *SolidMaterial(), discretization, lm, params, &elevec1, &elemat1, &elemat2);
      DRT::ELEMENTS::LumpMatrix(elemat2);
      return 0;
    }
    case DRT::ELEMENTS::struct_calc_update_istep:
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, this->GetEleTech(), this->GetKinemType(), this->GetEAStype())
          ->Update(*this, *SolidMaterial(), discretization, lm, params);
      return 0;
    case DRT::ELEMENTS::struct_calc_recover:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, this->GetEleTech(), this->GetKinemType(), this->GetEAStype())
          ->Recover(*this, discretization, lm, params);
      return 0;
    }
    case struct_calc_stress:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, this->GetEleTech(), this->GetKinemType(), this->GetEAStype())
          ->CalculateStress(*this, *SolidMaterial(),
              StressIO{GetIOStressType(*this, params), GetMutableStressData(*this, params)},
              StrainIO{GetIOStrainType(*this, params), GetMutableStrainData(*this, params)},
              discretization, lm, params);
      return 0;
    }
    case DRT::ELEMENTS::struct_calc_energy:
    {
      double int_energy =
          DRT::ELEMENTS::SolidFactory::ProvideImpl(
              this, this->GetEleTech(), this->GetKinemType(), this->GetEAStype())
              ->CalculateInternalEnergy(*this, *SolidMaterial(), discretization, lm, params);

      if (IsParamsInterface())
      {
        // new structural time integration
        ParamsInterface().AddContributionToEnergyType(int_energy, STR::internal_energy);
      }
      else
      {
        // old structural time integration
        // check length of elevec1
        if (elevec1.Length() < 1) dserror("The given result vector is too short.");

        elevec1(0) = int_energy;
      }
      return 0;
    }
    case struct_init_gauss_point_data_output:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, this->GetEleTech(), this->GetKinemType(), this->GetEAStype())
          ->InitializeGaussPointDataOutput(
              *this, *SolidMaterial(), *ParamsInterface().MutableGaussPointDataOutputManagerPtr());
      return 0;
    }
    case struct_gauss_point_data_output:
    {
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, this->GetEleTech(), this->GetKinemType(), this->GetEAStype())
          ->EvaluateGaussPointDataOutput(
              *this, *SolidMaterial(), *ParamsInterface().MutableGaussPointDataOutputManagerPtr());
      return 0;
    }
    case ELEMENTS::struct_calc_reset_all:
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, this->GetEleTech(), this->GetKinemType(), this->GetEAStype())
          ->ResetAll(*this, *SolidMaterial());
      return 0;
    case ELEMENTS::struct_calc_reset_istep:
      DRT::ELEMENTS::SolidFactory::ProvideImpl(
          this, this->GetEleTech(), this->GetKinemType(), this->GetEAStype())
          ->ResetToLastConverged(*this, *SolidMaterial());
      return 0;
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
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
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
