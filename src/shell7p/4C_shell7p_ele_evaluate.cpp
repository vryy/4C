/*! \file
 *
\brief Nonlinear Shell 7-Parameter Model Finite Element evaluation

\level 3
*/

#include "4C_global_data.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_shell7p_ele.hpp"
#include "4C_shell7p_ele_calc_interface.hpp"
#include "4C_shell7p_ele_neumann_evaluator.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_function_of_time.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace
{
  inline std::vector<char>& GetMutableStressData(
      const DRT::ELEMENTS::Shell7p& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return *ele.StrParamsInterface().StressDataPtr();
    }
    else
    {
      return *params.get<Teuchos::RCP<std::vector<char>>>("stress");
    }
  }

  inline std::vector<char>& GetMutableStrainData(
      const DRT::ELEMENTS::Shell7p& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return *ele.StrParamsInterface().StrainDataPtr();
    }
    else
    {
      return *params.get<Teuchos::RCP<std::vector<char>>>("strain");
    }
  }

  inline INPAR::STR::StressType GetIOStressType(
      const DRT::ELEMENTS::Shell7p& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return ele.StrParamsInterface().GetStressOutputType();
    }
    else
    {
      return CORE::UTILS::GetAsEnum<INPAR::STR::StressType>(params, "iostress");
    }
  }

  inline INPAR::STR::StrainType GetIOStrainType(
      const DRT::ELEMENTS::Shell7p& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return ele.StrParamsInterface().GetStrainOutputType();
    }
    else
    {
      return CORE::UTILS::GetAsEnum<INPAR::STR::StrainType>(params, "iostrain");
    }
  }
}  // namespace

int DRT::ELEMENTS::Shell7p::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& dof_index_array,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  // get params interface pointer
  SetParamsInterfacePtr(params);

  const ELEMENTS::ActionType action = std::invoke(
      [&]()
      {
        if (IsParamsInterface())
          return StrParamsInterface().GetActionType();
        else
          return String2ActionType(params.get<std::string>("action", "none"));
      });

  // what should the element do
  switch (action)
  {
    // nonlinear stiffness and internal force vector
    case ELEMENTS::struct_calc_nlnstiff:
    {
      shell_interface_->EvaluateNonlinearForceStiffnessMass(*this, *SolidMaterial(), discretization,
          nodal_directors_, dof_index_array, params, &elevec1, &elemat1, nullptr);
    }
    break;
    case ELEMENTS::struct_calc_linstiff:
    {
      shell_interface_->EvaluateNonlinearForceStiffnessMass(*this, *SolidMaterial(), discretization,
          nodal_directors_, dof_index_array, params, &elevec1, &elemat1, &elemat2);
    }
    break;
    case ELEMENTS::struct_calc_internalforce:
    {
      shell_interface_->EvaluateNonlinearForceStiffnessMass(*this, *SolidMaterial(), discretization,
          nodal_directors_, dof_index_array, params, &elevec1, nullptr, nullptr);
    }
    break;
    case ELEMENTS::struct_calc_linstiffmass:
    {
      FOUR_C_THROW("Case not yet implemented: struct_calc_linstiffmass");
    }
    case ELEMENTS::struct_calc_nlnstiffmass:   // do mass, stiffness and internal forces
    case ELEMENTS::struct_calc_nlnstifflmass:  // do lump mass, stiffness and internal forces
    {
      shell_interface_->EvaluateNonlinearForceStiffnessMass(*this, *SolidMaterial(), discretization,
          nodal_directors_, dof_index_array, params, &elevec1, &elemat1, &elemat2);
      if (action == ELEMENTS::struct_calc_nlnstifflmass) STR::UTILS::SHELL::LumpMassMatrix(elemat2);
    }
    break;
    case ELEMENTS::struct_calc_recover:
    {
      shell_interface_->Recover(
          *this, discretization, dof_index_array, params, StrParamsInterface());
    }
    break;
    case ELEMENTS::struct_calc_stress:
    {
      shell_interface_->CalculateStressesStrains(*this, *SolidMaterial(),
          ShellStressIO{GetIOStressType(*this, params), GetMutableStressData(*this, params)},
          ShellStrainIO{GetIOStrainType(*this, params), GetMutableStrainData(*this, params)},
          discretization, nodal_directors_, dof_index_array, params);
    }
    break;
    case DRT::ELEMENTS::struct_calc_energy:
    {
      double int_energy = shell_interface_->CalculateInternalEnergy(
          *this, *SolidMaterial(), discretization, nodal_directors_, dof_index_array, params);

      if (IsParamsInterface())
      {
        // new structural time integration
        StrParamsInterface().AddContributionToEnergyType(int_energy, STR::internal_energy);
      }
      else
      {
        // old structural time integration
        // check length of elevec1
        if (elevec1.length() < 1) FOUR_C_THROW("The given result vector is too short.");
        elevec1(0) = int_energy;
      }
    }
    break;
    case ELEMENTS::struct_calc_update_istep:
    {
      shell_interface_->Update(
          *this, *SolidMaterial(), discretization, nodal_directors_, dof_index_array, params);
    }
    break;
    case ELEMENTS::struct_calc_reset_istep:
    {
      // Reset of history (if needed)
      shell_interface_->ResetToLastConverged(*this, *SolidMaterial());
    }
    break;
    case ELEMENTS::struct_calc_predict:
    case ELEMENTS::struct_create_backup:
    case ELEMENTS::struct_recover_from_backup:
    {
      // do nothing for now
    }
    break;
    default:
      FOUR_C_THROW("The element action %s is not yet implemented for the Shell element yet",
          ActionType2String(action).c_str());
  }
  return 0;
}

// Integrate a Surface Neumann boundary condition
int DRT::ELEMENTS::Shell7p::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition,
    std::vector<int>& dof_index_array, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseMatrix* elemat1)
{
  SetParamsInterfacePtr(params);

  const double time = std::invoke(
      [&]()
      {
        if (IsParamsInterface())
          return StrParamsInterface().GetTotalTime();
        else
          return params.get("total time", -1.0);
      });

  DRT::ELEMENTS::SHELL::EvaluateNeumannByElement(
      *this, discretization, condition, dof_index_array, elevec1, elemat1, time);
  return 0;
}

FOUR_C_NAMESPACE_CLOSE
