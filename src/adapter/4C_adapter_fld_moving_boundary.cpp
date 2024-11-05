// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_fld_fbi_movingboundary.hpp"
#include "4C_adapter_fld_fluid_ale.hpp"
#include "4C_adapter_fld_fluid_ale_xfem.hpp"
#include "4C_adapter_fld_fluid_immersed.hpp"
#include "4C_adapter_fld_fluid_xfem.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::FluidMovingBoundaryBaseAlgorithm::FluidMovingBoundaryBaseAlgorithm(
    const Teuchos::ParameterList& prbdyn, std::string condname)
{
  const Core::ProblemType probtyp = Global::Problem::instance()->get_problem_type();

  // switch between moving domain fluid implementations
  switch (probtyp)
  {
    case Core::ProblemType::fsi:
    case Core::ProblemType::fluid_ale:
    case Core::ProblemType::fsi_redmodels:
    {
      // std::cout << "using FluidAle as FluidMovingBoundary" << std::endl;
      fluid_ = std::make_shared<FluidAle>(prbdyn, condname);
      break;
    }
    case Core::ProblemType::fluid_xfem:
    case Core::ProblemType::fsi_xfem:
    {
      const Teuchos::ParameterList xfluid = Global::Problem::instance()->x_fluid_dynamic_params();
      bool alefluid = xfluid.sublist("GENERAL").get<bool>("ALE_XFluid");
      if (!alefluid)  // xfluid
      {
        // std::cout << "using FluidXFEM as FluidMovingBoundary" << endl;
        fluid_ = std::make_shared<FluidXFEM>(prbdyn, condname);
      }
      else  // xafluid
      {
        fluid_ = std::make_shared<FluidAleXFEM>(prbdyn, condname);
      }
      break;
    }
    case Core::ProblemType::immersed_fsi:
    {
      fluid_ = std::make_shared<FluidImmersed>(prbdyn, condname);
      break;
    }
    case Core::ProblemType::fbi:
    {
      fluid_ = std::make_shared<FBIFluidMB>(prbdyn, condname);
      break;
    }
    default:
      FOUR_C_THROW("fsi type not supported");
      break;
  }
}

FOUR_C_NAMESPACE_CLOSE
