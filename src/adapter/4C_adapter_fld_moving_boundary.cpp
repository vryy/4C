/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for moving boundary problems


\level 2
*/
/*----------------------------------------------------------------------*/


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
    case Core::ProblemType::freesurf:
    case Core::ProblemType::fsi_redmodels:
    {
      // std::cout << "using FluidAle as FluidMovingBoundary" << std::endl;
      fluid_ = Teuchos::rcp(new FluidAle(prbdyn, condname));
      break;
    }
    case Core::ProblemType::fluid_xfem:
    case Core::ProblemType::fsi_xfem:
    {
      const Teuchos::ParameterList xfluid = Global::Problem::instance()->x_fluid_dynamic_params();
      bool alefluid = Core::UTILS::IntegralValue<bool>((xfluid.sublist("GENERAL")), "ALE_XFluid");
      if (!alefluid)  // xfluid
      {
        // std::cout << "using FluidXFEM as FluidMovingBoundary" << endl;
        fluid_ = Teuchos::rcp(new FluidXFEM(prbdyn, condname));
      }
      else  // xafluid
      {
        fluid_ = Teuchos::rcp(new FluidAleXFEM(prbdyn, condname));
      }
      break;
    }
    case Core::ProblemType::immersed_fsi:
    {
      fluid_ = Teuchos::rcp(new FluidImmersed(prbdyn, condname));
      break;
    }
    case Core::ProblemType::fbi:
    {
      fluid_ = Teuchos::rcp(new FBIFluidMB(prbdyn, condname));
      break;
    }
    default:
      FOUR_C_THROW("fsi type not supported");
      break;
  }
}

FOUR_C_NAMESPACE_CLOSE
