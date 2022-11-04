/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for moving boundary problems


\level 2
*/
/*----------------------------------------------------------------------*/


#include "drt_globalproblem.H"
#include "drt_validparameters.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "ad_fld_fbi_movingboundary.H"
#include "ad_fld_fluid_ale.H"
#include "ad_fld_fluid_xfem.H"
#include "ad_fld_fluid_immersed.H"
#include "ad_fld_fluid_ale_xfem.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidMovingBoundaryBaseAlgorithm::FluidMovingBoundaryBaseAlgorithm(
    const Teuchos::ParameterList& prbdyn, std::string condname)
{
  const ProblemType probtyp = DRT::Problem::Instance()->GetProblemType();

  // switch between moving domain fluid implementations
  switch (probtyp)
  {
    case ProblemType::fsi:
    case ProblemType::fluid_ale:
    case ProblemType::freesurf:
    case ProblemType::fsi_redmodels:
    {
      // std::cout << "using FluidAle as FluidMovingBoundary" << std::endl;
      fluid_ = Teuchos::rcp(new FluidAle(prbdyn, condname));
      break;
    }
    case ProblemType::fluid_xfem:
    case ProblemType::fsi_xfem:
    {
      const Teuchos::ParameterList xfluid = DRT::Problem::Instance()->XFluidDynamicParams();
      bool alefluid = DRT::INPUT::IntegralValue<bool>((xfluid.sublist("GENERAL")), "ALE_XFluid");
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
    case ProblemType::immersed_fsi:
    {
      fluid_ = Teuchos::rcp(new FluidImmersed(prbdyn, condname));
      break;
    }
    case ProblemType::fbi:
    {
      fluid_ = Teuchos::rcp(new FBIFluidMB(prbdyn, condname));
      break;
    }
    default:
      dserror("fsi type not supported");
      break;
  }
}
