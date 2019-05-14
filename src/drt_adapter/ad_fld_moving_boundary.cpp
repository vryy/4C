/*----------------------------------------------------------------------*/
/*!

\brief Fluid field adapter for moving boundary problems

\maintainer  Christoph Ager

\level 2
*/
/*----------------------------------------------------------------------*/


#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "ad_fld_fluid_ale.H"
#include "ad_fld_fluid_ale_immersed.H"
#include "ad_fld_fluid_xfem.H"
#include "ad_fld_fluid_immersed.H"
#include "ad_fld_fluid_ale_immersed.H"
#include "ad_fld_fluid_ale_xfem.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidMovingBoundaryBaseAlgorithm::FluidMovingBoundaryBaseAlgorithm(
    const Teuchos::ParameterList& prbdyn, std::string condname)
{
  const PROBLEM_TYP probtyp = DRT::Problem::Instance()->ProblemType();

  // switch between moving domain fluid implementations
  switch (probtyp)
  {
    case prb_fsi:
    case prb_fluid_ale:
    case prb_freesurf:
    case prb_fsi_redmodels:
    {
      // std::cout << "using FluidAle as FluidMovingBoundary" << std::endl;
      fluid_ = Teuchos::rcp(new FluidAle(prbdyn, condname));
      break;
    }
    case prb_fluid_xfem:
    case prb_fsi_xfem:
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
    case prb_immersed_fsi:
    case prb_immersed_membrane_fsi:
    {
      fluid_ = Teuchos::rcp(new FluidImmersed(prbdyn, condname));
      break;
    }
    case prb_immersed_ale_fsi:
    {
      fluid_ = Teuchos::rcp(new FluidAleImmersed(prbdyn, condname));
      break;
    }
    default:
      dserror("fsi type not supported");
      break;
  }
}
