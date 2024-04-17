/*--------------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for fsi with internal mesh tying or mesh sliding


\level 3
*/
/*--------------------------------------------------------------------------*/

#include "baci_adapter_fld_fluid_fsi_msht.hpp"

#include "baci_fluid_utils_mapextractor.hpp"

FOUR_C_NAMESPACE_OPEN

/*======================================================================*/
/* constructor */
ADAPTER::FluidFSIMsht::FluidFSIMsht(Teuchos::RCP<Fluid> fluid,
    Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<CORE::LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params, Teuchos::RCP<IO::DiscretizationWriter> output,
    bool isale, bool dirichletcond)
    : FluidFSI(fluid, dis, solver, params, output, isale, dirichletcond),
      fsiinterface_(Teuchos::rcp(new FLD::UTILS::FsiMapExtractor()))
{
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSIMsht::Init()
{
  // call base class init
  FluidFSI::Init();

  // create fluid map extractor
  SetupFsiInterface();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSIMsht::SetupFsiInterface() { fsiinterface_->Setup(*dis_); }

FOUR_C_NAMESPACE_CLOSE
