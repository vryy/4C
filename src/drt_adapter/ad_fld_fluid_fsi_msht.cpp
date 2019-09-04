/*--------------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for fsi with internal mesh tying or mesh sliding

\maintainer Matthias Mayr

\level 3
*/
/*--------------------------------------------------------------------------*/

#include "ad_fld_fluid_fsi_msht.H"

#include "../drt_fluid/fluid_utils_mapextractor.H"

/*======================================================================*/
/* constructor */
ADAPTER::FluidFSIMsht::FluidFSIMsht(Teuchos::RCP<Fluid> fluid,
    Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<LINALG::Solver> solver,
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
