/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for AC-FS3I problems



\level 3
*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
#include "adapter_fld_fluid_ac_fsi.H"
#include "fluid_implicit_integration.H"
#include "fluid_impedancecondition.H"

/*======================================================================*/
/* constructor */
ADAPTER::FluidACFSI::FluidACFSI(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<IO::DiscretizationWriter> output, bool isale, bool dirichletcond)
    : FluidFSI(fluid, dis, solver, params, output, isale, dirichletcond)
{
  return;
}

std::vector<double> ADAPTER::FluidACFSI::GetWindkesselErrors()
{
  if (fluidimpl_->ImpedanceBC_() == Teuchos::null)  // iff there is no windkessel condition
  {
    // dserror("fluid field has no Windkessel!");
    std::vector<double> tmp(1, true);
    tmp[0] = 1000.0;
    return tmp;
  }

  return fluidimpl_->ImpedanceBC_()->getWKrelerrors();
}
