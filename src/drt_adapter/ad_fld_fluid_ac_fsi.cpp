/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_ac_fsi.cpp

\brief Fluid field adapter for AC-FS3I problems

\date 2015-07-29

\maintainer Moritz Thon
            thon@mhpc.mw.tum.de
            089/289-10364

\level 3
*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
#include "ad_fld_fluid_ac_fsi.H"
#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid/fluidimpedancecondition.H"

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
