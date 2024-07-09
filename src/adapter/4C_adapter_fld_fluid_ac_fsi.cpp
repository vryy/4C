/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for AC-FS3I problems



\level 3
*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
#include "4C_adapter_fld_fluid_ac_fsi.hpp"

#include "4C_fluid_impedancecondition.hpp"
#include "4C_fluid_implicit_integration.hpp"

FOUR_C_NAMESPACE_OPEN

/*======================================================================*/
/* constructor */
Adapter::FluidACFSI::FluidACFSI(Teuchos::RCP<Fluid> fluid,
    Teuchos::RCP<Core::FE::Discretization> dis, Teuchos::RCP<Core::LinAlg::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Core::IO::DiscretizationWriter> output, bool isale, bool dirichletcond)
    : FluidFSI(fluid, dis, solver, params, output, isale, dirichletcond)
{
  return;
}

std::vector<double> Adapter::FluidACFSI::get_windkessel_errors()
{
  if (fluidimpl_->impedance_bc() == Teuchos::null)  // iff there is no windkessel condition
  {
    // FOUR_C_THROW("fluid field has no Windkessel!");
    std::vector<double> tmp(1, true);
    tmp[0] = 1000.0;
    return tmp;
  }

  return fluidimpl_->impedance_bc()->get_w_krelerrors();
}

FOUR_C_NAMESPACE_CLOSE
