/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for fluid beam interaction

\level 2

\maintainer Nora Hagmeyer
*/
/*----------------------------------------------------------------------*/
#include "ad_fld_fbi_wrapper.H"
#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_io/io_control.H"
#include "../linalg/linalg_utils.H"
#include <Teuchos_RCP.hpp>

/*======================================================================*/
/* constructor */
ADAPTER::FluidFBI::FluidFBI(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<IO::DiscretizationWriter> output, bool isale, bool dirichletcond)
    : FluidFSI(fluid, dis, solver, params, output, isale, dirichletcond)
{
  // make sure
  if (Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_, true) == Teuchos::null)
    dserror("Failed to create the correct underlying fluid adapter");
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void ADAPTER::FluidFBI::SetCouplingContributions(Teuchos::RCP<LINALG::SparseMatrix> matrix)
{
  Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_, true)
      ->SetCouplingContributions(matrix);
}
