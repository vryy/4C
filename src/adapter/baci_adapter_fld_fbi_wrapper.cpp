/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for fluid beam interaction

\level 2

*/
/*----------------------------------------------------------------------*/
#include "baci_adapter_fld_fbi_wrapper.H"
#include "baci_fluid_implicit_integration.H"
#include "baci_io_control.H"
#include "baci_linalg_sparseoperator.H"
#include <Teuchos_RCP.hpp>

/*======================================================================*/
/* constructor */
ADAPTER::FluidFBI::FluidFBI(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
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

void ADAPTER::FluidFBI::SetCouplingContributions(
    Teuchos::RCP<const CORE::LINALG::SparseOperator> matrix)
{
  Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_, true)
      ->SetCouplingContributions(matrix);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void ADAPTER::FluidFBI::ResetExternalForces()
{
  Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_, true)->ResetExternalForces();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

Teuchos::RCP<const FLD::Meshtying> ADAPTER::FluidFBI::GetMeshtying()
{
  return Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_, true)->GetMeshtying();
}
