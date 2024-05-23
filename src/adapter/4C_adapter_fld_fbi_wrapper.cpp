/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for fluid beam interaction

\level 2

*/
/*----------------------------------------------------------------------*/
#include "4C_adapter_fld_fbi_wrapper.hpp"

#include "4C_fluid_implicit_integration.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_sparseoperator.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*======================================================================*/
/* constructor */
ADAPTER::FluidFBI::FluidFBI(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<IO::DiscretizationWriter> output, bool isale, bool dirichletcond)
    : FluidFSI(fluid, dis, solver, params, output, isale, dirichletcond)
{
  // make sure
  if (Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_, true) == Teuchos::null)
    FOUR_C_THROW("Failed to create the correct underlying fluid adapter");
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void ADAPTER::FluidFBI::set_coupling_contributions(
    Teuchos::RCP<const CORE::LINALG::SparseOperator> matrix)
{
  Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_, true)
      ->set_coupling_contributions(matrix);
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

FOUR_C_NAMESPACE_CLOSE
