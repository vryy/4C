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
Adapter::FluidFBI::FluidFBI(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<Core::FE::Discretization> dis,
    Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Core::IO::DiscretizationWriter> output, bool isale, bool dirichletcond)
    : FluidFSI(fluid, dis, solver, params, output, isale, dirichletcond)
{
  // make sure
  if (Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_, true) == Teuchos::null)
    FOUR_C_THROW("Failed to create the correct underlying fluid adapter");
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void Adapter::FluidFBI::set_coupling_contributions(
    Teuchos::RCP<const Core::LinAlg::SparseOperator> matrix)
{
  Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_, true)
      ->set_coupling_contributions(matrix);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void Adapter::FluidFBI::reset_external_forces()
{
  Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_, true)->reset_external_forces();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

Teuchos::RCP<const FLD::Meshtying> Adapter::FluidFBI::get_meshtying()
{
  return Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_, true)->get_meshtying();
}

FOUR_C_NAMESPACE_CLOSE
