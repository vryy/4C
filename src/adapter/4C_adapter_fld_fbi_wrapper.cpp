// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_fld_fbi_wrapper.hpp"

#include "4C_fluid_implicit_integration.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_sparseoperator.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

/*======================================================================*/
/* constructor */
Adapter::FluidFBI::FluidFBI(std::shared_ptr<Fluid> fluid,
    std::shared_ptr<Core::FE::Discretization> dis, std::shared_ptr<Core::LinAlg::Solver> solver,
    std::shared_ptr<Teuchos::ParameterList> params,
    std::shared_ptr<Core::IO::DiscretizationWriter> output, bool isale, bool dirichletcond)
    : FluidFSI(fluid, dis, solver, params, output, isale, dirichletcond)
{
  // make sure
  if (std::dynamic_pointer_cast<FLD::FluidImplicitTimeInt>(fluid_) == nullptr)
    FOUR_C_THROW("Failed to create the correct underlying fluid adapter");
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void Adapter::FluidFBI::set_coupling_contributions(
    std::shared_ptr<const Core::LinAlg::SparseOperator> matrix)
{
  std::dynamic_pointer_cast<FLD::FluidImplicitTimeInt>(fluid_)->set_coupling_contributions(matrix);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void Adapter::FluidFBI::reset_external_forces()
{
  std::dynamic_pointer_cast<FLD::FluidImplicitTimeInt>(fluid_)->reset_external_forces();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

std::shared_ptr<const FLD::Meshtying> Adapter::FluidFBI::get_meshtying()
{
  return std::dynamic_pointer_cast<FLD::FluidImplicitTimeInt>(fluid_)->get_meshtying();
}

FOUR_C_NAMESPACE_CLOSE
