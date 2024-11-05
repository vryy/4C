// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_fld_fluid_fluid_xfsi.hpp"

#include "4C_fluid_xfluid_fluid.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_xfem_condition_manager.hpp"
#include "4C_xfem_discretization.hpp"

#include <Epetra_Map.h>

#include <memory>
#include <set>
#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::FluidFluidXFSI::FluidFluidXFSI(std::shared_ptr<Fluid> fluid,  // the XFluid object
    const std::string coupling_name_xfsi, std::shared_ptr<Core::LinAlg::Solver> solver,
    std::shared_ptr<Teuchos::ParameterList> params,
    std::shared_ptr<Core::IO::DiscretizationWriter> output)
    : XFluidFSI(fluid, coupling_name_xfsi, solver, params, output)
{
  // make sure
  if (fluid_ == nullptr) FOUR_C_THROW("Failed to create the underlying fluid adapter");
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidFluidXFSI::init()
{
  // call base class init
  XFluidFSI::init();

  // cast fluid to fluidimplicit
  xfluidfluid_ = std::dynamic_pointer_cast<FLD::XFluidFluid>(xfluid_);

  // use block matrix for fluid-fluid, do nothing otherwise
  xfluidfluid_->use_block_matrix();
}

FOUR_C_NAMESPACE_CLOSE
