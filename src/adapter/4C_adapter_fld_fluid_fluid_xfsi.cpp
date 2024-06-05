/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for XFSI allowing multiple fluid discretizations .Can only be used in
conjunction with XFluidFluid!

\level 2


*/
/*----------------------------------------------------------------------*/

#include "4C_adapter_fld_fluid_fluid_xfsi.hpp"

#include "4C_fluid_xfluid_fluid.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_xfem_condition_manager.hpp"
#include "4C_xfem_discretization.hpp"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <set>
#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidFluidXFSI::FluidFluidXFSI(Teuchos::RCP<Fluid> fluid,  // the XFluid object
    const std::string coupling_name_xfsi, Teuchos::RCP<CORE::LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<CORE::IO::DiscretizationWriter> output)
    : XFluidFSI(fluid, coupling_name_xfsi, solver, params, output)
{
  // make sure
  if (fluid_ == Teuchos::null) FOUR_C_THROW("Failed to create the underlying fluid adapter");
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFluidXFSI::Init()
{
  // call base class init
  XFluidFSI::Init();

  // cast fluid to fluidimplicit
  xfluidfluid_ = Teuchos::rcp_dynamic_cast<FLD::XFluidFluid>(xfluid_, true);

  // use block matrix for fluid-fluid, do nothing otherwise
  xfluidfluid_->use_block_matrix();
}

FOUR_C_NAMESPACE_CLOSE
