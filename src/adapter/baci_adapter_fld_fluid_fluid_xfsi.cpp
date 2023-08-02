/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for XFSI allowing multiple fluid discretizations .Can only be used in
conjunction with XFluidFluid!

\level 2


*/
/*----------------------------------------------------------------------*/

#include "baci_adapter_fld_fluid_fluid_xfsi.H"

#include "baci_fluid_xfluid_fluid.H"
#include "baci_lib_discret_xfem.H"
#include "baci_linalg_mapextractor.H"
#include "baci_linalg_utils_sparse_algebra_math.H"
#include "baci_xfem_condition_manager.H"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

#include <set>
#include <vector>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidFluidXFSI::FluidFluidXFSI(Teuchos::RCP<Fluid> fluid,  // the XFluid object
    const std::string coupling_name_xfsi, Teuchos::RCP<CORE::LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params, Teuchos::RCP<IO::DiscretizationWriter> output)
    : XFluidFSI(fluid, coupling_name_xfsi, solver, params, output)
{
  // make sure
  if (fluid_ == Teuchos::null) dserror("Failed to create the underlying fluid adapter");
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
  xfluidfluid_->UseBlockMatrix();
}
