/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_fsi.cpp

\brief Fluid field adapter for fsi

Can only be used in conjunction with XFluid!

<pre>
Maintainer:  Benedikt Schott
             schott@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#include "ad_fld_xfluid_fsi.H"

#include "../drt_adapter/ad_fld_fluid.H"
#include "../drt_fluid/xfluid.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_utils.H"

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <vector>
#include <set>
/*======================================================================*/
/* constructor */
ADAPTER::XFluidFSI::XFluidFSI(Teuchos::RCP< Fluid> fluid,
    Teuchos::RCP<DRT::Discretization> xfluiddis,
    Teuchos::RCP<DRT::Discretization> soliddis,
    Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<IO::DiscretizationWriter> output )
: FluidWrapper(fluid),
  xfluiddis_(xfluiddis),
  soliddis_(soliddis),
  solver_(solver),
  params_(params)
{
  // make sure
  if (fluid_ == Teuchos::null)
    dserror("Failed to create the underlying fluid adapter");

  // cast fluid to fluidimplicit
  xfluid_ = Teuchos::rcp_dynamic_cast<FLD::XFluid>(fluid_);
  if (xfluid_ == Teuchos::null)
    dserror("Failed to cast ADAPTER::Fluid to FLD::FluidImplicitTimeInt.");

  interface_ = Teuchos::rcp(new FLD::UTILS::MapExtractor());
  meshmap_   = Teuchos::rcp(new LINALG::MapExtractor());

  //Assign after calling the Xfluid-constructor
  boundarydis_ = xfluid_->Boundary_Dis();

  // the solid mesh has to match the interface mesh
  // so we have to compute a interface true residual vector itrueresidual_
  interface_->Setup(*boundarydis_);
  xfluid_->SetSurfaceSplitter(&(*interface_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidFSI::ExtractInterfaceForces()
{
  cout << "ExtractInterfaceForces (itrueresnp)" << endl;

  // the trueresidual vector has to match the solid dis
  // it contains the forces acting on the structural surface
  return interface_->ExtractFSICondVector(xfluid_->ITrueResidual());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidFSI::ExtractInterfaceVeln()
{
  cout << "call ExtractInterfaceVeln() "<< endl;

  // it depends, when this method is called, and when velnp is updated
  // the FSI algorithm expects first an time update and then asks for the old time step velocity
  // meaning that it gets the velocity from the new time step
  // not clear? exactly! thats why the FSI time update should be more clear about it
  // needs discussion with the FSI people
  return interface_->ExtractFSICondVector(xfluid_->IVelnp());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidFSI::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  cout << "ApplyInterfaceVelocities" << endl;

  interface_->InsertFSICondVector(ivel,xfluid_->IVelnp());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidFSI::ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> idisp)
{
  cout << "ApplyMeshDisplacement" << endl;

  interface_->InsertFSICondVector(idisp,xfluid_->IDispnp());

}
