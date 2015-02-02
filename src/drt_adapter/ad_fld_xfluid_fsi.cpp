/*----------------------------------------------------------------------*/
/*!
\file ad_fld_xfluid_fsi.cpp

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
#include "../drt_fluid_xfluid/xfluid.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_lib/drt_discret_xfem.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_utils.H"

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <vector>
#include <set>
/*======================================================================*/
/* constructor */
ADAPTER::XFluidFSI::XFluidFSI(Teuchos::RCP< Fluid> fluid,   // the XFluid object
    Teuchos::RCP<DRT::Discretization> xfluiddis,
    Teuchos::RCP<DRT::Discretization> soliddis,
    Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<IO::DiscretizationWriter> output )
: FluidWrapper(fluid),                                      // the XFluid object is set as fluid_ in the FluidWrapper
  fpsiinterface_(Teuchos::rcp(new FLD::UTILS::MapExtractor())),
  xfluiddis_(xfluiddis),
  soliddis_(soliddis),
  solver_(solver),
  params_(params)
{
  // make sure
  if (fluid_ == Teuchos::null)
    dserror("Failed to create the underlying fluid adapter");
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidFSI::Init()
{
  // call base class init
  FluidWrapper::Init();

  // cast fluid to fluidimplicit
  xfluid_ = Teuchos::rcp_dynamic_cast<FLD::XFluid>(fluid_);
  if (xfluid_ == Teuchos::null)
    dserror("Failed to cast ADAPTER::Fluid to FLD::XFluid.");

  structinterface_ = Teuchos::rcp(new FLD::UTILS::MapExtractor());
  interface_ = Teuchos::rcp(new FLD::UTILS::MapExtractor());
  meshmap_   = Teuchos::rcp(new LINALG::MapExtractor());


  // the solid mesh has to match the interface mesh
  // so we have to compute a interface true residual vector itrueresidual_
  structinterface_->Setup(*xfluid_->BoundaryDiscretization());

  interface_->Setup(*xfluiddis_,false, true); //Always Create overlapping FSI/FPSI Interface

  fpsiinterface_->Setup(*xfluiddis_,true, true); //Always Create overlapping FSI/FPSI Interface
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::XFluidFSI::TimeScaling() const
{
  // second order (OST(0.5) except for the first starting step, otherwise 1st order BackwardEuler
  if (params_->get<bool>("interface second order"))
    return 2./xfluid_->Dt();
  else
    return 1./xfluid_->Dt();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidFSI::ExtractStructInterfaceForces()
{
  // the trueresidual vector has to match the solid dis
  // it contains the forces acting on the structural surface
  return structinterface_->ExtractFSICondVector(xfluid_->ITrueResidual());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidFSI::ExtractStructInterfaceVeln()
{
  // it depends, when this method is called, and when velnp is updated
  // the FSI algorithm expects first an time update and then asks for the old time step velocity
  // meaning that it gets the velocity from the new time step
  // not clear? exactly! thats why the FSI time update should be more clear about it
  // needs discussion with the FSI people
  return structinterface_->ExtractFSICondVector(xfluid_->IVeln());
}


/*----------------------------------------------------------------------*/
// apply the interface velocities to the fluid
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidFSI::ApplyStructInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  structinterface_->InsertFSICondVector(ivel,xfluid_->IVelnp());
}


/*----------------------------------------------------------------------*/
//  apply the interface displacements to the fluid
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidFSI::ApplyStructMeshDisplacement(Teuchos::RCP<const Epetra_Vector> idisp)
{
  structinterface_->InsertFSICondVector(idisp,xfluid_->IDispnp());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::XFluidFSI::SetMeshMap(Teuchos::RCP<const Epetra_Map> mm)
{
  meshmap_->Setup(*xfluid_->DiscretisationXFEM()->InitialDofRowMap(),mm,LINALG::SplitMap(*xfluid_->DiscretisationXFEM()->InitialDofRowMap(),*mm));
}

/*----------------------------------------------------------------------*/
//  apply the ale displacements to the fluid
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidFSI::ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp)
{
  meshmap_->InsertCondVector(fluiddisp,xfluid_->WriteAccessDispnp());

  // new grid velocity
  xfluid_->UpdateGridv();
}


/*----------------------------------------------------------------------*
 * convert increment of displacement to increment in velocity
 * Delta d = d^(n+1,i+1)-d^n is converted to the interface velocity increment
 * Delta u = u^(n+1,i+1)-u^n
 * via first order or second order OST-discretization of d/dt d(t) = u(t)
 *----------------------------------------------------------------------*/
void ADAPTER::XFluidFSI::DisplacementToVelocity(
    Teuchos::RCP<Epetra_Vector> fcx         /// Delta d = d^(n+1,i+1)-d^n
)
{

  // get interface velocity at t(n)
  const Teuchos::RCP<const Epetra_Vector> veln = structinterface_->ExtractFSICondVector(xfluid_->IVeln());

#ifdef DEBUG
  // check, whether maps are the same
  if (! fcx->Map().PointSameAs(veln->Map()))  { dserror("Maps do not match, but they have to."); }
#endif

  /*
   * Delta u(n+1,i+1) = fac * (Delta d(n+1,i+1) - dt * u(n))
   *
   *             / = 2 / dt   if interface time integration is second order
   * with fac = |
   *             \ = 1 / dt   if interface time integration is first order
   */
  const double timescale = TimeScaling();
  fcx->Update(-timescale*xfluid_->Dt(),*veln,timescale);
}


/// return xfluid coupling matrix between structure and fluid as sparse matrices
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::XFluidFSI::C_Struct_Fluid_Matrix()
{
  return xfluid_->C_Struct_Fluid_Matrix();
}

/// return xfluid coupling matrix between fluid and structure as sparse matrices
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::XFluidFSI::C_Fluid_Struct_Matrix()
{
  return xfluid_->C_Fluid_Struct_Matrix();
}

/// return xfluid coupling matrix between structure and structure as sparse matrices
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::XFluidFSI::C_Struct_Struct_Matrix()
{
  return xfluid_->C_Struct_Struct_Matrix();
}

/// return xfluid coupling matrix between structure and structure as sparse matrices
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidFSI::RHS_Struct_Vec()
{
  return xfluid_->RHS_Struct_Vec();
}

/*----------------------------------------------------------------------*
 * Rebuild FSI interface in case of crack-FSI problem               sudhakar 03/14
 * This is needed when we add new nodes to the FSI interface
 *----------------------------------------------------------------------*/
void ADAPTER::XFluidFSI::RebuildFSIStructInterface()
{
  StructInterface()->Setup(*xfluid_->BoundaryDiscretization());
}

/// GmshOutput for background mesh and cut mesh
void ADAPTER::XFluidFSI::GmshOutput(
    const std::string & name,            ///< name for output file
    const int step,                      ///< step number
    const int count,                     ///< counter for iterations within a global time step
    Teuchos::RCP<Epetra_Vector> vel,     ///< vector holding velocity and pressure dofs
    Teuchos::RCP<Epetra_Vector> acc      ///< vector holding accelerations
)
{
  xfluid_->GmshOutput(name, step, count, vel, acc);
}
