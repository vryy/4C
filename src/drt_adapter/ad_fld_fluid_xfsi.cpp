/*----------------------------------------------------------------------*/
/*!

\brief Fluid field adapter for XFSI. Can only be used in conjunction with XFluid!

\level 1

\maintainer  Christoph Ager

*/
/*----------------------------------------------------------------------*/

#include "ad_fld_fluid_xfsi.H"

#include "../drt_adapter/ad_fld_fluid.H"
#include "../drt_fluid_xfluid/xfluid.H"
#include "../drt_fluid_xfluid/xfluidfluid.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_lib/drt_discret_xfem.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_utils.H"
#include "../drt_xfem/xfem_condition_manager.H"

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <vector>
#include <set>
/*======================================================================*/
/* constructor */
ADAPTER::XFluidFSI::XFluidFSI(Teuchos::RCP<Fluid> fluid,  // the XFluid object
    const std::string coupling_name,                      // name of the FSI coupling condition
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<IO::DiscretizationWriter> output)
    : FluidWrapper(fluid),  // the XFluid object is set as fluid_ in the FluidWrapper
      fpsiinterface_(Teuchos::rcp(new FLD::UTILS::MapExtractor())),
      coupling_name_(coupling_name),
      solver_(solver),
      params_(params)
{
  // make sure
  if (fluid_ == Teuchos::null) dserror("Failed to create the underlying fluid adapter");
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
  if (xfluid_ == Teuchos::null) dserror("Failed to cast ADAPTER::Fluid to FLD::XFluid.");

  // NOTE: currently we are using the XFluidFSI adapter also for pure ALE-fluid problems with
  // level-set boundary in this case no mesh coupling object is available and no interface objects
  // can be created
  Teuchos::RCP<XFEM::MeshCoupling> mc = xfluid_->GetMeshCoupling(coupling_name_);

  if (mc != Teuchos::null)  // classical mesh coupling case for FSI
  {
    // get the mesh coupling object
    mesh_coupling_fsi_ = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFSI>(mc, true);

    structinterface_ = Teuchos::rcp(new FLD::UTILS::MapExtractor());

    // the solid mesh has to match the interface mesh
    // so we have to compute a interface true residual vector itrueresidual_
    structinterface_->Setup(*mesh_coupling_fsi_->GetCutterDis());
  }

  interface_ = Teuchos::rcp(new FLD::UTILS::MapExtractor());

  interface_->Setup(
      *xfluid_->Discretization(), false, true);  // Always Create overlapping FSI/FPSI Interface

  fpsiinterface_->Setup(
      *xfluid_->Discretization(), true, true);  // Always Create overlapping FSI/FPSI Interface

  meshmap_ = Teuchos::rcp(new LINALG::MapExtractor());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::XFluidFSI::TimeScaling() const
{
  // second order (OST(0.5) except for the first starting step, otherwise 1st order BackwardEuler
  if (params_->get<bool>("interface second order"))
    return 2. / xfluid_->Dt();
  else
    return 1. / xfluid_->Dt();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const Teuchos::RCP<DRT::Discretization> ADAPTER::XFluidFSI::BoundaryDiscretization()
{
  return mesh_coupling_fsi_->GetCutterDis();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidFSI::ExtractStructInterfaceForces()
{
  // the trueresidual vector has to match the solid dis
  // it contains the forces acting on the structural surface
  return structinterface_->ExtractFSICondVector(mesh_coupling_fsi_->ITrueResidual());
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
  return structinterface_->ExtractFSICondVector(mesh_coupling_fsi_->IVeln());
}


/*----------------------------------------------------------------------*/
// apply the interface velocities to the fluid
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidFSI::ApplyStructInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  structinterface_->InsertFSICondVector(ivel, mesh_coupling_fsi_->IVelnp());
}


/*----------------------------------------------------------------------*/
//  apply the interface displacements to the fluid
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidFSI::ApplyStructMeshDisplacement(Teuchos::RCP<const Epetra_Vector> idisp)
{
  // update last increment, before we set new idispnp
  mesh_coupling_fsi_->UpdateDisplacementIterationVectors();

  // set new idispnp
  structinterface_->InsertFSICondVector(idisp, mesh_coupling_fsi_->IDispnp());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::XFluidFSI::SetMeshMap(Teuchos::RCP<const Epetra_Map> mm)
{
  meshmap_->Setup(*xfluid_->DiscretisationXFEM()->InitialDofRowMap(), mm,
      LINALG::SplitMap(*xfluid_->DiscretisationXFEM()->InitialDofRowMap(), *mm));
}

/*----------------------------------------------------------------------*/
//  apply the ale displacements to the fluid
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidFSI::ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp)
{
  meshmap_->InsertCondVector(fluiddisp, xfluid_->WriteAccessDispnp());

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
    Teuchos::RCP<Epetra_Vector> fcx  /// Delta d = d^(n+1,i+1)-d^n
)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<const Epetra_Vector> veln =
      structinterface_->ExtractFSICondVector(mesh_coupling_fsi_->IVeln());

#ifdef DEBUG
  // check, whether maps are the same
  if (!fcx->Map().PointSameAs(veln->Map()))
  {
    dserror("Maps do not match, but they have to.");
  }
#endif

  /*
   * Delta u(n+1,i+1) = fac * (Delta d(n+1,i+1) - dt * u(n))
   *
   *             / = 2 / dt   if interface time integration is second order
   * with fac = |
   *             \ = 1 / dt   if interface time integration is first order
   */
  const double timescale = TimeScaling();
  fcx->Update(-timescale * xfluid_->Dt(), *veln, timescale);
}


/// return xfluid coupling matrix between structure and fluid as sparse matrices
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::XFluidFSI::C_Struct_Fluid_Matrix()
{
  return xfluid_->C_sx_Matrix(coupling_name_);
}

/// return xfluid coupling matrix between fluid and structure as sparse matrices
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::XFluidFSI::C_Fluid_Struct_Matrix()
{
  return xfluid_->C_xs_Matrix(coupling_name_);
}

/// return xfluid coupling matrix between structure and structure as sparse matrices
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::XFluidFSI::C_Struct_Struct_Matrix()
{
  return xfluid_->C_ss_Matrix(coupling_name_);
}

/// return xfluid coupling matrix between structure and structure as sparse matrices
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidFSI::RHS_Struct_Vec()
{
  return xfluid_->RHS_s_Vec(coupling_name_);
}

/// GmshOutput for background mesh and cut mesh
void ADAPTER::XFluidFSI::GmshOutput(const std::string& name,  ///< name for output file
    const int step,                                           ///< step number
    const int count,                  ///< counter for iterations within a global time step
    Teuchos::RCP<Epetra_Vector> vel,  ///< vector holding velocity and pressure dofs
    Teuchos::RCP<Epetra_Vector> acc   ///< vector holding accelerations
)
{
  // TODO (kruse): find a substitute!
  // xfluid_->GmshOutput(name, step, count, vel, acc);
  dserror("Gmsh output for XFSI during Newton currently not available.");
}
