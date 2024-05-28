/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for xfem-fluids with moving boundaries

\level 1


*/
/*----------------------------------------------------------------------*/
#include "4C_adapter_fld_fluid_xfem.hpp"

#include "4C_adapter_fld_fluid_xfsi.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidXFEM::FluidXFEM(const Teuchos::ParameterList& prbdyn, std::string condname)
    : fluid_(Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(prbdyn,
                              GLOBAL::Problem::Instance()->FluidDynamicParams(), "fluid", false))
                 ->fluid_field())
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidXFEM::discretization()
{
  // returns the boundary discretization
  // REMARK:
  // the returned discretization has to match the structure discretization at the interface coupling
  // (see FSI::Partitioned::Partitioned(const Epetra_Comm& comm) ) therefore return the boundary dis
  // this is similar to the matching of fluid dis and ale dis in case of ADAPTER::FluidALE
  return fluid_field()->discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidXFEM::boundary_discretization()
{
  // returns the boundary discretization
  // REMARK:
  // the returned discretization has to match the structure discretization at the interface coupling
  // (see FSI::Partitioned::Partitioned(const Epetra_Comm& comm) ) therefore return the boundary dis
  // this is similar to the matching of fluid dis and ale dis in case of ADAPTER::FluidALE

  Teuchos::RCP<XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<XFluidFSI>(fluid_field(), true);

  return xfluid->boundary_discretization();
}


/*----------------------------------------------------------------------*/
/// communication object at the struct interface
/*----------------------------------------------------------------------*/
Teuchos::RCP<FLD::UTILS::MapExtractor> const& ADAPTER::FluidXFEM::StructInterface()
{
  // returns the boundary discretization
  // REMARK:
  // the returned discretization has to match the structure discretization at the interface coupling
  // (see FSI::Partitioned::Partitioned(const Epetra_Comm& comm) ) therefore return the boundary dis
  // this is similar to the matching of fluid dis and ale dis in case of ADAPTER::FluidALE

  Teuchos::RCP<XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<XFluidFSI>(fluid_field(), true);

  return xfluid->StructInterface();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::prepare_time_step() { fluid_field()->prepare_time_step(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::Update() { fluid_field()->Update(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::Output() { fluid_field()->StatisticsAndOutput(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidXFEM::read_restart(int step)
{
  fluid_field()->read_restart(step);
  return fluid_field()->Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::nonlinear_solve(
    Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel)
{
  // if we have values at the interface we need to apply them

  // REMARK: for XFLUID idisp = Teuchos::null, ivel = Teuchos::null (called by fsi_fluid_xfem with
  // default Teuchos::null)
  //         for XFSI   idisp != Teuchos::null

  Teuchos::RCP<XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<XFluidFSI>(fluid_field(), true);

  // set idispnp in Xfluid
  if (idisp != Teuchos::null) xfluid->apply_struct_mesh_displacement(idisp);

  // set ivelnp in Xfluid
  if (ivel != Teuchos::null) xfluid->apply_struct_interface_velocities(ivel);

  fluid_field()->Solve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::RelaxationSolve(
    Teuchos::RCP<Epetra_Vector> idisp, double dt)
{
  std::cout << "WARNING: RelaxationSolve for XFEM useful?" << std::endl;

  // the displacement -> velocity conversion at the interface
  idisp->Scale(1. / dt);

  return fluid_field()->RelaxationSolve(idisp);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::extract_interface_forces()
{
  Teuchos::RCP<XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<XFluidFSI>(fluid_field(), true);
  return xfluid->extract_struct_interface_forces();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::extract_interface_velnp()
{
  FOUR_C_THROW("Robin stuff");
  Teuchos::RCP<XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<XFluidFSI>(fluid_field(), true);
  return xfluid->extract_struct_interface_velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::extract_interface_veln()
{
  Teuchos::RCP<XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<XFluidFSI>(fluid_field(), true);
  return xfluid->extract_struct_interface_veln();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::integrate_interface_shape()
{
  return fluid_field()->integrate_interface_shape();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<CORE::UTILS::ResultTest> ADAPTER::FluidXFEM::CreateFieldTest()
{
  return fluid_field()->CreateFieldTest();
}

FOUR_C_NAMESPACE_CLOSE
