/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for xfem-fluids with moving boundaries

\level 1


*/
/*----------------------------------------------------------------------*/
#include "baci_adapter_fld_fluid_xfem.H"

#include "baci_adapter_fld_fluid_xfsi.H"
#include "baci_coupling_adapter.H"
#include "baci_fluid_utils_mapextractor.H"
#include "baci_global_data.H"
#include "baci_inpar_validparameters.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidXFEM::FluidXFEM(const Teuchos::ParameterList& prbdyn, std::string condname)
    : fluid_(Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(prbdyn,
                              DRT::Problem::Instance()->FluidDynamicParams(), "fluid", false))
                 ->FluidField())
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidXFEM::Discretization()
{
  // returns the boundary discretization
  // REMARK:
  // the returned discretization has to match the structure discretization at the interface coupling
  // (see FSI::Partitioned::Partitioned(const Epetra_Comm& comm) ) therefore return the boundary dis
  // this is similar to the matching of fluid dis and ale dis in case of ADAPTER::FluidALE
  return FluidField()->Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidXFEM::BoundaryDiscretization()
{
  // returns the boundary discretization
  // REMARK:
  // the returned discretization has to match the structure discretization at the interface coupling
  // (see FSI::Partitioned::Partitioned(const Epetra_Comm& comm) ) therefore return the boundary dis
  // this is similar to the matching of fluid dis and ale dis in case of ADAPTER::FluidALE

  Teuchos::RCP<XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<XFluidFSI>(FluidField(), true);

  return xfluid->BoundaryDiscretization();
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

  Teuchos::RCP<XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<XFluidFSI>(FluidField(), true);

  return xfluid->StructInterface();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::PrepareTimeStep() { FluidField()->PrepareTimeStep(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::Update() { FluidField()->Update(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::Output() { FluidField()->StatisticsAndOutput(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidXFEM::ReadRestart(int step)
{
  FluidField()->ReadRestart(step);
  return FluidField()->Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::NonlinearSolve(
    Teuchos::RCP<Epetra_Vector> idisp, Teuchos::RCP<Epetra_Vector> ivel)
{
  // if we have values at the interface we need to apply them

  // REMARK: for XFLUID idisp = Teuchos::null, ivel = Teuchos::null (called by fsi_fluid_xfem with
  // default Teuchos::null)
  //         for XFSI   idisp != Teuchos::null

  Teuchos::RCP<XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<XFluidFSI>(FluidField(), true);

  // set idispnp in Xfluid
  if (idisp != Teuchos::null) xfluid->ApplyStructMeshDisplacement(idisp);

  // set ivelnp in Xfluid
  if (ivel != Teuchos::null) xfluid->ApplyStructInterfaceVelocities(ivel);

  FluidField()->Solve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::RelaxationSolve(
    Teuchos::RCP<Epetra_Vector> idisp, double dt)
{
  std::cout << "WARNING: RelaxationSolve for XFEM useful?" << std::endl;

  // the displacement -> velocity conversion at the interface
  idisp->Scale(1. / dt);

  return FluidField()->RelaxationSolve(idisp);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::ExtractInterfaceForces()
{
  Teuchos::RCP<XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<XFluidFSI>(FluidField(), true);
  return xfluid->ExtractStructInterfaceForces();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::ExtractInterfaceVelnp()
{
  dserror("Robin stuff");
  Teuchos::RCP<XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<XFluidFSI>(FluidField(), true);
  return xfluid->ExtractStructInterfaceVelnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::ExtractInterfaceVeln()
{
  Teuchos::RCP<XFluidFSI> xfluid = Teuchos::rcp_dynamic_cast<XFluidFSI>(FluidField(), true);
  return xfluid->ExtractStructInterfaceVeln();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::IntegrateInterfaceShape()
{
  return FluidField()->IntegrateInterfaceShape();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::FluidXFEM::CreateFieldTest()
{
  return FluidField()->CreateFieldTest();
}

BACI_NAMESPACE_CLOSE
