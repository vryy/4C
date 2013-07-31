/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_xfem.cpp

\brief Fluid field adapter for xfem-fluids with moving boundaries

<pre>
Maintainer:  Benedikt Schott
             schott@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "adapter_coupling.H"
#include "ad_fld_fluid_xfem.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidXFEM::FluidXFEM(
        const Teuchos::ParameterList& prbdyn,
        std::string condname)
  : fluid_(prbdyn,false)
{
  icoupsf_ = Teuchos::rcp(new Coupling());
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidXFEM::Discretization()
{
  // returns the boundary discretization
  // REMARK:
  // the returned discretization has to match the structure discretization at the interface coupling (see FSI::Partitioned::Partitioned(const Epetra_Comm& comm) )
  // therefore return the boundary dis
  // this is similar to the matching of fluid dis and ale dis in case of ADAPTER::FluidALE
  return FluidField().Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const FLD::UTILS::MapExtractor& ADAPTER::FluidXFEM::Interface() const
{
  return *FluidField().Interface();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::PrepareTimeStep()
{
  FluidField().PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::Update()
{
  FluidField().Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::Output()
{
  FluidField().StatisticsAndOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidXFEM::ReadRestart(int step)
{
  FluidField().ReadRestart(step);
  return FluidField().Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::NonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                        Teuchos::RCP<Epetra_Vector> ivel)
{
  // if we have values at the interface we need to apply them

  // REMARK: for XFLUID idisp = Teuchos::null, ivel = Teuchos::null (called by fsi_fluid_xfem with default Teuchos::null)
  //         for XFSI   idisp != Teuchos::null

  // set idispnp in Xfluid
  if (idisp != Teuchos::null)
    FluidField().ApplyMeshDisplacement(idisp);

  // set ivelnp in Xfluid
  if (ivel != Teuchos::null)
    FluidField().ApplyInterfaceVelocities(ivel);


  FluidField().PrepareSolve();
  FluidField().Solve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::RelaxationSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                                                      double dt)
{

  std::cout << "WARNING: RelaxationSolve for XFEM useful?" << std::endl;

  // the displacement -> velocity conversion at the interface
  idisp->Scale(1./dt);

  return FluidField().RelaxationSolve(idisp);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::ExtractInterfaceForces()
{
  return FluidField().ExtractInterfaceForces();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::ExtractInterfaceVelnp()
{
  dserror("Robin stuff");
  return FluidField().ExtractInterfaceVelnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::ExtractInterfaceVeln()
{
  return FluidField().ExtractInterfaceVeln();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::IntegrateInterfaceShape()
{
  return FluidField().IntegrateInterfaceShape();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::FluidXFEM::CreateFieldTest()
{
  return FluidField().CreateFieldTest();
}
