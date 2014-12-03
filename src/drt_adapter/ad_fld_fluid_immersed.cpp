/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_immersed.cpp

\brief Fluid field adapter for immersed fluids

<pre>
Maintainer:  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "adapter_coupling.H"
#include "ad_fld_fluid_immersed.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_adapter/ad_fld_fluid_fsi.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidImmersed::FluidImmersed(
        const Teuchos::ParameterList& prbdyn,
        std::string condname)
{
  fluid_=Teuchos::rcp(new FluidBaseAlgorithm(prbdyn,DRT::Problem::Instance()->FluidDynamicParams(),"fluid",false));
  fluidadapter_=fluid_->FluidField();
  icoupsf_ = Teuchos::rcp(new Coupling());
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidImmersed::Discretization()
{
  // returns the boundary discretization
  // REMARK:
  // the returned discretization has to match the structure discretization at the interface coupling (see FSI::Partitioned::Partitioned(const Epetra_Comm& comm) )
  // therefore return the boundary dis
  // this is similar to the matching of fluid dis and ale dis in case of ADAPTER::FluidALE
  return FluidField()->Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<FLD::UTILS::MapExtractor>const& ADAPTER::FluidImmersed::Interface() const
{
  return fluidadapter_->Interface();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImmersed::PrepareTimeStep()
{
  FluidField()->PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImmersed::Update()
{
  FluidField()->Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImmersed::Output()
{
  FluidField()->StatisticsAndOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidImmersed::ReadRestart(int step)
{
  FluidField()->ReadRestart(step);
  return FluidField()->Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImmersed::NonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                            Teuchos::RCP<Epetra_Vector> ivel)
{
  FluidField()->Solve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidImmersed::RelaxationSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                                                      double dt)
{

  std::cout << "WARNING: RelaxationSolve for XFEM useful?" << std::endl;

  // the displacement -> velocity conversion at the interface
  idisp->Scale(1./dt);

  return FluidField()->RelaxationSolve(idisp);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidImmersed::ExtractInterfaceForces()
{
  return FluidField()->ExtractInterfaceForces();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidImmersed::ExtractInterfaceVelnp()
{
  dserror("Robin stuff");
  return FluidField()->ExtractInterfaceVelnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidImmersed::ExtractInterfaceVeln()
{
  return FluidField()->ExtractInterfaceVeln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidImmersed::IntegrateInterfaceShape()
{
  return FluidField()->IntegrateInterfaceShape();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::FluidImmersed::CreateFieldTest()
{
  return FluidField()->CreateFieldTest();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImmersed::AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd)
{
  Teuchos::rcp_dynamic_cast<ADAPTER::FluidFSI>(FluidField())->FluidImplTimeInt()->AddDirichCond(maptoadd);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidImmersed::RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove)
{
  Teuchos::rcp_dynamic_cast<ADAPTER::FluidFSI>(FluidField())->FluidImplTimeInt()->RemoveDirichCond(maptoremove);
}
