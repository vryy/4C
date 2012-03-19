/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_xfem.cpp

\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_fluid/xfluid.H"

#include "adapter_fluid_xfem.H"
#include "adapter_coupling.H"


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
  return FluidField().Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const FLD::UTILS::MapExtractor& ADAPTER::FluidXFEM::Interface() const
{
  return FluidField().Interface();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::PrepareTimeStep()
{
  FLD::XFluid* ffield = dynamic_cast<FLD::XFluid*>(&(FluidField()));
  // update velocity n-1
  ffield->ivelnm_->Update(1.0,*ffield->iveln_,0.0);

  // update velocity n
  ffield->iveln_->Update(1.0,*ffield->ivelnp_,0.0);

  // update displacement n
  ffield->idispn_->Update(1.0,*ffield->idispnp_,0.0);

  FluidField().PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::Evaluate(
    Teuchos::RCP<Epetra_Vector> idispstepinc,
    Teuchos::RCP<const Epetra_Vector> fluidstepinc)
{
  if (idispstepinc == Teuchos::null)
    dserror("idispstepinc == Teuchos::null");
  if (fluidstepinc == Teuchos::null)
    dserror("fluidstepinc == Teuchos::null");

  FluidField().ApplyMeshDisplacementIncrement(idispstepinc);
  FluidField().Evaluate(fluidstepinc);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::Update()
{
  FluidField().Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidXFEM::DofRowMap()
{
  return FluidField().DofRowMap();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidXFEM::InitialGuess()
{
  return FluidField().InitialGuess();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidXFEM::RHS()
{
  return FluidField().RHS();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::FluidXFEM::SystemMatrix()
{
  return FluidField().SystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::MapExtractor> ADAPTER::FluidXFEM::GetDBCMapExtractor()
{
  return FluidField().GetDBCMapExtractor();
}


///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
//std::map<std::string,Teuchos::RCP<LINALG::SparseMatrix> > ADAPTER::FluidXFEM::CouplingMatrices()
//{
////  return XFluidField().CouplingMatrices();
//}
//
//
///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
//std::map<std::string,Teuchos::RCP<Epetra_Vector> > ADAPTER::FluidXFEM::CouplingVectors()
//{
////  return XFluidField().CouplingVectors();
//}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::Output()
{
  FluidField().StatisticsAndOutput();

  FluidField().LiftDrag();
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


  //if (FluidField().FreeSurface().Relevant())
  //{
  //  Teuchos::RCP<const Epetra_Vector> dispnp = FluidField().Dispnp();
  //  Teuchos::RCP<Epetra_Vector> fsdispnp = FluidField().FreeSurface().ExtractCondVector(dispnp);
  //  AleField().ApplyFreeSurfaceDisplacements(fscoupfa_.MasterToSlave(fsdispnp));
  //}

  // Note: We do not look for moving ale boundaries (outside the coupling
  // interface) on the fluid side. Thus if you prescribe time variable ale
  // Dirichlet conditions the according fluid Dirichlet conditions will not
  // notice.

  //AleField().Solve();
  //Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(AleField().ExtractDisplacement());
  //FluidField().ApplyMeshDisplacement(fluiddisp);


  FluidField().PrepareSolve();
  FluidField().NonlinearSolve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::RelaxationSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                                                      double dt)
{
  // Here we have a mesh position independent of the
  // given trial vector, but still the grid velocity depends on the
  // trial vector only.

  // grid velocity
  //AleField().ApplyInterfaceDisplacements(FluidToAle(idisp));

  //AleField().Solve();
  //Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(AleField().ExtractDisplacement());
  //fluiddisp->Scale(1./dt);

  //FluidField().ApplyMeshVelocity(fluiddisp);

  // grid position is done inside RelaxationSolve

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
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::ExtractInterfaceFluidVelocity()
{
  dserror("Robin stuff");
  return Teuchos::null;
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
  //dserror("not implemented yet!");
  return FluidField().CreateFieldTest();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
{
  FluidField().DisplacementToVelocity(fcx);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
{
  FluidField().VelocityToDisplacement(fcx);
}


#endif
