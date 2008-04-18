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
#include "../drt_lib/drt_validparameters.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "adapter_fluid_xfem.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidXFEM::FluidXFEM(const Teuchos::ParameterList& prbdyn,
                                          std::string condname)
  : fluid_(prbdyn,false),
    boundary_()
{
//  icoupfa_.SetupConditionCoupling(*FluidField().Discretization(),
//                                   FluidField().Interface(),
//                                  *AleField().Discretization(),
//                                   AleField().Interface(),
//                                   condname);

  //FSI::Coupling& coupfa = FluidAleFieldCoupling();

  // the fluid-ale coupling always matches
  //const Epetra_Map* fluidnodemap = FluidField().Discretization()->NodeRowMap();
  //const Epetra_Map* alenodemap   = AleField().Discretization()->NodeRowMap();

//  coupfa_.SetupCoupling(*FluidField().Discretization(),
//                        *AleField().Discretization(),
//                        *fluidnodemap,
//                        *alenodemap);

  //FluidField().SetMeshMap(coupfa_.MasterDofMap());

  // the ale matrix is build just once
  //AleField().BuildSystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidXFEM::Discretization()
{
  return FluidField().Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const LINALG::MapExtractor& ADAPTER::FluidXFEM::Interface() const
{
  return FluidField().Interface();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::PrepareTimeStep()
{
  FluidField().PrepareTimeStep();
  //AleField().PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::Update()
{
  FluidField().Update();
  //AleField().Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::Output()
{
  FluidField().Output();
  //AleField().Output();

  FluidField().LiftDrag();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidXFEM::ReadRestart(int step)
{
  FluidField().ReadRestart(step);
  //AleField().ReadRestart(step);
  return FluidField().Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::NonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                          Teuchos::RCP<Epetra_Vector> ivel)
{
  if (idisp!=Teuchos::null)
  {
    // if we have values at the interface we need to apply them
    //AleField().ApplyInterfaceDisplacements(FluidToAle(idisp));
    FluidField().ApplyInterfaceVelocities(ivel);
  }

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


///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
//Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::AleToFluidField(Teuchos::RCP<Epetra_Vector> iv) const
//{
//  return coupfa_.SlaveToMaster(iv);
//}
//
//
///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
//Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::AleToFluidField(Teuchos::RCP<const Epetra_Vector> iv) const
//{
//  return coupfa_.SlaveToMaster(iv);
//}
//
//
///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
//Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const
//{
//  return icoupfa_.MasterToSlave(iv);
//}
//
//
///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
//Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::FluidToAle(Teuchos::RCP<const Epetra_Vector> iv) const
//{
//  return icoupfa_.MasterToSlave(iv);
//}


#endif
