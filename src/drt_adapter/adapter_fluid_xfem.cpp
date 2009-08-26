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

#include "adapter_fluid_xfem.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidXFEM::FluidXFEM(
        const Teuchos::ParameterList& prbdyn,
        std::string condname)
  : fluid_(prbdyn,false)
{
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
  FluidField().PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::Evaluate(Teuchos::RCP<const Epetra_Vector> vel)
{
  std::cout << "ADAPTER::FluidXFEM::Evaluate()" << endl;
  FluidField().Evaluate(vel);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::Update()
{
  std::cout << "ADAPTER::FluidXFEM::Update()" << endl;
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
std::map<std::string,Teuchos::RCP<LINALG::SparseMatrix> > ADAPTER::FluidXFEM::CouplingMatrices()
{
  return FluidField().CouplingMatrices();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidXFEM::Output()
{
  FluidField().Output();

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
#ifdef DEBUG
  std::cout << "ADAPTER::FluidXFEM::NonlinearSolve" << endl;
#endif
  if (idisp!=Teuchos::null)
  {
    // if we have values at the interface we need to apply them
    FluidField().ApplyMeshDisplacement(idisp);
    FluidField().ApplyInterfaceVelocities(ivel);
  }
#ifdef DEBUG
  std::cout << "applied interface displacement and velocity" << endl;
#endif

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
void ADAPTER::FluidXFEM::RobinNonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                             Teuchos::RCP<Epetra_Vector> ivel,
                                             Teuchos::RCP<Epetra_Vector> iforce)
{
  dserror("you should not use robin-BC with XFEM, yet");
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
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::ExtractInterfaceForcesRobin()
{
  dserror("no Robin around here");
  return Teuchos::null;
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


#endif
