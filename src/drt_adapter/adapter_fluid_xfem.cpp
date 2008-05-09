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
#include "adapter_utils.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidXFEM::FluidXFEM(
        const Teuchos::ParameterList& prbdyn,
        std::string condname,
        RCP<DRT::Discretization> soliddis)
  : fluid_(prbdyn,false)
{
  boundarydis_solidparalleldistrib_ = ADAPTER::UTILS::CreateDiscretizationFromCondition(soliddis, "FSICoupling", "Boundary", "BELE3");
  UTILS::SetupNDimExtractor(*boundarydis_solidparalleldistrib_,"FSICoupling",interface_);
  //UTILS::SetupNDimExtractor(*boundarydis_solidparalleldistrib_,"FREESURFCoupling",freesurface_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidXFEM::Discretization()
{
  return boundarydis_solidparalleldistrib_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const LINALG::MapExtractor& ADAPTER::FluidXFEM::Interface() const
{
  return interface_;
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
  cout << "ADAPTER::FluidXFEM::NonlinearSolve" << endl;
  // the solid-fluid coupling always matches geometrically
  const Epetra_Map* ifluidnodemap = FluidField().Discretization()->NodeRowMap();
  const Epetra_Map* isolidnodemap = boundarydis_solidparalleldistrib_->NodeRowMap();

  icoupsf_.SetupCouplingCheap(*FluidField().Discretization(),
                              *boundarydis_solidparalleldistrib_,
                              *ifluidnodemap,
                              *isolidnodemap);


  if (idisp!=Teuchos::null)
  {
    // if we have values at the interface we need to apply them
    FluidField().ApplyMeshDisplacement(SolidToFluid(idisp));
    FluidField().ApplyInterfaceVelocities(SolidToFluid(ivel));
  }
  cout << "applied interface displacement and velocity" << endl;

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

  return FluidField().RelaxationSolve(SolidToFluid(idisp));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::ExtractInterfaceForces()
{
  return FluidToSolid(FluidField().ExtractInterfaceForces());
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
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::IntegrateInterfaceShape()
{
  return FluidToSolid(FluidField().IntegrateInterfaceShape());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::FluidXFEM::CreateFieldTest()
{
  dserror("not implemented yet!");
  return FluidField().CreateFieldTest();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::SolidToFluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return icoupsf_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::SolidToFluid(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupsf_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::FluidToSolid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return icoupsf_.MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidXFEM::FluidToSolid(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupsf_.MasterToSlave(iv);
}


#endif
