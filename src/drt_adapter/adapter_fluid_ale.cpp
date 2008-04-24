/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_ale.cpp

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "adapter_fluid_ale.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidAle::FluidAle(const Teuchos::ParameterList& prbdyn,
                                          std::string condname)
  : fluid_(prbdyn,true),
    ale_()
{
  icoupfa_.SetupConditionCoupling(*FluidField().Discretization(),
                                   FluidField().Interface(),
                                  *AleField().Discretization(),
                                   AleField().Interface(),
                                   condname);

  fscoupfa_.SetupConditionCoupling(*FluidField().Discretization(),
                                    FluidField().FreeSurface(),
                                   *AleField().Discretization(),
                                    AleField().FreeSurface(),
                                    "FREESURFCoupling");

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = FluidField().Discretization()->NodeRowMap();
  const Epetra_Map* alenodemap   = AleField().Discretization()->NodeRowMap();

  coupfa_.SetupCoupling(*FluidField().Discretization(),
                        *AleField().Discretization(),
                        *fluidnodemap,
                        *alenodemap);

  FluidField().SetMeshMap(coupfa_.MasterDofMap());

  // the ale matrix might be build just once
  AleField().BuildSystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidAle::Discretization()
{
  return FluidField().Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const LINALG::MapExtractor& ADAPTER::FluidAle::Interface() const
{
  return FluidField().Interface();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAle::PrepareTimeStep()
{
  FluidField().PrepareTimeStep();
  AleField().PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAle::Update()
{
  FluidField().Update();
  AleField().Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAle::Output()
{
  FluidField().Output();
  AleField().Output();

  FluidField().LiftDrag();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidAle::ReadRestart(int step)
{
  FluidField().ReadRestart(step);
  AleField().ReadRestart(step);
  return FluidField().Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAle::NonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                              Teuchos::RCP<Epetra_Vector> ivel)
{

  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  if (idisp!=Teuchos::null)
  {
    // if we have values at the interface we need to apply them
    AleField().ApplyInterfaceDisplacements(FluidToAle(idisp));
    if (Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO") != fsi_pseudo_structureale)
    {
      FluidField().ApplyInterfaceVelocities(ivel);
    }
  }

  if (FluidField().FreeSurface().Relevant())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = FluidField().Dispnp();
    Teuchos::RCP<Epetra_Vector> fsdispnp = FluidField().FreeSurface().ExtractCondVector(dispnp);
    AleField().ApplyFreeSurfaceDisplacements(fscoupfa_.MasterToSlave(fsdispnp));
  }

  // Note: We do not look for moving ale boundaries (outside the coupling
  // interface) on the fluid side. Thus if you prescribe time variable ale
  // Dirichlet conditions the according fluid Dirichlet conditions will not
  // notice.

  AleField().Solve();
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(AleField().ExtractDisplacement());
  FluidField().ApplyMeshDisplacement(fluiddisp);

  // no computation of fluid velocities in case only structure and ALE are to compute
  if (Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO") != fsi_pseudo_structureale)
  {
    FluidField().NonlinearSolve();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAle::RobinNonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                            Teuchos::RCP<Epetra_Vector> ivel,
                                            Teuchos::RCP<Epetra_Vector> iforce)
{
  // if we have values at the interface we need to apply them
  AleField().ApplyInterfaceDisplacements(FluidToAle(idisp));

  // pass coupling values for subsequent application at the interface of the
  // fluid domain
  FluidField().ApplyInterfaceRobinValue(ivel, iforce);

  // Note: We do not look for moving ale boundaries (outside the coupling
  // interface) on the fluid side. Thus if you prescribe time variable ale
  // Dirichlet conditions the according fluid Dirichlet conditions will not
  // notice.

  AleField().Solve();
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(AleField().ExtractDisplacement());
  FluidField().ApplyMeshDisplacement(fluiddisp);
  FluidField().NonlinearSolve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::RelaxationSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                                               double dt)
{
  // Here we have a mesh position independent of the
  // given trial vector, but still the grid velocity depends on the
  // trial vector only.

  // grid velocity
  AleField().ApplyInterfaceDisplacements(FluidToAle(idisp));

  AleField().Solve();
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(AleField().ExtractDisplacement());
  fluiddisp->Scale(1./dt);

  FluidField().ApplyMeshVelocity(fluiddisp);

  // grid position is done inside RelaxationSolve

  // the displacement -> velocity conversion at the interface
  idisp->Scale(1./dt);

  return FluidField().RelaxationSolve(idisp);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::ExtractInterfaceForces()
{
  return FluidField().ExtractInterfaceForces();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::ExtractInterfaceFluidVelocity()
{
  return FluidField().ExtractInterfaceFluidVelocity();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::IntegrateInterfaceShape()
{
  return FluidField().IntegrateInterfaceShape();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::FluidAle::CreateFieldTest()
{
  return FluidField().CreateFieldTest();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::AleToFluidField(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupfa_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::AleToFluidField(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const
{
  return icoupfa_.MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidAle::FluidToAle(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupfa_.MasterToSlave(iv);
}


#endif
