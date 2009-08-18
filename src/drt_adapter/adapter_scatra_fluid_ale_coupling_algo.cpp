/*----------------------------------------------------------------------*/
/*!
\file adapter_scatra_fluid_ale_coupling_algo.cpp

\brief Basis of all algorithms that perform a coupling between Navier-Stokes
       and scalar transport equations including deforming meshes

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "adapter_scatra_fluid_ale_coupling_algo.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraFluidAleCouplingAlgorithm::ScaTraFluidAleCouplingAlgorithm(
    Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn,
    const std::string condname
)
:  ScaTraFluidCouplingAlgorithm(comm, prbdyn, true), // yes, we need the ALE formulation
   AleBaseAlgorithm(prbdyn) // construct ale base algorithm as well
{
   // set up couplings
   icoupfa_.SetupConditionCoupling(*FluidField().Discretization(),
                                    FluidField().Interface().FSICondMap(),
                                   *AleField().Discretization(),
                                    AleField().Interface().FSICondMap(),
                                    condname);

   fscoupfa_.SetupConditionCoupling(*FluidField().Discretization(),
                                     FluidField().Interface().FSCondMap(),
                                    *AleField().Discretization(),
                                     AleField().Interface().FSCondMap(),
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

   return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraFluidAleCouplingAlgorithm::~ScaTraFluidAleCouplingAlgorithm()
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ScaTraFluidAleCouplingAlgorithm::ReadRestart(int step)
{
  ScaTraFluidCouplingAlgorithm::ReadRestart(step);
  AleField().ReadRestart(step); // add reading of ALE restart data

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ScaTraFluidAleCouplingAlgorithm::FluidAleNonlinearSolve(
    Teuchos::RCP<Epetra_Vector> idisp,
    Teuchos::RCP<Epetra_Vector> ivel)
{


  if (idisp!=Teuchos::null)
  {
    // if we have values at the interface we need to apply them
    AleField().ApplyInterfaceDisplacements(FluidToAle(idisp));
    FluidField().ApplyInterfaceVelocities(ivel);
  }

  if (FluidField().Interface().FSCondRelevant())
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = FluidField().Dispnp();
    Teuchos::RCP<Epetra_Vector> fsdispnp = FluidField().Interface().ExtractFSCondVector(dispnp);
    AleField().ApplyFreeSurfaceDisplacements(fscoupfa_.MasterToSlave(fsdispnp));
  }

  // Note: We do not look for moving ale boundaries (outside the coupling
  // interface) on the fluid side. Thus if you prescribe time variable ale
  // Dirichlet conditions the according fluid Dirichlet conditions will not
  // notice.

  AleField().Solve();
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(AleField().ExtractDisplacement());
  FluidField().ApplyMeshDisplacement(fluiddisp);

  // computation of fluid solution
  FluidField().NonlinearSolve();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::ScaTraFluidAleCouplingAlgorithm::AleToFluidField(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupfa_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::ScaTraFluidAleCouplingAlgorithm::AleToFluidField(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::ScaTraFluidAleCouplingAlgorithm::FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const
{
  return icoupfa_.MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::ScaTraFluidAleCouplingAlgorithm::FluidToAle(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupfa_.MasterToSlave(iv);
}

#endif // CCADISCRET
