/*----------------------------------------------------------------------*/
/*!
\file fpsi_robinneumann.cpp

\brief Solve FSI problems using a Robin-Neumann partitioning approach

\maintainer Andreas Rauch
            rauch@lnm.mw.tum.de

\level 3
*/
/*----------------------------------------------------------------------*/


#include "fpsi_robinneumann.H"
#include "../drt_adapter/ad_str_fpsiwrapper.H"
#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FPSI::RobinNeumann::RobinNeumann(const Epetra_Comm& comm,
                                 const Teuchos::ParameterList& fpsidynparams)
    : Partitioned(comm,fpsidynparams)
{

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::RobinNeumann::FSIOp(const Epetra_Vector &x, Epetra_Vector &F)
{
  const Teuchos::RCP<Epetra_Vector> idispn = Teuchos::rcp(new Epetra_Vector(x));

  const Teuchos::RCP<Epetra_Vector> iforce  = FluidOp(idispn);
  //const Teuchos::RCP<Epetra_Vector> idispnp = PoroOp(iforce);

  //F.Update(1.0, *idispnp, -1.0, *idispn, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
FPSI::RobinNeumann::FluidOp(Teuchos::RCP<Epetra_Vector> idisp)
{
  // Calculate the interface velocity to be applied as dirichlet bc to
  // the fluid field
  /*const Teuchos::RCP<Epetra_Vector> ivel = InterfaceVelocity(idisp);

  FluidField()->NonlinearSolve(StructToFluid(idisp),StructToFluid(ivel));*/

  return FluidToStruct(FluidField()->ExtractInterfaceForces());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void
FPSI::RobinNeumann::PoroOp(Teuchos::RCP<Epetra_Vector> iforce)
{
  //PoroField()->StructureField()->ApplyInterfaceForces(iforce);
  PoroField()->Solve();

  /* Return the interface velocity of the fluid in the porous medium.
     The normal component of the dirichlet boundary condition for the
     fluid at the interface is determined by both, fluid velocity and
     structural velocity and furthermore by the porosity.
     The tangential component is determined by the Beavers-Joseph (BJ)
     interface condition. To evaluate the BJ condition we need to ask
     the porous structure for its interface displacements (to be able
     to calculate the interface velocity of the structure). We also
     need to know some Parameters like the Beavers-Joseph-coefficient,
     the Permeability and the dynamic viscosity.
     In addition the strain rate tensor and thus the spatial derivatives
     of the fluid velocity at the interface are necessary.
     All this is assembled and evaluated in the method "InterfaceVelocity"
   */
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FPSI::RobinNeumann::InitialGuess()
{
  if (displacementcoupling_)
  {
    // predict displacement
    return PoroField()->StructureField()->PredictInterfaceDispnp();
  }
  return Teuchos::null;
}


