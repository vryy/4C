/*----------------------------------------------------------------------*/
/*!
\file adapter_structure_constrained.cpp

\brief Structure field adapter

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "../drt_structure/strtimint_create.H"
#include "adapter_structure_constrained.H"
#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// further includes for StructureBaseAlgorithm:
#include "../drt_lib/drt_validparameters.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

/*======================================================================*/
/* constructor */
ADAPTER::StructureConstrained::StructureConstrained
(
  Teuchos::RCP<Structure> stru
)
: structure_(stru)
{
    
  // make sure
  if (structure_ == Teuchos::null)
    dserror("Failed to create structural integrator");

}


  
/*----------------------------------------------------------------------*/
/* */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureConstrained::InitialGuess()
{
  //TODO: make it big
  return structure_->InitialGuess();
  
}


/*----------------------------------------------------------------------*/
/* right-hand side alias the dynamic force residual */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureConstrained::RHS()
{
  //TODO: make it big
  return structure_->RHS();
}


/*----------------------------------------------------------------------*/
/* get current displacements D_{n+1} */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureConstrained::Dispnp()
{
  //TODO: make it big
  return structure_->Dispnp();
}


/*----------------------------------------------------------------------*/
/* get last converged displacements D_{n} */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureConstrained::Dispn()
{
  //TODO: make it big
  return structure_->Dispn();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// UNWANTED
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureConstrained::Dispnm()
{
  return structure_->Dispnm();
}


/*----------------------------------------------------------------------*/
/* non-overlapping DOF map */
Teuchos::RCP<const Epetra_Map> ADAPTER::StructureConstrained::DofRowMap()
{
  //TODO: := Mergedmap
  return structure_->DofRowMap();
}


/*----------------------------------------------------------------------*/
/* stiffness, i.e. force residual R_{n+1} differentiated
 * by displacements D_{n+1} */
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::StructureConstrained::SystemMatrix()
{
  //TODO: make it big
  return structure_->SystemMatrix();
}


/*----------------------------------------------------------------------*/
/* get discretisation */
Teuchos::RCP<DRT::Discretization> ADAPTER::StructureConstrained::Discretization()
{
  return structure_->Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// UNWANTED
double ADAPTER::StructureConstrained::DispIncrFactor()
{
  return -1; //structure_->DispIncrFactor();
}

/*----------------------------------------------------------------------*/
/* */
Teuchos::RCP<const Epetra_Vector> ADAPTER::StructureConstrained::FRobin()
{
  //return structure_->GetForceRobinFSI();
  return LINALG::CreateVector(*discret_->DofRowMap(), true);
}


/*----------------------------------------------------------------------*/
/* prepare time step */
void ADAPTER::StructureConstrained::PrepareTimeStep()
{
  // Note: MFSI requires a constant predictor. Otherwise the fields will get
  // out of sync.

  // predict
  structure_->PrepareTimeStep();

  // initialise incremental displacements
  if (disinc_ != Teuchos::null)
    disinc_->PutScalar(0.0);

}


/*----------------------------------------------------------------------*/
/* build linear system stiffness matrix and rhs/force residual
 *
 * Monolithic FSI accesses the linearised structure problem. */
void ADAPTER::StructureConstrained::Evaluate(
  Teuchos::RCP<const Epetra_Vector> disp
)
{
  //TODO: correct?
  structure_->Evaluate(disp);
}

/*----------------------------------------------------------------------*/
/* update time step */
void ADAPTER::StructureConstrained::Update()
{
  structure_->Update();
}


/*----------------------------------------------------------------------*/
/* output */
void ADAPTER::StructureConstrained::Output()
{
  structure_->Output();
}


/*----------------------------------------------------------------------*/
/* domain map */
const Epetra_Map& ADAPTER::StructureConstrained::DomainMap()
{
  return structure_->DomainMap();
}


/*----------------------------------------------------------------------*/
/* read restart */
void ADAPTER::StructureConstrained::ReadRestart(int step)
{
  structure_->ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/* find iteratively solution */
void ADAPTER::StructureConstrained::Solve()
{
  structure_->Solve();
}


/*----------------------------------------------------------------------*/
/* */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureConstrained::RelaxationSolve(
  Teuchos::RCP<Epetra_Vector> iforce
)
{
  return structure_->RelaxationSolve(iforce);
}

/*----------------------------------------------------------------------*/
/* extract interface displacements D_{n} */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureConstrained::ExtractInterfaceDispn()
{
  return structure_->ExtractInterfaceDispn();
}

/*----------------------------------------------------------------------*/
/* extract interface displacements D_{n+1} */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureConstrained::ExtractInterfaceDispnp()
{
  return structure_->ExtractInterfaceDispnp();
}

/*----------------------------------------------------------------------*/
/* extract external forces at interface F_{ext,n+1} */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureConstrained::ExtractInterfaceForces()
{
  return structure_->ExtractInterfaceForces();
}

/*----------------------------------------------------------------------*/
/* */
Teuchos::RCP<Epetra_Vector> ADAPTER::StructureConstrained::PredictInterfaceDispnp()
{
   return structure_->PredictInterfaceDispnp();
}


/*----------------------------------------------------------------------*/
/* */
void ADAPTER::StructureConstrained::ApplyInterfaceForces(
  Teuchos::RCP<Epetra_Vector> iforce
)
{
/*
  // Play it save. In the first iteration everything is already set up
  // properly. However, all following iterations need to calculate the
  // stiffness matrix here. Furthermore we are bound to reset fextm_
  // before we add our special contribution.
  // So we calculate the stiffness anyway (and waste the available
  // stiffness in the first iteration).
  structure_->ApplyExternalForce(interface_,iforce);
*/
  // This will add the provided interface force onto the residual forces
  // The sign convention of the interface force is external-force-like.
  structure_->ApplyInterfaceForces(iforce);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureConstrained::ApplyInterfaceRobinValue(
  Teuchos::RCP<Epetra_Vector> iforce,
  Teuchos::RCP<Epetra_Vector> ifluidvel
)
{
  dserror("Not impl.");
/*
  // get robin parameter and timestep
  double alphas  = params_->get<double>("alpha s",-1.);
  double dt      = params_->get<double>("delta time",-1.);
  double alphaf  = params_->get<double>("alpha f", 0.459);

  if (alphas<0. or dt<0.)
    dserror("couldn't get robin parameter alpha_s or time step size");

  // the RobinRHS is going to be:
  //
  // RobinRHS =
  //     - (alpha_s/dt)*(dis(n))
  //     - alpha_s*(1-alpha_f)*(fluidvel(n+1))
  //     + (1-alpha_f)*(iforce(n+1))

  // Attention: We must not change iforce here, because we would
  // implicitely change fextn_, too. fextn_ is needed to set fext_
  // after successfully reaching timestep end.
  // This is why an additional robin force vector is needed.

  Teuchos::RCP<Epetra_Vector> idisn  = interface_.ExtractCondVector(structure_->Disp());
  Teuchos::RCP<Epetra_Vector> frobin = interface_.ExtractCondVector(structure_->FRobin());

  // save robin coupling values in frobin vector (except iforce which
  // is passed separately)
  frobin->Update(alphas/dt,*idisn,alphas*(1-alphaf),*ifluidvel,0.0);

  interface_.InsertCondVector(frobin,structure_->FRobin());
  structure_->ApplyExternalForce(interface_,iforce);
*/
}


/*----------------------------------------------------------------------*/
/* structural result test */
Teuchos::RCP<DRT::ResultTest> ADAPTER::StructureConstrained::CreateFieldTest()
{
  return structure_->CreateFieldTest();
}


/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
