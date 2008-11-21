/*----------------------------------------------------------------------*/
/*!
\file adapter_structure_constrained.cpp

\brief Adapter Layer for Structures with Algebraic Constraints

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15267
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
#include "../drt_inpar/drt_validparameters.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

/*======================================================================*/
/* constructor */
ADAPTER::StructureConstrained::StructureConstrained
(
  RCP<Structure> stru
)
: structure_(stru)
{    
  // make sure
  if (structure_ == null)
    dserror("Failed to create the underlying structural adapter");
  
  // build merged dof row map
  dofrowmap_ = LINALG::MergeMap(*(structure_->DofRowMap()),
                                *(structure_->GetConstraintManager()->GetConstraintMap()),
                                false);

  // set up interface between merged and single maps
  conmerger_.Setup(*dofrowmap_, 
                    structure_->DofRowMap(), 
                    structure_->GetConstraintManager()->GetConstraintMap());

  // initialise displacement increments to 0 (in words zero)
  // this variable in only used in monolithic FSI
  disinc_ = Teuchos::rcp(new Epetra_Vector(*dofrowmap_,true));
  
  //setup fsi-Interface
  DRT::UTILS::SetupNDimExtractor(*(structure_->Discretization()),"FSICoupling",dofrowmap_,interface_);
}


  
/*----------------------------------------------------------------------*/
/* */
RCP<const Epetra_Vector> ADAPTER::StructureConstrained::InitialGuess()
{
  //get initial guesses from structure and constraintmanager
  RCP<const Epetra_Vector> strucGuess = structure_->InitialGuess(); 
  RCP<const Epetra_Vector> lagrGuess = rcp(new Epetra_Vector(*(structure_->GetConstraintManager()->GetConstraintMap()),true));
  
  //merge stuff together
  RCP<Epetra_Vector> mergedGuess = rcp(new Epetra_Vector(*dofrowmap_,true));
  conmerger_.AddCondVector(strucGuess,mergedGuess);
  conmerger_.AddOtherVector(lagrGuess,mergedGuess);
  
  return mergedGuess;
}

/*----------------------------------------------------------------------*/
/* right-hand side alias the dynamic force residual */
RCP<const Epetra_Vector> ADAPTER::StructureConstrained::RHS()
{
  //get rhs-vectors from structure and constraintmanager
  RCP<const Epetra_Vector> struRHS = structure_->RHS(); 
  RCP<const Epetra_Vector> lagrRHS = structure_->GetConstraintManager()->GetError();
  
  //merge stuff together
  RCP<Epetra_Vector> mergedRHS = rcp(new Epetra_Vector(*dofrowmap_,true));
  conmerger_.AddCondVector(struRHS,mergedRHS);
  conmerger_.AddOtherVector(-1.0,lagrRHS,mergedRHS);
  
  return mergedRHS;
}


/*----------------------------------------------------------------------*/
/* get current displacements D_{n+1} */
RCP<const Epetra_Vector> ADAPTER::StructureConstrained::Dispnp()
{
  //get current state from structure and constraintmanager
  RCP<const Epetra_Vector> strudis = structure_->Dispnp(); 
  RCP<const Epetra_Vector> lagrmult = structure_->GetConstraintManager()->GetLagrMultVector();
  
  //merge stuff together
  RCP<Epetra_Vector> mergedstat = rcp(new Epetra_Vector(*dofrowmap_,true));
  conmerger_.AddCondVector(strudis,mergedstat);
  conmerger_.AddOtherVector(lagrmult,mergedstat);
  
  return mergedstat;
}


/*----------------------------------------------------------------------*/
/* get last converged displacements D_{n} */
RCP<const Epetra_Vector> ADAPTER::StructureConstrained::Dispn()
{
  //get last converged state from structure and constraintmanager
  RCP<const Epetra_Vector> strudis = structure_->Dispn(); 
  RCP<const Epetra_Vector> lagrmult = structure_->GetConstraintManager()->GetLagrMultVectorOld();
   
  //merge stuff together
   RCP<Epetra_Vector> mergedstat = rcp(new Epetra_Vector(*dofrowmap_,true));
   conmerger_.AddCondVector(strudis,mergedstat);
   conmerger_.AddOtherVector(lagrmult,mergedstat);
  
  return mergedstat;
}

/*----------------------------------------------------------------------*/
/* non-overlapping DOF map */
RCP<const Epetra_Map> ADAPTER::StructureConstrained::DofRowMap()
{
  return dofrowmap_;
}


/*----------------------------------------------------------------------*/
/* stiffness, i.e. force residual R_{n+1} differentiated
 * by displacements D_{n+1} */
RCP<LINALG::SparseMatrix> ADAPTER::StructureConstrained::SystemMatrix()
{
  //create empty large matrix and get small ones from structure and constraints
  RCP<LINALG::SparseMatrix> mergedmatrix = rcp(new LINALG::SparseMatrix(*dofrowmap_,dofrowmap_->NumMyElements()));
  RCP<LINALG::SparseMatrix> strustiff = structure_->SystemMatrix();
  RCP<LINALG::SparseMatrix> constiff = structure_->GetConstraintManager()->GetConstrMatrix();
  
  // Add matrices together
  mergedmatrix -> Add(*strustiff,false,1.0,0.0);
  mergedmatrix -> Add(*constiff,false,1.0,1.0);
  mergedmatrix -> Add(*constiff,true,1.0,1.0);
  mergedmatrix -> Complete(*dofrowmap_,*dofrowmap_);
  
  mergedmatrix->ApplyDirichlet( *(structure_->GetDBCMapExtractor()->CondMap()));
  
  return mergedmatrix;
}


/*----------------------------------------------------------------------*/
/* get discretisation */
RCP<DRT::Discretization> ADAPTER::StructureConstrained::Discretization()
{
  return structure_->Discretization();
}


/*----------------------------------------------------------------------*/
/* */
RCP<const Epetra_Vector> ADAPTER::StructureConstrained::FRobin()
{
  //return structure_->GetForceRobinFSI();
  return LINALG::CreateVector(*dofrowmap_, true);
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
  if (disinc_ != null)
    disinc_->PutScalar(0.0);

}


/*----------------------------------------------------------------------*/
/* build linear system stiffness matrix and rhs/force residual
 *
 * Monolithic FSI accesses the linearised structure problem. */
void ADAPTER::StructureConstrained::Evaluate(
  RCP<const Epetra_Vector> disp
)
{
  // 'initialize' structural displacement as null-pointer
  RCP<Epetra_Vector> dispstruct = Teuchos::null;
  
  // Compute residual increments, update total increments and update lagrange multipliers  
  if (disp != Teuchos::null)
  {
    // residual displacements (or iteration increments or iteratively incremental displacements)
    Teuchos::RCP<Epetra_Vector> disi = Teuchos::rcp(new Epetra_Vector(*disp));
    disi->Update(-1.0, *disinc_, 1.0);

    // update incremental displacement member to provided step increments
    // shortly: disinc_^<i> := disp^<i+1>
    disinc_->Update(1.0, *disp, 0.0);
    
    // Extract increments for lagr multipliers and do update
    RCP<Epetra_Vector> lagrincr = conmerger_.ExtractOtherVector(disi);
    structure_->UpdateIterIncrConstr(lagrincr);
    dispstruct = conmerger_.ExtractCondVector(disp);
  }
  // Hand down incremental displacements, 
  // structure_ will compute the residual increments on its own
  structure_->Evaluate(dispstruct);
  
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
  return *(LINALG::MergeMap(structure_->DomainMap(),
                            *(structure_->GetConstraintManager()->GetConstraintMap()),
                            false));
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
RCP<Epetra_Vector> ADAPTER::StructureConstrained::RelaxationSolve(
  RCP<Epetra_Vector> iforce
)
{
  return structure_->RelaxationSolve(iforce);
}

/*----------------------------------------------------------------------*/
/* extract interface displacements D_{n} */
RCP<Epetra_Vector> ADAPTER::StructureConstrained::ExtractInterfaceDispn()
{
  return structure_->ExtractInterfaceDispn();
}

/*----------------------------------------------------------------------*/
/* extract interface displacements D_{n+1} */
RCP<Epetra_Vector> ADAPTER::StructureConstrained::ExtractInterfaceDispnp()
{
  return structure_->ExtractInterfaceDispnp();
}

/*----------------------------------------------------------------------*/
/* extract external forces at interface F_{ext,n+1} */
RCP<Epetra_Vector> ADAPTER::StructureConstrained::ExtractInterfaceForces()
{
  return structure_->ExtractInterfaceForces();
}

/*----------------------------------------------------------------------*/
/* */
RCP<Epetra_Vector> ADAPTER::StructureConstrained::PredictInterfaceDispnp()
{
   return structure_->PredictInterfaceDispnp();
}


/*----------------------------------------------------------------------*/
/* */
void ADAPTER::StructureConstrained::ApplyInterfaceForces(
  RCP<Epetra_Vector> iforce
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
  RCP<Epetra_Vector> iforce,
  RCP<Epetra_Vector> ifluidvel
)
{
  dserror("Not implemented!");
}


/*----------------------------------------------------------------------*/
/* structural result test */
RCP<DRT::ResultTest> ADAPTER::StructureConstrained::CreateFieldTest()
{
  return structure_->CreateFieldTest();
}


/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
