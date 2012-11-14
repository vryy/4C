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
/* headers */
#include "../drt_structure/strtimint_create.H"
#include "ad_str_constr_merged.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_constraint/constraint_manager.H"
#include "../drt_structure/stru_aux.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// further includes for StructureBaseAlgorithm:
#include "../drt_inpar/drt_validparameters.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

/*======================================================================*/
/* constructor */
ADAPTER::StructureConstrMerged::StructureConstrMerged
(
  RCP<Structure> stru
)
: FSIStructureWrapper(stru)
{
  // make sure
  if (structure_ == null)
    dserror("Failed to create the underlying structural adapter");

  // build merged dof row map
  dofrowmap_ = LINALG::MergeMap(*(structure_->DofRowMap()),
                                *(structure_->GetConstraintManager()->GetConstraintMap()),
                                false);

  // set up interface between merged and single maps
  conmerger_ = Teuchos::rcp(new LINALG::MapExtractor);
  conmerger_->Setup(*dofrowmap_,
                    structure_->DofRowMap(),
                    structure_->GetConstraintManager()->GetConstraintMap());

  //setup fsi-Interface
  interface_ = Teuchos::rcp(new STR::AUX::MapExtractor);
  interface_->Setup(*Discretization(), *dofrowmap_);
}



/*----------------------------------------------------------------------*/
/* */
RCP<const Epetra_Vector> ADAPTER::StructureConstrMerged::InitialGuess()
{
  //get initial guesses from structure and constraintmanager
  RCP<const Epetra_Vector> strucGuess = structure_->InitialGuess();
  RCP<const Epetra_Vector> lagrGuess = Teuchos::rcp(new Epetra_Vector(*(structure_->GetConstraintManager()->GetConstraintMap()),true));

  //merge stuff together
  RCP<Epetra_Vector> mergedGuess = Teuchos::rcp(new Epetra_Vector(*dofrowmap_,true));
  conmerger_->AddCondVector(strucGuess,mergedGuess);
  conmerger_->AddOtherVector(lagrGuess,mergedGuess);

  return mergedGuess;
}

/*----------------------------------------------------------------------*/
/* right-hand side alias the dynamic force residual */
RCP<const Epetra_Vector> ADAPTER::StructureConstrMerged::RHS()
{
  //get rhs-vectors from structure and constraintmanager
  RCP<const Epetra_Vector> struRHS = structure_->RHS();
  RCP<const Epetra_Vector> lagrRHS = structure_->GetConstraintManager()->GetError();

  //merge stuff together
  RCP<Epetra_Vector> mergedRHS = Teuchos::rcp(new Epetra_Vector(*dofrowmap_,true));
  conmerger_->AddCondVector(struRHS,mergedRHS);
  conmerger_->AddOtherVector(-1.0,lagrRHS,mergedRHS);

  return mergedRHS;
}


/*----------------------------------------------------------------------*/
/* get current displacements D_{n+1} */
RCP<const Epetra_Vector> ADAPTER::StructureConstrMerged::Dispnp()
{
  //get current state from structure and constraintmanager
  RCP<const Epetra_Vector> strudis = structure_->Dispnp();
  RCP<const Epetra_Vector> lagrmult = structure_->GetConstraintManager()->GetLagrMultVector();

  //merge stuff together
  RCP<Epetra_Vector> mergedstat = Teuchos::rcp(new Epetra_Vector(*dofrowmap_,true));
  conmerger_->AddCondVector(strudis,mergedstat);
  conmerger_->AddOtherVector(lagrmult,mergedstat);

  return mergedstat;
}


/*----------------------------------------------------------------------*/
/* get last converged displacements D_{n} */
RCP<const Epetra_Vector> ADAPTER::StructureConstrMerged::Dispn()
{
  //get last converged state from structure and constraintmanager
  RCP<const Epetra_Vector> strudis = structure_->Dispn();
  RCP<const Epetra_Vector> lagrmult = structure_->GetConstraintManager()->GetLagrMultVectorOld();

  //merge stuff together
   RCP<Epetra_Vector> mergedstat = Teuchos::rcp(new Epetra_Vector(*dofrowmap_,true));
   conmerger_->AddCondVector(strudis,mergedstat);
   conmerger_->AddOtherVector(lagrmult,mergedstat);

  return mergedstat;
}

/*----------------------------------------------------------------------*/
/* non-overlapping DOF map */
RCP<const Epetra_Map> ADAPTER::StructureConstrMerged::DofRowMap()
{
  return dofrowmap_;
}


/*----------------------------------------------------------------------*/
/* stiffness, i.e. force residual R_{n+1} differentiated
 * by displacements D_{n+1} */
RCP<LINALG::SparseMatrix> ADAPTER::StructureConstrMerged::SystemMatrix()
{
  //create empty large matrix and get small ones from structure and constraints
  RCP<LINALG::SparseMatrix> mergedmatrix = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_, 81));
  RCP<LINALG::SparseMatrix> strustiff = structure_->SystemMatrix();
  strustiff->Complete();
  
  RCP<LINALG::SparseOperator> constiff = structure_->GetConstraintManager()->GetConstrMatrix();
  constiff->Complete();

  // Add matrices together
  mergedmatrix -> Add(*strustiff,false,1.0,0.0);
  mergedmatrix -> Add(*constiff,false,1.0,1.0);
  mergedmatrix -> Add(*constiff,true,1.0,1.0);
  mergedmatrix -> Complete(*dofrowmap_,*dofrowmap_);

  mergedmatrix -> ApplyDirichlet( *(structure_->GetDBCMapExtractor()->CondMap()));

  return mergedmatrix;
}


/*----------------------------------------------------------------------*/
RCP<LINALG::BlockSparseMatrixBase> ADAPTER::StructureConstrMerged::BlockSystemMatrix()
{
  dserror("constrained BlockSparseMatrix never to be implemented");
  return null;
}


/*----------------------------------------------------------------------*/
/* build linear system stiffness matrix and rhs/force residual
 *
 * Monolithic FSI accesses the linearised structure problem. */
void ADAPTER::StructureConstrMerged::Evaluate(
  RCP<const Epetra_Vector> dispstepinc
)
{
  // 'initialize' structural displacement as null-pointer
  RCP<Epetra_Vector> dispstructstepinc = Teuchos::null;

  // Compute residual increments, update total increments and update lagrange multipliers
  if (dispstepinc != Teuchos::null)
  {
    // Extract increments for lagr multipliers and do update
    RCP<Epetra_Vector> lagrincr = conmerger_->ExtractOtherVector(dispstepinc);
    structure_->UpdateIterIncrConstr(lagrincr);
    dispstructstepinc = conmerger_->ExtractCondVector(dispstepinc);
  }
  // Hand down incremental displacements,
  // structure_ will compute the residual increments on its own
  structure_->Evaluate(dispstructstepinc);

}


/*----------------------------------------------------------------------*/
/* domain map */
const Epetra_Map& ADAPTER::StructureConstrMerged::DomainMap()
{
  return *(LINALG::MergeMap(structure_->DomainMap(),
                            *(structure_->GetConstraintManager()->GetConstraintMap()),
                            false));
}
