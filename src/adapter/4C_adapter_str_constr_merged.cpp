/*----------------------------------------------------------------------*/
/*! \file

\brief Adapter Layer for Structures with Algebraic Constraints

\level 2


*/

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_adapter_str_constr_merged.hpp"

#include "4C_constraint_manager.hpp"
#include "4C_discretization_condition_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_structure_aux.hpp"
#include "4C_structure_timint_create.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*======================================================================*/
/* constructor */
Adapter::StructureConstrMerged::StructureConstrMerged(Teuchos::RCP<Structure> stru)
    : FSIStructureWrapper(stru), issetup_(false)
{
  // do nothing
}


/*----------------------------------------------------------------------*/
/* */
void Adapter::StructureConstrMerged::Setup()
{
  // call setup on time integrator
  StructureWrapper::Setup();

  // make sure
  if (structure_ == Teuchos::null)
    FOUR_C_THROW("Failed to create the underlying structural adapter");

  // build merged dof row map
  dofrowmap_ = Core::LinAlg::MergeMap(*(structure_->dof_row_map()),
      *(structure_->get_constraint_manager()->GetConstraintMap()), false);

  // set up interface between merged and single maps
  conmerger_ = Teuchos::rcp(new Core::LinAlg::MapExtractor);
  conmerger_->Setup(*dofrowmap_, structure_->dof_row_map(),
      structure_->get_constraint_manager()->GetConstraintMap());

  // setup fsi-Interface
  interface_ = Teuchos::rcp(new STR::MapExtractor);
  interface_->Setup(*discretization(), *dofrowmap_);

  issetup_ = true;
}


/*----------------------------------------------------------------------*/
/* */
Teuchos::RCP<const Epetra_Vector> Adapter::StructureConstrMerged::initial_guess()
{
  if (not issetup_) FOUR_C_THROW("Call Setup() first!");

  // get initial guesses from structure and constraintmanager
  Teuchos::RCP<const Epetra_Vector> strucGuess = structure_->initial_guess();
  Teuchos::RCP<const Epetra_Vector> lagrGuess = Teuchos::rcp(
      new Epetra_Vector(*(structure_->get_constraint_manager()->GetConstraintMap()), true));

  // merge stuff together
  Teuchos::RCP<Epetra_Vector> mergedGuess = Teuchos::rcp(new Epetra_Vector(*dofrowmap_, true));
  conmerger_->AddCondVector(strucGuess, mergedGuess);
  conmerger_->AddOtherVector(lagrGuess, mergedGuess);

  return mergedGuess;
}

/*----------------------------------------------------------------------*/
/* right-hand side alias the dynamic force residual */
Teuchos::RCP<const Epetra_Vector> Adapter::StructureConstrMerged::RHS()
{
  // get rhs-vectors from structure and constraintmanager
  Teuchos::RCP<const Epetra_Vector> struRHS = structure_->RHS();
  Teuchos::RCP<const Epetra_Vector> lagrRHS = structure_->get_constraint_manager()->GetError();

  // merge stuff together
  Teuchos::RCP<Epetra_Vector> mergedRHS = Teuchos::rcp(new Epetra_Vector(*dofrowmap_, true));
  conmerger_->AddCondVector(struRHS, mergedRHS);
  conmerger_->AddOtherVector(-1.0, lagrRHS, mergedRHS);

  return mergedRHS;
}


/*----------------------------------------------------------------------*/
/* get current displacements D_{n+1} */
Teuchos::RCP<const Epetra_Vector> Adapter::StructureConstrMerged::Dispnp() const
{
  // get current state from structure and constraintmanager
  Teuchos::RCP<const Epetra_Vector> strudis = structure_->Dispnp();
  Teuchos::RCP<const Epetra_Vector> lagrmult =
      structure_->get_constraint_manager()->GetLagrMultVector();

  // merge stuff together
  Teuchos::RCP<Epetra_Vector> mergedstat = Teuchos::rcp(new Epetra_Vector(*dofrowmap_, true));
  conmerger_->AddCondVector(strudis, mergedstat);
  conmerger_->AddOtherVector(lagrmult, mergedstat);

  return mergedstat;
}


/*----------------------------------------------------------------------*/
/* get last converged displacements D_{n} */
Teuchos::RCP<const Epetra_Vector> Adapter::StructureConstrMerged::Dispn() const
{
  // get last converged state from structure and constraintmanager
  Teuchos::RCP<const Epetra_Vector> strudis = structure_->Dispn();
  Teuchos::RCP<const Epetra_Vector> lagrmult =
      structure_->get_constraint_manager()->get_lagr_mult_vector_old();

  // merge stuff together
  Teuchos::RCP<Epetra_Vector> mergedstat = Teuchos::rcp(new Epetra_Vector(*dofrowmap_, true));
  conmerger_->AddCondVector(strudis, mergedstat);
  conmerger_->AddOtherVector(lagrmult, mergedstat);

  return mergedstat;
}

/*----------------------------------------------------------------------*/
/* get last converged velocities V_{n} with zeroed Lagrange multiplier */
Teuchos::RCP<const Epetra_Vector> Adapter::StructureConstrMerged::Veln() const
{
  // get last converged state from structure and constraintmanager
  Teuchos::RCP<const Epetra_Vector> strudis = structure_->Veln();
  Teuchos::RCP<const Epetra_Vector> lagrmult = Teuchos::rcp(new Epetra_Vector(
      structure_->get_constraint_manager()->get_lagr_mult_vector_old()->Map(), true));

  // merge stuff together
  Teuchos::RCP<Epetra_Vector> mergedstat = Teuchos::rcp(new Epetra_Vector(*dofrowmap_, true));
  conmerger_->AddCondVector(strudis, mergedstat);
  conmerger_->AddOtherVector(lagrmult, mergedstat);

  return mergedstat;
}

/*----------------------------------------------------------------------*/
/* get last converged accelerations A_{n} with zeroed Lagrange multiplier */
Teuchos::RCP<const Epetra_Vector> Adapter::StructureConstrMerged::Accn() const
{
  // get last converged state from structure and constraintmanager
  Teuchos::RCP<const Epetra_Vector> strudis = structure_->Accn();
  Teuchos::RCP<const Epetra_Vector> lagrmult = Teuchos::rcp(new Epetra_Vector(
      structure_->get_constraint_manager()->get_lagr_mult_vector_old()->Map(), true));

  // merge stuff together
  Teuchos::RCP<Epetra_Vector> mergedstat = Teuchos::rcp(new Epetra_Vector(*dofrowmap_, true));
  conmerger_->AddCondVector(strudis, mergedstat);
  conmerger_->AddOtherVector(lagrmult, mergedstat);

  return mergedstat;
}

/*----------------------------------------------------------------------*/
/* non-overlapping DOF map */
Teuchos::RCP<const Epetra_Map> Adapter::StructureConstrMerged::dof_row_map() { return dofrowmap_; }


/*----------------------------------------------------------------------*/
/* stiffness, i.e. force residual R_{n+1} differentiated
 * by displacements D_{n+1} */
Teuchos::RCP<Core::LinAlg::SparseMatrix> Adapter::StructureConstrMerged::SystemMatrix()
{
  // create empty large matrix and get small ones from structure and constraints
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mergedmatrix =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dofrowmap_, 81));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> strustiff = structure_->SystemMatrix();
  strustiff->Complete();

  Teuchos::RCP<Core::LinAlg::SparseOperator> constiff =
      structure_->get_constraint_manager()->GetConstrMatrix();
  constiff->Complete();

  // Add matrices together
  mergedmatrix->Add(*strustiff, false, 1.0, 0.0);
  mergedmatrix->Add(*constiff, false, 1.0, 1.0);
  mergedmatrix->Add(*constiff, true, 1.0, 1.0);
  mergedmatrix->Complete(*dofrowmap_, *dofrowmap_);

  mergedmatrix->ApplyDirichlet(*(structure_->GetDBCMapExtractor()->CondMap()));

  return mergedmatrix;
}


/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>
Adapter::StructureConstrMerged::BlockSystemMatrix()
{
  FOUR_C_THROW("constrained BlockSparseMatrix never to be implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/* build linear system stiffness matrix and rhs/force residual
 *
 * Monolithic FSI accesses the linearised structure problem. */
void Adapter::StructureConstrMerged::Evaluate(Teuchos::RCP<const Epetra_Vector> dispstepinc)
{
  // 'initialize' structural displacement as null-pointer
  Teuchos::RCP<Epetra_Vector> dispstructstepinc = Teuchos::null;

  // Compute residual increments, update total increments and update lagrange multipliers
  if (dispstepinc != Teuchos::null)
  {
    // Extract increments for lagr multipliers and do update
    Teuchos::RCP<Epetra_Vector> lagrincr = conmerger_->ExtractOtherVector(dispstepinc);
    structure_->update_iter_incr_constr(lagrincr);
    dispstructstepinc = conmerger_->ExtractCondVector(dispstepinc);
  }
  // Hand down incremental displacements,
  // structure_ will compute the residual increments on its own
  structure_->Evaluate(dispstructstepinc);
}


/*----------------------------------------------------------------------*/
/* domain map */
const Epetra_Map& Adapter::StructureConstrMerged::DomainMap() const
{
  return *(Core::LinAlg::MergeMap(
      structure_->DomainMap(), *(structure_->get_constraint_manager()->GetConstraintMap()), false));
}

/*----------------------------------------------------------------------*/
// Apply interface forces
void Adapter::StructureConstrMerged::apply_interface_forces_temporary_deprecated(
    Teuchos::RCP<Epetra_Vector> iforce)
{
  // create vector with displacement and constraint DOFs
  Teuchos::RCP<Epetra_Vector> fifc = Core::LinAlg::CreateVector(*dof_row_map(), true);

  // insert interface forces
  interface_->AddFSICondVector(iforce, fifc);

  // extract the force values from the displacement DOFs only
  Teuchos::RCP<Epetra_Vector> fifcdisp = Core::LinAlg::CreateVector(*conmerger_->CondMap(), true);
  conmerger_->ExtractCondVector(fifc, fifcdisp);

  // set interface forces within the structural time integrator
  SetForceInterface(fifcdisp);

  prepare_partition_step();

  return;
}

FOUR_C_NAMESPACE_CLOSE
