/*----------------------------------------------------------------------*/
/*! \file

\brief Adapter Layer for Structures with Algebraic Constraints

\level 2


*/

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_adapter_str_constr_merged.hpp"

#include "4C_constraint_manager.hpp"
#include "4C_fem_condition_utils.hpp"
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
void Adapter::StructureConstrMerged::setup()
{
  // call setup on time integrator
  StructureWrapper::setup();

  // make sure
  if (structure_ == Teuchos::null)
    FOUR_C_THROW("Failed to create the underlying structural adapter");

  // build merged dof row map
  dofrowmap_ = Core::LinAlg::merge_map(*(structure_->dof_row_map()),
      *(structure_->get_constraint_manager()->get_constraint_map()), false);

  // set up interface between merged and single maps
  conmerger_ = Teuchos::RCP(new Core::LinAlg::MapExtractor);
  conmerger_->setup(*dofrowmap_, structure_->dof_row_map(),
      structure_->get_constraint_manager()->get_constraint_map());

  // setup fsi-Interface
  interface_ = Teuchos::RCP(new Solid::MapExtractor);
  interface_->setup(*discretization(), *dofrowmap_);

  issetup_ = true;
}


/*----------------------------------------------------------------------*/
/* */
Teuchos::RCP<const Core::LinAlg::Vector<double>> Adapter::StructureConstrMerged::initial_guess()
{
  if (not issetup_) FOUR_C_THROW("Call setup() first!");

  // get initial guesses from structure and constraintmanager
  Teuchos::RCP<const Core::LinAlg::Vector<double>> strucGuess = structure_->initial_guess();
  Teuchos::RCP<const Core::LinAlg::Vector<double>> lagrGuess =
      Teuchos::RCP(new Core::LinAlg::Vector<double>(
          *(structure_->get_constraint_manager()->get_constraint_map()), true));

  // merge stuff together
  Teuchos::RCP<Core::LinAlg::Vector<double>> mergedGuess =
      Teuchos::RCP(new Core::LinAlg::Vector<double>(*dofrowmap_, true));
  conmerger_->add_cond_vector(*strucGuess, *mergedGuess);
  conmerger_->add_other_vector(*lagrGuess, *mergedGuess);

  return mergedGuess;
}

/*----------------------------------------------------------------------*/
/* right-hand side alias the dynamic force residual */
Teuchos::RCP<const Core::LinAlg::Vector<double>> Adapter::StructureConstrMerged::rhs()
{
  // get rhs-vectors from structure and constraintmanager
  Teuchos::RCP<const Core::LinAlg::Vector<double>> struRHS = structure_->rhs();
  Teuchos::RCP<const Core::LinAlg::Vector<double>> lagrRHS =
      structure_->get_constraint_manager()->get_error();

  // merge stuff together
  Teuchos::RCP<Core::LinAlg::Vector<double>> mergedRHS =
      Teuchos::RCP(new Core::LinAlg::Vector<double>(*dofrowmap_, true));
  conmerger_->add_cond_vector(*struRHS, *mergedRHS);
  conmerger_->add_other_vector(-1.0, *lagrRHS, *mergedRHS);

  return mergedRHS;
}


/*----------------------------------------------------------------------*/
/* get current displacements D_{n+1} */
Teuchos::RCP<const Core::LinAlg::Vector<double>> Adapter::StructureConstrMerged::dispnp() const
{
  // get current state from structure and constraintmanager
  Teuchos::RCP<const Core::LinAlg::Vector<double>> strudis = structure_->dispnp();
  Teuchos::RCP<const Core::LinAlg::Vector<double>> lagrmult =
      structure_->get_constraint_manager()->get_lagr_mult_vector();

  // merge stuff together
  Teuchos::RCP<Core::LinAlg::Vector<double>> mergedstat =
      Teuchos::RCP(new Core::LinAlg::Vector<double>(*dofrowmap_, true));
  conmerger_->add_cond_vector(*strudis, *mergedstat);
  conmerger_->add_other_vector(*lagrmult, *mergedstat);

  return mergedstat;
}


/*----------------------------------------------------------------------*/
/* get last converged displacements D_{n} */
Teuchos::RCP<const Core::LinAlg::Vector<double>> Adapter::StructureConstrMerged::dispn() const
{
  // get last converged state from structure and constraintmanager
  Teuchos::RCP<const Core::LinAlg::Vector<double>> strudis = structure_->dispn();
  Teuchos::RCP<const Core::LinAlg::Vector<double>> lagrmult =
      structure_->get_constraint_manager()->get_lagr_mult_vector_old();

  // merge stuff together
  Teuchos::RCP<Core::LinAlg::Vector<double>> mergedstat =
      Teuchos::RCP(new Core::LinAlg::Vector<double>(*dofrowmap_, true));
  conmerger_->add_cond_vector(*strudis, *mergedstat);
  conmerger_->add_other_vector(*lagrmult, *mergedstat);

  return mergedstat;
}

/*----------------------------------------------------------------------*/
/* get last converged velocities V_{n} with zeroed Lagrange multiplier */
Teuchos::RCP<const Core::LinAlg::Vector<double>> Adapter::StructureConstrMerged::veln() const
{
  // get last converged state from structure and constraintmanager
  Teuchos::RCP<const Core::LinAlg::Vector<double>> strudis = structure_->veln();
  Teuchos::RCP<const Core::LinAlg::Vector<double>> lagrmult =
      Teuchos::RCP(new Core::LinAlg::Vector<double>(
          structure_->get_constraint_manager()->get_lagr_mult_vector_old()->Map(), true));

  // merge stuff together
  Teuchos::RCP<Core::LinAlg::Vector<double>> mergedstat =
      Teuchos::RCP(new Core::LinAlg::Vector<double>(*dofrowmap_, true));
  conmerger_->add_cond_vector(*strudis, *mergedstat);
  conmerger_->add_other_vector(*lagrmult, *mergedstat);

  return mergedstat;
}

/*----------------------------------------------------------------------*/
/* get last converged accelerations A_{n} with zeroed Lagrange multiplier */
Teuchos::RCP<const Core::LinAlg::Vector<double>> Adapter::StructureConstrMerged::accn() const
{
  // get last converged state from structure and constraintmanager
  Teuchos::RCP<const Core::LinAlg::Vector<double>> strudis = structure_->accn();
  Teuchos::RCP<const Core::LinAlg::Vector<double>> lagrmult =
      Teuchos::RCP(new Core::LinAlg::Vector<double>(
          structure_->get_constraint_manager()->get_lagr_mult_vector_old()->Map(), true));

  // merge stuff together
  Teuchos::RCP<Core::LinAlg::Vector<double>> mergedstat =
      Teuchos::RCP(new Core::LinAlg::Vector<double>(*dofrowmap_, true));
  conmerger_->add_cond_vector(*strudis, *mergedstat);
  conmerger_->add_other_vector(*lagrmult, *mergedstat);

  return mergedstat;
}

/*----------------------------------------------------------------------*/
/* non-overlapping DOF map */
Teuchos::RCP<const Epetra_Map> Adapter::StructureConstrMerged::dof_row_map() { return dofrowmap_; }


/*----------------------------------------------------------------------*/
/* stiffness, i.e. force residual R_{n+1} differentiated
 * by displacements D_{n+1} */
Teuchos::RCP<Core::LinAlg::SparseMatrix> Adapter::StructureConstrMerged::system_matrix()
{
  // create empty large matrix and get small ones from structure and constraints
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mergedmatrix =
      Teuchos::RCP(new Core::LinAlg::SparseMatrix(*dofrowmap_, 81));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> strustiff = structure_->system_matrix();
  strustiff->complete();

  Teuchos::RCP<Core::LinAlg::SparseOperator> constiff =
      structure_->get_constraint_manager()->get_constr_matrix();
  constiff->complete();

  // Add matrices together
  mergedmatrix->add(*strustiff, false, 1.0, 0.0);
  mergedmatrix->add(*constiff, false, 1.0, 1.0);
  mergedmatrix->add(*constiff, true, 1.0, 1.0);
  mergedmatrix->complete(*dofrowmap_, *dofrowmap_);

  mergedmatrix->apply_dirichlet(*(structure_->get_dbc_map_extractor()->cond_map()));

  return mergedmatrix;
}


/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>
Adapter::StructureConstrMerged::block_system_matrix()
{
  FOUR_C_THROW("constrained BlockSparseMatrix never to be implemented");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/* build linear system stiffness matrix and rhs/force residual
 *
 * Monolithic FSI accesses the linearised structure problem. */
void Adapter::StructureConstrMerged::evaluate(
    Teuchos::RCP<const Core::LinAlg::Vector<double>> dispstepinc)
{
  // 'initialize' structural displacement as null-pointer
  Teuchos::RCP<Core::LinAlg::Vector<double>> dispstructstepinc = Teuchos::null;

  // Compute residual increments, update total increments and update lagrange multipliers
  if (dispstepinc != Teuchos::null)
  {
    // Extract increments for lagr multipliers and do update
    Teuchos::RCP<Core::LinAlg::Vector<double>> lagrincr =
        conmerger_->extract_other_vector(*dispstepinc);
    structure_->update_iter_incr_constr(lagrincr);
    dispstructstepinc = conmerger_->extract_cond_vector(*dispstepinc);
  }
  // Hand down incremental displacements,
  // structure_ will compute the residual increments on its own
  structure_->evaluate(dispstructstepinc);
}


/*----------------------------------------------------------------------*/
/* domain map */
const Epetra_Map& Adapter::StructureConstrMerged::domain_map() const
{
  return *(Core::LinAlg::merge_map(structure_->domain_map(),
      *(structure_->get_constraint_manager()->get_constraint_map()), false));
}

/*----------------------------------------------------------------------*/
// Apply interface forces
void Adapter::StructureConstrMerged::apply_interface_forces_temporary_deprecated(
    Teuchos::RCP<Core::LinAlg::Vector<double>> iforce)
{
  // create vector with displacement and constraint DOFs
  Teuchos::RCP<Core::LinAlg::Vector<double>> fifc =
      Core::LinAlg::create_vector(*dof_row_map(), true);

  // insert interface forces
  interface_->add_fsi_cond_vector(*iforce, *fifc);

  // extract the force values from the displacement DOFs only
  Teuchos::RCP<Core::LinAlg::Vector<double>> fifcdisp =
      Core::LinAlg::create_vector(*conmerger_->cond_map(), true);
  conmerger_->extract_cond_vector(*fifc, *fifcdisp);

  // set interface forces within the structural time integrator
  set_force_interface(fifcdisp->get_ptr_of_Epetra_MultiVector());

  prepare_partition_step();

  return;
}

FOUR_C_NAMESPACE_CLOSE
