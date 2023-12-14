/*---------------------------------------------------------------------------*/
/*! \file
\brief model evaluator for structure part of partitioned pasi
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_pasi_str_model_evaluator_partitioned.H"

#include "baci_linalg_utils_sparse_algebra_assemble.H"
#include "baci_solver_nonlin_nox_group.H"
#include "baci_structure_new_dbc.H"
#include "baci_structure_new_impl_generic.H"
#include "baci_structure_new_nln_solver_generic.H"
#include "baci_structure_new_timint_basedataglobalstate.H"
#include "baci_structure_new_timint_implicit.H"

#include <Epetra_Comm.h>

BACI_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
STR::MODELEVALUATOR::PartitionedPASI::PartitionedPASI()
{
  // empty constructor
}

void STR::MODELEVALUATOR::PartitionedPASI::Setup()
{
  // pasi interface force at t_{n+1}
  interface_force_np_ptr_ = Teuchos::rcp(new Epetra_Vector(*GState().DofRowMap(), true));

  // set flag
  issetup_ = true;
}

Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::PartitionedPASI::GetBlockDofRowMapPtr() const
{
  CheckInitSetup();
  return GState().DofRowMap();
}

Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::PartitionedPASI::GetCurrentSolutionPtr()
    const
{
  CheckInit();
  return GState().GetDisNp();
}

Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::PartitionedPASI::GetLastTimeStepSolutionPtr()
    const
{
  CheckInit();
  return GState().GetDisN();
}

bool STR::MODELEVALUATOR::PartitionedPASI::AssembleForce(
    Epetra_Vector& f, const double& timefac_np) const
{
  CORE::LINALG::AssembleMyVector(1.0, f, -timefac_np, *interface_force_np_ptr_);

  return true;
}

void STR::MODELEVALUATOR::PartitionedPASI::UpdateStepState(const double& timefac_n) { return; }

BACI_NAMESPACE_CLOSE
