/*---------------------------------------------------------------------------*/
/*! \file
\brief model evaluator for structure part of partitioned pasi
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "pasi_str_model_evaluator_partitioned.H"

#include "../drt_structure_new/str_dbc.H"
#include "../drt_structure_new/str_impl_generic.H"
#include "../drt_structure_new/str_timint_implicit.H"
#include "../drt_structure_new/str_nln_solver_generic.H"
#include "../drt_structure_new/str_timint_basedataglobalstate.H"

#include "../solver_nonlin_nox/nox_nln_group.H"

#include "../linalg/linalg_utils_sparse_algebra_assemble.H"

#include "Epetra_Comm.h"

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
  LINALG::AssembleMyVector(1.0, f, -timefac_np, *interface_force_np_ptr_);

  return true;
}

void STR::MODELEVALUATOR::PartitionedPASI::UpdateStepState(const double& timefac_n) { return; }
