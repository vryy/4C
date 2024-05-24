/*---------------------------------------------------------------------------*/
/*! \file
\brief model evaluator for structure part of partitioned pasi
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_pasi_str_model_evaluator_partitioned.hpp"

#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_solver_nonlin_nox_group.hpp"
#include "4C_structure_new_dbc.hpp"
#include "4C_structure_new_impl_generic.hpp"
#include "4C_structure_new_nln_solver_generic.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"
#include "4C_structure_new_timint_implicit.hpp"

#include <Epetra_Comm.h>

FOUR_C_NAMESPACE_OPEN

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
  interface_force_np_ptr_ = Teuchos::rcp(new Epetra_Vector(*GState().dof_row_map(), true));

  // set flag
  issetup_ = true;
}

Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::PartitionedPASI::get_block_dof_row_map_ptr()
    const
{
  check_init_setup();
  return GState().dof_row_map();
}

Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::PartitionedPASI::get_current_solution_ptr()
    const
{
  check_init();
  return GState().GetDisNp();
}

Teuchos::RCP<const Epetra_Vector>
STR::MODELEVALUATOR::PartitionedPASI::get_last_time_step_solution_ptr() const
{
  check_init();
  return GState().GetDisN();
}

bool STR::MODELEVALUATOR::PartitionedPASI::assemble_force(
    Epetra_Vector& f, const double& timefac_np) const
{
  CORE::LINALG::AssembleMyVector(1.0, f, -timefac_np, *interface_force_np_ptr_);

  return true;
}

void STR::MODELEVALUATOR::PartitionedPASI::UpdateStepState(const double& timefac_n) { return; }

FOUR_C_NAMESPACE_CLOSE
