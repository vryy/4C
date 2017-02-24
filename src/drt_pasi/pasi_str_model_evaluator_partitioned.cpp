/*!----------------------------------------------------------------------
\file pasi_str_model_evaluator_partitioned.cpp

\brief model evaluator for structure part of partitioned pasi

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                               sfuchs 01/2017 |
 *----------------------------------------------------------------------*/
#include "pasi_str_model_evaluator_partitioned.H"

#include "../drt_structure_new/str_dbc.H"
#include "../drt_structure_new/str_utils.H"
#include "../drt_structure_new/str_impl_generic.H"
#include "../drt_structure_new/str_timint_implicit.H"
#include "../drt_structure_new/str_nln_solver_generic.H"
#include "../drt_structure_new/str_timint_basedataglobalstate.H"

#include "../solver_nonlin_nox/nox_nln_group.H"

#include "../linalg/linalg_utils.H"

#include "Epetra_Comm.h"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::PartitionedPASI::PartitionedPASI()
{
  // empty
} // STR::MODELEVALUATOR::PartitionedPASI::PartitionedPASI()

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::PartitionedPASI::Setup()
{
  // set flag
  issetup_ = true;

  return;
} // STR::MODELEVALUATOR::PartitionedPASI::Setup()

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::PartitionedPASI::
    GetBlockDofRowMapPtr() const
{
  CheckInitSetup();
  return GState().DofRowMap();
} // STR::MODELEVALUATOR::PartitionedPASI::GetBlockDofRowMapPtr()

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::PartitionedPASI::
    GetCurrentSolutionPtr() const
{
  CheckInit();
  return GState().GetDisNp();
} // STR::MODELEVALUATOR::PartitionedPASI::GetCurrentSolutionPtr()

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::PartitionedPASI::
    GetLastTimeStepSolutionPtr() const
{
  CheckInit();
  return GState().GetDisN();
} // STR::MODELEVALUATOR::PartitionedPASI::GetLastTimeStepSolutionPtr()

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::PartitionedPASI::
    AssembleForce(Epetra_Vector& f,
      const double & timefac_np) const
{
  return true;
} // STR::MODELEVALUATOR::PartitionedPASI::AssembleForce()

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::PartitionedPASI::
     UpdateStepState(
    const double& timefac_n)
{
  return;
} // STR::MODELEVALUATOR::PartitionedPASI::UpdateStepState()
