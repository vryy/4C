/*-----------------------------------------------------------*/
/*!
\file ssi_str_model_evaluator_partitioned.cpp

\brief Model evaluator for structure part of partitioned ssi.

\maintainer Andreas Rauch

\date Nov 11, 2016

\level 3

*/
/*-----------------------------------------------------------*/


#include "ssi_str_model_evaluator_partitioned.H"

#include "../drt_structure_new/str_dbc.H"
#include "../drt_structure_new/str_impl_generic.H"
#include "../drt_structure_new/str_timint_implicit.H"
#include "../drt_structure_new/str_nln_solver_generic.H"
#include "../drt_structure_new/str_timint_basedataglobalstate.H"

#include "../solver_nonlin_nox/nox_nln_group.H"

#include "../linalg/linalg_utils.H"

#include "Epetra_Comm.h"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::PartitionedSSI::PartitionedSSI()
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::PartitionedSSI::
    Setup()
{
  // set flag
  issetup_ = true;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::PartitionedSSI::
    GetBlockDofRowMapPtr() const
{
  CheckInitSetup();
  return GState().DofRowMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::PartitionedSSI::
    GetCurrentSolutionPtr() const
{
  CheckInit();
  return GState().GetDisNp();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::PartitionedSSI::
    GetLastTimeStepSolutionPtr() const
{
  CheckInit();
  return GState().GetDisN();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::PartitionedSSI::
    AssembleForce(Epetra_Vector& f,
      const double & timefac_np) const
{
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::PartitionedSSI::
     UpdateStepState(
    const double& timefac_n)
{
  return;
}
