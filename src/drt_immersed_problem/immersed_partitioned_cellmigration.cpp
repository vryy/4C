/*!----------------------------------------------------------------------

\brief partitioned immersed cell migration algorithm

\level 2

\maintainer Jonas Eichinger

*----------------------------------------------------------------------*/
#include "immersed_partitioned_cellmigration.H"


IMMERSED::ImmersedPartitionedCellMigration::ImmersedPartitionedCellMigration(
    const Teuchos::ParameterList& params, const Epetra_Comm& comm)
    : ImmersedPartitioned(comm)
{
  // get pointer to global problem
  globalproblem_ = DRT::Problem::Instance();

  // construct immersed exchange manager. singleton class that makes immersed variables comfortably
  // accessible from everywhere in the code
  exchange_manager_ = DRT::ImmersedFieldExchangeManager::Instance();

  // initialize the interaction switches
  fluid_interaction_ =
      DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams(), "FLUID_INTERACTION");
  ecm_interaction_ =
      DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams(), "ECM_INTERACTION");
  adhesion_dynamics_ =
      DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams(), "ADHESION_DYNAMICS");

  exchange_manager_->SetIsFluidInteraction(fluid_interaction_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int IMMERSED::ImmersedPartitionedCellMigration::Init(const Teuchos::ParameterList& params)
{
  // reset the setup flag
  SetIsSetup(false);

  // do all init stuff here

  // set isinit_ flag true
  SetIsInit(true);

  return 0;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::Setup()
{
  // make sure Init(...) was called first
  CheckIsInit();

  // do all setup stuff here

  // set flag issetup true
  SetIsSetup(true);
}
