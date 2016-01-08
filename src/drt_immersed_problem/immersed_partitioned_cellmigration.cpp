/*!----------------------------------------------------------------------
\file immersed_partitioned_cellmigration.cpp

\brief partitioned immersed cell migration algorithm

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#include "immersed_partitioned_cellmigration.H"


IMMERSED::ImmersedPartitionedCellMigration::ImmersedPartitionedCellMigration(const Teuchos::ParameterList& params, const Epetra_Comm& comm)
  : ImmersedPartitioned(comm)
{
  // get pointer to global problem
  globalproblem_ = DRT::Problem::Instance();

  // construct immersed exchange manager. singleton class that makes immersed variables comfortably accessible from everywhere in the code
  exchange_manager_ = DRT::ImmersedFieldExchangeManager::Instance();

  // initialize the interaction switches
  fluid_interaction_=DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams(),"FLUID_INTERACTION");
  ecm_interaction_=DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams(),"ECM_INTERACTION");
  adhesion_dynamics_=DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams(),"ADHESION_DYNAMICS");

  exchange_manager_->SetIsFluidInteraction(fluid_interaction_);
}
