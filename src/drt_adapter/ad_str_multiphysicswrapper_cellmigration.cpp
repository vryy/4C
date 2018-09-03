/*----------------------------------------------------------------------*/
/*!
\file ad_str_multiphysicswrapper_cellmigration.cpp

\brief Adapter wrapper for multiphysics in cell migration.

\maintainer Andreas Rauch

\level 1
*/
/*----------------------------------------------------------------------*/


#include "ad_str_multiphysicswrapper_cellmigration.H"
#include "ad_str_fsiwrapper_immersed.H"
#include "ad_str_structalewrapper.H"
#include "ad_str_ssiwrapper.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::MultiphysicsStructureWrapperCellMigration::MultiphysicsStructureWrapperCellMigration(
    Teuchos::RCP<Structure> structure)
    : StructureWrapper(structure),
      ssi_structure_wrapper_(Teuchos::null),
      fsi_structure_wrapper_(Teuchos::null),
      struct_ale_structure_wrapper_(Teuchos::null),
      issetup_(false),
      isinit_(false)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::MultiphysicsStructureWrapperCellMigration::Init(Teuchos::RCP<Structure> ti_strategy)
{
  // reset the setup flag
  SetIsSetup(false);

  // construct the individual structural wrappers
  ssi_structure_wrapper_ = Teuchos::rcp(new SSIStructureWrapper(ti_strategy));
  fsi_structure_wrapper_ = Teuchos::rcp(new FSIStructureWrapperImmersed(ti_strategy));
  struct_ale_structure_wrapper_ = Teuchos::rcp(new StructAleWrapper(ti_strategy));

  // set isinit_ flag true
  SetIsInit(true);

  return 0;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::MultiphysicsStructureWrapperCellMigration::Setup()
{
  // make sure Init(...) was called first
  CheckIsInit();

  // set flag issetup true
  SetIsSetup(true);
}
