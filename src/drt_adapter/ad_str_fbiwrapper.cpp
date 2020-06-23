/*----------------------------------------------------------------------*/
/*! \file

\brief Structural adapter for Fluid-beam interaction problems

\maintainer Nora Hagmeyer

\level 2
*/

#include "ad_str_fbiwrapper.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_beaminteraction/beaminteraction_calc_utils.H"
#include "../drt_fsi/fsi_str_model_evaluator_partitioned.H"
#include "../drt_structure_new/str_timint_basedataio_runtime_vtk_output.H"
#include "../drt_structure_new/str_timint_basedataio.H"
#include "../drt_lib/prestress_service.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FBIStructureWrapper::FBIStructureWrapper(Teuchos::RCP<Structure> structure)
    : FSIStructureWrapper(structure)
{
  if (::UTILS::PRESTRESS::IsAny())
  {
    dserror("Prestressing for fluid-beam interaction not tested yet.");
  }
  eletypeextractor_ = Teuchos::rcp(new BEAMINTERACTION::UTILS::MapExtractor);
  BEAMINTERACTION::UTILS::SetupEleTypeMapExtractor(structure_->Discretization(), eletypeextractor_);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIStructureWrapper::ExtractInterfaceVeln()
{
  Teuchos::RCP<Epetra_Vector> veli = Teuchos::rcp(new Epetra_Vector(Veln()->Map()));
  veli->Update(1.0, *Veln(), 0.0);
  return veli;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIStructureWrapper::ExtractInterfaceVelnp()
{
  Teuchos::RCP<Epetra_Vector> veli = Teuchos::rcp(new Epetra_Vector(Velnp()->Map()));
  veli->Update(1.0, *Velnp(), 0.0);
  return veli;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIStructureWrapper::PredictInterfaceVelnp()
{
  Teuchos::RCP<Epetra_Vector> veli = Teuchos::rcp(new Epetra_Vector(Veln()->Map()));
  veli->Update(1.0, *Veln(), 0.0);
  return veli;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIStructureWrapper::RelaxationSolve(
    Teuchos::RCP<Epetra_Vector> iforce)
{
  dserror("RelaxationSolve not implemented for immersed fluid-beam interaction\n");
  return Teuchos::null;
}
/*------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------*/
void ADAPTER::FBIStructureWrapper::RebuildInterface() { dserror("Not implemented yet"); }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIStructureWrapper::PredictInterfaceDispnp()
{
  Teuchos::RCP<Epetra_Vector> disi = Teuchos::rcp(new Epetra_Vector(Dispn()->Map()));
  disi->Update(1.0, *Dispnp(), 0.0);
  return disi;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIStructureWrapper::ExtractInterfaceDispnp()
{
  Teuchos::RCP<Epetra_Vector> disi = Teuchos::rcp(new Epetra_Vector(Dispnp()->Map()));
  disi->Update(1.0, *Dispnp(), 0.0);
  return disi;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIStructureWrapper::ExtractInterfaceDispn()
{
  Teuchos::RCP<Epetra_Vector> disi = Teuchos::rcp(new Epetra_Vector(Dispn()->Map()));
  disi->Update(1.0, *Dispn(), 0.0);
  return disi;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Apply interface forces
void ADAPTER::FBIStructureWrapper::ApplyInterfaceForces(Teuchos::RCP<Epetra_Vector> iforce)
{
  FSIModelEvaluator()->GetInterfaceForceNpPtr()->Update(
      1.0, *iforce, 0.0);  // todo This has to be changed for mixed structure
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIStructureWrapper::SetupMultiMapExtractor()
{
  FSIModelEvaluator()->SetupMultiMapExtractor();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

Teuchos::RCP<const STR::TIMINT::ParamsRuntimeVtkOutput> ADAPTER::FBIStructureWrapper::GetIOData()
{
  return FSIModelEvaluator()->GetInOutput().GetRuntimeVtkOutputParams();
}
