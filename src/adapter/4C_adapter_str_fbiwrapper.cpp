/*----------------------------------------------------------------------*/
/*! \file

\brief Structural adapter for Fluid-beam interaction problems


\level 2
*/

#include "4C_adapter_str_fbiwrapper.hpp"

#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_fsi_str_model_evaluator_partitioned.hpp"
#include "4C_global_data.hpp"
#include "4C_structure_aux.hpp"
#include "4C_structure_new_timint_basedataio.hpp"
#include "4C_structure_new_timint_basedataio_runtime_vtk_output.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::FBIStructureWrapper::FBIStructureWrapper(Teuchos::RCP<Structure> structure)
    : FSIStructureWrapper(structure)
{
  const bool is_prestress = Teuchos::getIntegralValue<Inpar::Solid::PreStress>(
                                Global::Problem::instance()->structural_dynamic_params(),
                                "PRESTRESS") != Inpar::Solid::PreStress::none;
  if (is_prestress)
  {
    FOUR_C_THROW("Prestressing for fluid-beam interaction not tested yet.");
  }
  eletypeextractor_ = Teuchos::rcp(new BEAMINTERACTION::UTILS::MapExtractor);
  BEAMINTERACTION::UTILS::SetupEleTypeMapExtractor(structure_->discretization(), eletypeextractor_);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FBIStructureWrapper::extract_interface_veln()
{
  Teuchos::RCP<Epetra_Vector> veli = Teuchos::rcp(new Epetra_Vector(veln()->Map()));
  veli->Update(1.0, *veln(), 0.0);
  return veli;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FBIStructureWrapper::extract_interface_velnp()
{
  Teuchos::RCP<Epetra_Vector> veli = Teuchos::rcp(new Epetra_Vector(velnp()->Map()));
  veli->Update(1.0, *velnp(), 0.0);
  return veli;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FBIStructureWrapper::predict_interface_velnp()
{
  Teuchos::RCP<Epetra_Vector> veli = Teuchos::rcp(new Epetra_Vector(veln()->Map()));
  veli->Update(1.0, *veln(), 0.0);
  return veli;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FBIStructureWrapper::relaxation_solve(
    Teuchos::RCP<Epetra_Vector> iforce)
{
  FOUR_C_THROW("RelaxationSolve not implemented for immersed fluid-beam interaction\n");
  return Teuchos::null;
}
/*------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------*/
void Adapter::FBIStructureWrapper::rebuild_interface() { FOUR_C_THROW("Not implemented yet"); }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FBIStructureWrapper::predict_interface_dispnp()
{
  Teuchos::RCP<Epetra_Vector> disi = Teuchos::rcp(new Epetra_Vector(dispn()->Map()));
  disi->Update(1.0, *dispnp(), 0.0);
  return disi;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FBIStructureWrapper::extract_interface_dispnp()
{
  Teuchos::RCP<Epetra_Vector> disi = Teuchos::rcp(new Epetra_Vector(dispnp()->Map()));
  disi->Update(1.0, *dispnp(), 0.0);
  return disi;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FBIStructureWrapper::extract_interface_dispn()
{
  Teuchos::RCP<Epetra_Vector> disi = Teuchos::rcp(new Epetra_Vector(dispn()->Map()));
  disi->Update(1.0, *dispn(), 0.0);
  return disi;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Apply interface forces
void Adapter::FBIStructureWrapper::apply_interface_forces(Teuchos::RCP<Epetra_Vector> iforce)
{
  fsi_model_evaluator()->get_interface_force_np_ptr()->Update(
      1.0, *iforce, 0.0);  // todo This has to be changed for mixed structure
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIStructureWrapper::setup_multi_map_extractor()
{
  fsi_model_evaluator()->setup_multi_map_extractor();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

Teuchos::RCP<const Solid::TimeInt::ParamsRuntimeOutput> Adapter::FBIStructureWrapper::get_io_data()
{
  return fsi_model_evaluator()->get_in_output().get_runtime_output_params();
}

FOUR_C_NAMESPACE_CLOSE
