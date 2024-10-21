// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_str_fbiwrapper.hpp"

#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_fsi_str_model_evaluator_partitioned.hpp"
#include "4C_global_data.hpp"
#include "4C_structure_aux.hpp"
#include "4C_structure_new_timint_basedataio.hpp"
#include "4C_structure_new_timint_basedataio_runtime_vtk_output.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

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
  eletypeextractor_ = Teuchos::make_rcp<BEAMINTERACTION::Utils::MapExtractor>();
  BEAMINTERACTION::Utils::setup_ele_type_map_extractor(
      structure_->discretization(), eletypeextractor_);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Vector<double>> Adapter::FBIStructureWrapper::extract_interface_veln()
{
  Teuchos::RCP<Core::LinAlg::Vector<double>> veli =
      Teuchos::make_rcp<Core::LinAlg::Vector<double>>(veln()->Map());
  veli->Update(1.0, *veln(), 0.0);
  return veli;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Vector<double>> Adapter::FBIStructureWrapper::extract_interface_velnp()
{
  Teuchos::RCP<Core::LinAlg::Vector<double>> veli =
      Teuchos::make_rcp<Core::LinAlg::Vector<double>>(velnp()->Map());
  veli->Update(1.0, *velnp(), 0.0);
  return veli;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Vector<double>> Adapter::FBIStructureWrapper::predict_interface_velnp()
{
  Teuchos::RCP<Core::LinAlg::Vector<double>> veli =
      Teuchos::make_rcp<Core::LinAlg::Vector<double>>(veln()->Map());
  veli->Update(1.0, *veln(), 0.0);
  return veli;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Vector<double>> Adapter::FBIStructureWrapper::relaxation_solve(
    Teuchos::RCP<Core::LinAlg::Vector<double>> iforce)
{
  FOUR_C_THROW("RelaxationSolve not implemented for immersed fluid-beam interaction\n");
  return Teuchos::null;
}
/*------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------*/
void Adapter::FBIStructureWrapper::rebuild_interface() { FOUR_C_THROW("Not implemented yet"); }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Vector<double>> Adapter::FBIStructureWrapper::predict_interface_dispnp()
{
  Teuchos::RCP<Core::LinAlg::Vector<double>> disi =
      Teuchos::make_rcp<Core::LinAlg::Vector<double>>(dispn()->Map());
  disi->Update(1.0, *dispnp(), 0.0);
  return disi;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Vector<double>> Adapter::FBIStructureWrapper::extract_interface_dispnp()
{
  Teuchos::RCP<Core::LinAlg::Vector<double>> disi =
      Teuchos::make_rcp<Core::LinAlg::Vector<double>>(dispnp()->Map());
  disi->Update(1.0, *dispnp(), 0.0);
  return disi;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::Vector<double>> Adapter::FBIStructureWrapper::extract_interface_dispn()
{
  Teuchos::RCP<Core::LinAlg::Vector<double>> disi =
      Teuchos::make_rcp<Core::LinAlg::Vector<double>>(dispn()->Map());
  disi->Update(1.0, *dispn(), 0.0);
  return disi;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Apply interface forces
void Adapter::FBIStructureWrapper::apply_interface_forces(
    Teuchos::RCP<Core::LinAlg::Vector<double>> iforce)
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
