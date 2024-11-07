// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_str_fsiwrapper_immersed.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fsi_str_model_evaluator_partitioned.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_structure_aux.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::FSIStructureWrapperImmersed::FSIStructureWrapperImmersed(
    std::shared_ptr<Structure> structure)
    : FPSIStructureWrapper(structure)
{
  // immersed_ale fsi part
  std::vector<std::shared_ptr<const Epetra_Map>> vecSpaces;

  vecSpaces.push_back(interface_->fsi_cond_map());       // fsi
  vecSpaces.push_back(interface_->immersed_cond_map());  // immersed

  combinedmap_ = Core::LinAlg::MultiMapExtractor::merge_maps(vecSpaces);

  // full blockmap
  Core::LinAlg::MultiMapExtractor blockrowdofmap;
  blockrowdofmap.setup(*combinedmap_, vecSpaces);

  combinedinterface_ = std::make_shared<Core::LinAlg::MapExtractor>(
      *combinedmap_, interface_->fsi_cond_map(), interface_->immersed_cond_map());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FSIStructureWrapperImmersed::apply_immersed_interface_forces(
    std::shared_ptr<Core::LinAlg::Vector<double>> iforce_fsi,
    std::shared_ptr<Core::LinAlg::Vector<double>> iforce_immersed)
{
  fsi_model_evaluator()->get_interface_force_np_ptr()->PutScalar(0.0);

  if (iforce_fsi != nullptr)
    interface_->add_fsi_cond_vector(
        *iforce_fsi, *fsi_model_evaluator()->get_interface_force_np_ptr());
  if (iforce_immersed != nullptr)
    interface_->add_immersed_cond_vector(
        *iforce_immersed, *fsi_model_evaluator()->get_interface_force_np_ptr());

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>>
Adapter::FSIStructureWrapperImmersed::extract_immersed_interface_dispnp()
{
  FOUR_C_ASSERT(interface_->full_map()->PointSameAs(dispnp()->Map()),
      "Full map of map extractor and Dispnp() do not match.");

  return interface_->extract_immersed_cond_vector(*dispnp());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FSIStructureWrapperImmersed::output(
    bool forced_writerestart, const int step, const double time)
{
  // always write velocity and displacement for extra output
  bool writevelacc_ = true;

  // write standard output if no arguments are provided (default -1)
  if (step == -1 and time == -1.0) structure_->output(forced_writerestart);
  // write extra output for specified step and time
  else
  {
    if (structure_->discretization()->get_comm().MyPID() == 0)
      std::cout << "\n   Write EXTRA STRUCTURE Output Step=" << step << " Time=" << time
                << " ...   \n"
                << std::endl;


    structure_->disc_writer()->new_step(step, time);
    structure_->disc_writer()->write_vector("displacement", structure_->dispnp());

    // for visualization of vel and acc do not forget to comment in corresponding lines in
    // StructureEnsightWriter
    if (writevelacc_)
    {
      structure_->disc_writer()->write_vector("velocity", structure_->velnp());
      structure_->disc_writer()->write_vector("acceleration", structure_->accnp());
    }
  }  // write extra output
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Solid::Dbc& Adapter::FSIStructureWrapperImmersed::get_dbc()
{
  return std::dynamic_pointer_cast<Solid::TimeInt::Base>(structure_)->get_dbc();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FSIStructureWrapperImmersed::add_dirich_dofs(
    const std::shared_ptr<const Epetra_Map> maptoadd)
{
  get_dbc().add_dirich_dofs(maptoadd);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FSIStructureWrapperImmersed::remove_dirich_dofs(
    const std::shared_ptr<const Epetra_Map> maptoremove)
{
  get_dbc().remove_dirich_dofs(maptoremove);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FSIStructureWrapperImmersed::set_state(
    const std::shared_ptr<Core::LinAlg::Vector<double>>& x)
{
  return std::dynamic_pointer_cast<Solid::TimeInt::Implicit>(structure_)->set_state(x);
}

FOUR_C_NAMESPACE_CLOSE
