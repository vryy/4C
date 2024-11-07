// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_fld_fluid_ale_xfem.hpp"

#include "4C_adapter_ale_fluid.hpp"
#include "4C_adapter_fld_fluid_xfsi.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_inpar_validparameters.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::FluidAleXFEM::FluidAleXFEM(const Teuchos::ParameterList& prbdyn, std::string condname)
    : FluidAle(prbdyn, condname)
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::FE::Discretization> Adapter::FluidAleXFEM::boundary_discretization()
{
  // returns the boundary discretization
  // REMARK:
  // the returned discretization has to match the structure discretization at the interface coupling
  // (see FSI::Partitioned::Partitioned(const Epetra_Comm& comm) ) therefore return the boundary dis
  // this is similar to the matching of fluid dis and ale dis in case of Adapter::FluidALE

  std::shared_ptr<XFluidFSI> xfluid = std::dynamic_pointer_cast<XFluidFSI>(fluid_field());

  return xfluid->boundary_discretization();
}


/*----------------------------------------------------------------------*/
/// communication object at the struct interface
/*----------------------------------------------------------------------*/
std::shared_ptr<FLD::Utils::MapExtractor> const& Adapter::FluidAleXFEM::struct_interface()
{
  std::shared_ptr<XFluidFSI> xfluid = std::dynamic_pointer_cast<XFluidFSI>(fluid_field());

  return xfluid->struct_interface();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidAleXFEM::nonlinear_solve(std::shared_ptr<Core::LinAlg::Vector<double>> idisp,
    std::shared_ptr<Core::LinAlg::Vector<double>> ivel)
{
  // if we have values at the interface we need to apply them

  // REMARK: for XFLUID idisp = nullptr, ivel = nullptr (called by fsi_fluid_xfem with
  // default nullptr)
  //         for XFSI   idisp != nullptr

  std::shared_ptr<XFluidFSI> xfluid = std::dynamic_pointer_cast<XFluidFSI>(fluid_field());

  // set idispnp in Xfluid
  if (idisp != nullptr) xfluid->apply_struct_mesh_displacement(idisp);

  // set ivelnp in Xfluid
  if (ivel != nullptr) xfluid->apply_struct_interface_velocities(ivel);

  // Update the ale update part
  if (fluid_field()->interface()->au_cond_relevant())
  {
    std::shared_ptr<const Core::LinAlg::Vector<double>> dispnp = fluid_field()->dispnp();
    std::shared_ptr<Core::LinAlg::Vector<double>> audispnp =
        fluid_field()->interface()->extract_au_cond_vector(*dispnp);
    ale_field()->apply_ale_update_displacements(aucoupfa_->master_to_slave(*audispnp));
  }

  // Note: We do not look for moving ale boundaries (outside the coupling
  // interface) on the fluid side. Thus if you prescribe time variable ale
  // Dirichlet conditions the according fluid Dirichlet conditions will not
  // notice.


  ale_field()->solve();
  std::shared_ptr<Core::LinAlg::Vector<double>> fluiddisp =
      ale_to_fluid_field(ale_field()->dispnp());
  fluid_field()->apply_mesh_displacement(fluiddisp);
  fluid_field()->solve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Adapter::FluidAleXFEM::relaxation_solve(
    std::shared_ptr<Core::LinAlg::Vector<double>> idisp, double dt)
{
  FOUR_C_THROW("RelaxationSolve for XFEM useful?");
  std::cout << "WARNING: RelaxationSolve for XFEM useful?" << std::endl;

  // the displacement -> velocity conversion at the interface
  idisp->Scale(1. / dt);

  return fluid_field()->relaxation_solve(idisp);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Adapter::FluidAleXFEM::extract_interface_forces()
{
  std::shared_ptr<XFluidFSI> xfluid = std::dynamic_pointer_cast<XFluidFSI>(fluid_field());
  return xfluid->extract_struct_interface_forces();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Adapter::FluidAleXFEM::extract_interface_velnp()
{
  FOUR_C_THROW("Robin stuff");
  std::shared_ptr<XFluidFSI> xfluid = std::dynamic_pointer_cast<XFluidFSI>(fluid_field());
  return xfluid->extract_struct_interface_velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Adapter::FluidAleXFEM::extract_interface_veln()
{
  std::shared_ptr<XFluidFSI> xfluid = std::dynamic_pointer_cast<XFluidFSI>(fluid_field());
  return xfluid->extract_struct_interface_veln();
}

FOUR_C_NAMESPACE_CLOSE
