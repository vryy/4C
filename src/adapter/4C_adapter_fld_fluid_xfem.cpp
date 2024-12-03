// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_fld_fluid_xfem.hpp"

#include "4C_adapter_fld_fluid_xfsi.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::FluidXFEM::FluidXFEM(const Teuchos::ParameterList& prbdyn, std::string condname)
    : fluid_(std::make_shared<Adapter::FluidBaseAlgorithm>(
          prbdyn, Global::Problem::instance()->fluid_dynamic_params(), "fluid", false)
                 ->fluid_field())
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::FE::Discretization> Adapter::FluidXFEM::discretization()
{
  // returns the boundary discretization
  // REMARK:
  // the returned discretization has to match the structure discretization at the interface coupling
  // (see FSI::Partitioned::Partitioned(MPI_Comm comm) ) therefore return the boundary dis
  // this is similar to the matching of fluid dis and ale dis in case of Adapter::FluidALE
  return fluid_field()->discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::FE::Discretization> Adapter::FluidXFEM::boundary_discretization()
{
  // returns the boundary discretization
  // REMARK:
  // the returned discretization has to match the structure discretization at the interface coupling
  // (see FSI::Partitioned::Partitioned(MPI_Comm comm) ) therefore return the boundary dis
  // this is similar to the matching of fluid dis and ale dis in case of Adapter::FluidALE

  std::shared_ptr<XFluidFSI> xfluid = std::dynamic_pointer_cast<XFluidFSI>(fluid_field());

  return xfluid->boundary_discretization();
}


/*----------------------------------------------------------------------*/
/// communication object at the struct interface
/*----------------------------------------------------------------------*/
std::shared_ptr<FLD::Utils::MapExtractor> const& Adapter::FluidXFEM::struct_interface()
{
  // returns the boundary discretization
  // REMARK:
  // the returned discretization has to match the structure discretization at the interface coupling
  // (see FSI::Partitioned::Partitioned(MPI_Comm comm) ) therefore return the boundary dis
  // this is similar to the matching of fluid dis and ale dis in case of Adapter::FluidALE

  std::shared_ptr<XFluidFSI> xfluid = std::dynamic_pointer_cast<XFluidFSI>(fluid_field());

  return xfluid->struct_interface();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidXFEM::prepare_time_step() { fluid_field()->prepare_time_step(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidXFEM::update() { fluid_field()->update(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidXFEM::output() { fluid_field()->statistics_and_output(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Adapter::FluidXFEM::read_restart(int step)
{
  fluid_field()->read_restart(step);
  return fluid_field()->time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidXFEM::nonlinear_solve(std::shared_ptr<Core::LinAlg::Vector<double>> idisp,
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

  fluid_field()->solve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Adapter::FluidXFEM::relaxation_solve(
    std::shared_ptr<Core::LinAlg::Vector<double>> idisp, double dt)
{
  std::cout << "WARNING: RelaxationSolve for XFEM useful?" << std::endl;

  // the displacement -> velocity conversion at the interface
  idisp->Scale(1. / dt);

  return fluid_field()->relaxation_solve(idisp);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Adapter::FluidXFEM::extract_interface_forces()
{
  std::shared_ptr<XFluidFSI> xfluid = std::dynamic_pointer_cast<XFluidFSI>(fluid_field());
  return xfluid->extract_struct_interface_forces();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Adapter::FluidXFEM::extract_interface_velnp()
{
  FOUR_C_THROW("Robin stuff");
  std::shared_ptr<XFluidFSI> xfluid = std::dynamic_pointer_cast<XFluidFSI>(fluid_field());
  return xfluid->extract_struct_interface_velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Adapter::FluidXFEM::extract_interface_veln()
{
  std::shared_ptr<XFluidFSI> xfluid = std::dynamic_pointer_cast<XFluidFSI>(fluid_field());
  return xfluid->extract_struct_interface_veln();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::LinAlg::Vector<double>> Adapter::FluidXFEM::integrate_interface_shape()
{
  return fluid_field()->integrate_interface_shape();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::Utils::ResultTest> Adapter::FluidXFEM::create_field_test()
{
  return fluid_field()->create_field_test();
}

FOUR_C_NAMESPACE_CLOSE
