// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fpsi.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fpsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_poroelast_utils.hpp"

FOUR_C_NAMESPACE_OPEN

FPSI::FpsiBase::FpsiBase(const Epetra_Comm& comm, const Teuchos::ParameterList& fpsidynparams)
    : AlgorithmBase(comm, fpsidynparams)
{
  // nothing to do ... so far
}


/*----------------------------------------------------------------------*
 | redistribute the FPSI interface                           thon 11/14 |
 *----------------------------------------------------------------------*/
void FPSI::FpsiBase::redistribute_interface()
{
  Global::Problem* problem = Global::Problem::instance();
  const Epetra_Comm& comm = problem->get_dis("structure")->get_comm();
  FPSI::InterfaceUtils* FPSI_UTILS = FPSI::InterfaceUtils::instance();

  if (Core::Communication::num_mpi_ranks(comm) >
      1)  // if we have more than one processor, we need to redistribute at the FPSI interface
  {
    std::shared_ptr<std::map<int, int>> Fluid_PoroFluid_InterfaceMap =
        FPSI_UTILS->get_fluid_poro_fluid_interface_map();
    std::shared_ptr<std::map<int, int>> PoroFluid_Fluid_InterfaceMap =
        FPSI_UTILS->get_poro_fluid_fluid_interface_map();

    FPSI_UTILS->redistribute_interface(
        *problem->get_dis("fluid"), "fpsi_coupling", *PoroFluid_Fluid_InterfaceMap);
    FPSI_UTILS->redistribute_interface(
        *problem->get_dis("ale"), "fpsi_coupling", *PoroFluid_Fluid_InterfaceMap);
    FPSI_UTILS->redistribute_interface(
        *problem->get_dis("porofluid"), "fpsi_coupling", *Fluid_PoroFluid_InterfaceMap);
    FPSI_UTILS->redistribute_interface(
        *problem->get_dis("structure"), "fpsi_coupling", *Fluid_PoroFluid_InterfaceMap);

    // Material pointers need to be reset after redistribution.
    PoroElast::Utils::set_material_pointers_matching_grid(
        *problem->get_dis("structure"), *problem->get_dis("porofluid"));
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
