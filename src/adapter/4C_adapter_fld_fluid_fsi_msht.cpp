// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_fld_fluid_fsi_msht.hpp"

#include "4C_fluid_utils_mapextractor.hpp"

FOUR_C_NAMESPACE_OPEN

/*======================================================================*/
/* constructor */
Adapter::FluidFSIMsht::FluidFSIMsht(Teuchos::RCP<Fluid> fluid,
    Teuchos::RCP<Core::FE::Discretization> dis, Teuchos::RCP<Core::LinAlg::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Core::IO::DiscretizationWriter> output, bool isale, bool dirichletcond)
    : FluidFSI(fluid, dis, solver, params, output, isale, dirichletcond),
      fsiinterface_(Teuchos::make_rcp<FLD::Utils::FsiMapExtractor>())
{
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FluidFSIMsht::init()
{
  // call base class init
  FluidFSI::init();

  // create fluid map extractor
  setup_fsi_interface();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Adapter::FluidFSIMsht::setup_fsi_interface() { fsiinterface_->setup(*dis_); }

FOUR_C_NAMESPACE_CLOSE
