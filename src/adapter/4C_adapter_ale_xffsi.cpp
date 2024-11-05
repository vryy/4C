// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_ale_xffsi.hpp"

#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Adapter::AleXFFsiWrapper::AleXFFsiWrapper(std::shared_ptr<Ale> ale) : AleFsiWrapper(ale)
{
  // create the FSI interface
  xff_interface_ = std::make_shared<ALE::Utils::XFluidFluidMapExtractor>();
  xff_interface_->setup(*discretization());
  setup_dbc_map_ex(ALE::Utils::MapExtractor::dbc_set_x_ff, interface(), xff_interface_);
  setup_dbc_map_ex(ALE::Utils::MapExtractor::dbc_set_x_fsi, interface());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::MapExtractor> Adapter::AleXFFsiWrapper::get_dbc_map_extractor()
{
  return AleWrapper::get_dbc_map_extractor(ALE::Utils::MapExtractor::dbc_set_x_ff);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Adapter::AleXFFsiWrapper::evaluate(std::shared_ptr<const Core::LinAlg::Vector<double>> stepinc)
{
  AleFsiWrapper::evaluate(stepinc, ALE::Utils::MapExtractor::dbc_set_x_ff);

  // set dispnp_ of xfem dofs to dispn_
  xff_interface_->insert_xfluid_fluid_cond_vector(
      *xff_interface_->extract_xfluid_fluid_cond_vector(*dispn()), *write_access_dispnp());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int Adapter::AleXFFsiWrapper::solve()
{
  AleFsiWrapper::evaluate(nullptr, ALE::Utils::MapExtractor::dbc_set_x_fsi);

  int err = AleFsiWrapper::solve();

  update_iter();

  return err;
}

FOUR_C_NAMESPACE_CLOSE
