// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_ale_fsi_msht.hpp"

#include "4C_ale_utils_mapextractor.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Adapter::AleFsiMshtWrapper::AleFsiMshtWrapper(Teuchos::RCP<Ale> ale) : AleFsiWrapper(ale)
{
  // create the FSI interface
  fsiinterface_ = Teuchos::make_rcp<ALE::Utils::FsiMapExtractor>();
  fsiinterface_->setup(*discretization());

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const ALE::Utils::FsiMapExtractor> Adapter::AleFsiMshtWrapper::fsi_interface() const
{
  return fsiinterface_;
}

FOUR_C_NAMESPACE_CLOSE
