// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_ale_wrapper.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Adapter::AleNOXCorrectionWrapper::prepare_time_step()
{
  AleWrapper::prepare_time_step();

  if (stepinc_ != Teuchos::null) stepinc_->PutScalar(0.0);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Adapter::AleNOXCorrectionWrapper::evaluate(
    Teuchos::RCP<const Core::LinAlg::Vector<double>> stepinc)
{
  if (stepinc != Teuchos::null)
  {
    // iteration increments
    Teuchos::RCP<Core::LinAlg::Vector<double>> iterinc =
        Teuchos::make_rcp<Core::LinAlg::Vector<double>>(*stepinc);
    if (stepinc_ != Teuchos::null)
    {
      iterinc->Update(-1.0, *stepinc_, 1.0);

      // update incremental displacement member to provided step increments
      // shortly: disinc_^<i> := disp^<i+1>
      stepinc_->Update(1.0, *stepinc, 0.0);
    }
    else
    {
      stepinc_ = Teuchos::make_rcp<Core::LinAlg::Vector<double>>(*stepinc);
    }

    // do structural update with provided residual displacements - iteration increment
    AleWrapper::evaluate(iterinc);
  }
  else
  {
    AleWrapper::evaluate(Teuchos::null);
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
