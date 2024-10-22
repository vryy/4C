// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_str_wrapper.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::StructureNOXCorrectionWrapper::prepare_time_step()
{
  StructureWrapper::prepare_time_step();
  if (disstepinc_ != Teuchos::null) disstepinc_->PutScalar(0.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::StructureNOXCorrectionWrapper::evaluate(
    Teuchos::RCP<const Core::LinAlg::Vector<double>> disstepinc)
{
  // The field solver always expects an iteration increment only. And
  // there are Dirichlet conditions that need to be preserved. So take
  // the sum of increments we get from NOX and apply the latest iteration
  // increment only.
  // Naming:
  //
  // x^n+1_i+1 = x^n+1_i + disiterinc  (sometimes referred to as residual increment), and
  //
  // x^n+1_i+1 = x^n     + disstepinc

  if (disstepinc != Teuchos::null)
  {
    // iteration increments
    Teuchos::RCP<Core::LinAlg::Vector<double>> disiterinc =
        Teuchos::make_rcp<Core::LinAlg::Vector<double>>(*disstepinc);
    if (disstepinc_ != Teuchos::null)
    {
      disiterinc->Update(-1.0, *disstepinc_, 1.0);

      // update incremental displacement member to provided step increments
      // shortly: disinc_^<i> := disp^<i+1>
      disstepinc_->Update(1.0, *disstepinc, 0.0);
    }
    else
    {
      disstepinc_ = Teuchos::make_rcp<Core::LinAlg::Vector<double>>(*disstepinc);
    }

    // do structural update with provided residual displacements - iteration increment
    StructureWrapper::evaluate(disiterinc);
  }
  else
  {
    StructureWrapper::evaluate(Teuchos::null);
  }
}

FOUR_C_NAMESPACE_CLOSE
