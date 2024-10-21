// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_str_timeada_zienxie.hpp"

#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::StructureTimeAdaZienXie::integrate_step_auxiliar()
{
  const Solid::TimeInt::Base& stm = *stm_;
  const Solid::TimeInt::BaseDataGlobalState& gstate = stm.data_global_state();

  // get state vectors of marching integrator
  Teuchos::RCP<const Core::LinAlg::Vector<double>> dis = gstate.get_dis_n();    // D_{n}^{A2}
  Teuchos::RCP<const Core::LinAlg::Vector<double>> vel = gstate.get_vel_n();    // V_{n}^{A2}
  Teuchos::RCP<const Core::LinAlg::Vector<double>> acc = gstate.get_acc_n();    // A_{n}^{A2}
  Teuchos::RCP<const Core::LinAlg::Vector<double>> accn = gstate.get_acc_np();  // A_{n+1}^{A2}

  // build ZX displacements D_{n+1}^{ZX}
  // using the second order (or lower) accurate new accelerations
  locerrdisn_->Update(1.0, *dis, stepsize_, *vel, 0.0);
  locerrdisn_->Update(stepsize_ * stepsize_ / 3.0, *acc, stepsize_ * stepsize_ / 6.0, *accn, 1.0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::StructureTimeAdaZienXie::update_auxiliar()
{
  // NOTHING TO UPDATE
}
FOUR_C_NAMESPACE_CLOSE
