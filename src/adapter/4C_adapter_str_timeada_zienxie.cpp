/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the structural time integration which gives fine grained
       access in the adaptive time marching loop


\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_adapter_str_timeada_zienxie.hpp"

#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::StructureTimeAdaZienXie::integrate_step_auxiliar()
{
  const STR::TimeInt::Base& stm = *stm_;
  const STR::TimeInt::BaseDataGlobalState& gstate = stm.data_global_state();

  // get state vectors of marching integrator
  Teuchos::RCP<const Epetra_Vector> dis = gstate.get_dis_n();    // D_{n}^{A2}
  Teuchos::RCP<const Epetra_Vector> vel = gstate.get_vel_n();    // V_{n}^{A2}
  Teuchos::RCP<const Epetra_Vector> acc = gstate.get_acc_n();    // A_{n}^{A2}
  Teuchos::RCP<const Epetra_Vector> accn = gstate.get_acc_np();  // A_{n+1}^{A2}

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
