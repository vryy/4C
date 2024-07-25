/*-----------------------------------------------------------*/
/*! \file

\brief implementation of predictor for either constant displacement, velocity or acceleration


\level 3

*/
/*-----------------------------------------------------------*/


#include "4C_structure_new_predict_constdisvelaccpress.hpp"

#include "4C_structure_new_impl_generic.hpp"
#include "4C_structure_new_model_evaluator_manager.hpp"
#include "4C_structure_new_predict_factory.hpp"
#include "4C_structure_new_timint_base.hpp"

#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Solid::Predict::ConstDisVelAccPress::ConstDisVelAccPress()
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::Predict::ConstDisVelAccPress::setup()
{
  check_init();

  // fallback predictor
  tangdis_ptr_ = Solid::Predict::build_predictor(Inpar::Solid::pred_tangdis);
  tangdis_ptr_->init(
      get_type(), impl_int_ptr(), dbc_ptr(), global_state_ptr(), io_data_ptr(), nox_params_ptr());
  tangdis_ptr_->setup();

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::Predict::ConstDisVelAccPress::compute(::NOX::Abstract::Group& grp)
{
  check_init_setup();

  Teuchos::RCP<Epetra_Vector>& disnp_ptr = global_state().get_dis_np();
  Teuchos::RCP<Epetra_Vector>& velnp_ptr = global_state().get_vel_np();
  Teuchos::RCP<Epetra_Vector>& accnp_ptr = global_state().get_acc_np();

  bool ok = true;
  switch (get_type())
  {
    case Inpar::Solid::pred_constdis:
    case Inpar::Solid::pred_constdispres:
    {
      impl_int().predict_const_dis_consist_vel_acc(*disnp_ptr, *velnp_ptr, *accnp_ptr);
      break;
    }
    case Inpar::Solid::pred_constvel:
    {
      ok = impl_int().predict_const_vel_consist_acc(*disnp_ptr, *velnp_ptr, *accnp_ptr);
      break;
    }
    case Inpar::Solid::pred_constacc:
    {
      ok = impl_int().predict_const_acc(*disnp_ptr, *velnp_ptr, *accnp_ptr);
      break;
    }
    case Inpar::Solid::pred_constdisvelacc:
    case Inpar::Solid::pred_constdisvelaccpres:
    {
      disnp_ptr->Update(1.0, *global_state().get_dis_n(), 0.0);
      velnp_ptr->Update(1.0, *global_state().get_vel_n(), 0.0);
      accnp_ptr->Update(1.0, *global_state().get_acc_n(), 0.0);
      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported ConstDisVelAccPress predictor!");
      break;
    }
  }
  impl_int().model_eval().predict(get_type());

  // If the const predictors failed e.g. due to too little history information,
  // we use the tangdis predictor as fallback predictor.
  if (not ok) tangdis_ptr_->compute(grp);
}

FOUR_C_NAMESPACE_CLOSE
