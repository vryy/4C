/*-----------------------------------------------------------*/
/*! \file

\brief implementation of predictor for either constant displacement, velocity or acceleration


\level 3

*/
/*-----------------------------------------------------------*/


#include "4C_structure_new_predict_constdisvelaccpress.hpp"

#include "4C_structure_new_impl_generic.hpp"
#include "4C_structure_new_model_evaluator.hpp"
#include "4C_structure_new_predict_factory.hpp"
#include "4C_structure_new_timint_base.hpp"

#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::PREDICT::ConstDisVelAccPress::ConstDisVelAccPress()
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::ConstDisVelAccPress::Setup()
{
  CheckInit();

  // fallback predictor
  tangdis_ptr_ = STR::PREDICT::BuildPredictor(INPAR::STR::pred_tangdis);
  tangdis_ptr_->Init(
      GetType(), ImplIntPtr(), DbcPtr(), GlobalStatePtr(), IODataPtr(), NoxParamsPtr());
  tangdis_ptr_->Setup();

  issetup_ = true;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::ConstDisVelAccPress::Compute(::NOX::Abstract::Group& grp)
{
  CheckInitSetup();

  Teuchos::RCP<Epetra_Vector>& disnp_ptr = GlobalState().GetDisNp();
  Teuchos::RCP<Epetra_Vector>& velnp_ptr = GlobalState().GetVelNp();
  Teuchos::RCP<Epetra_Vector>& accnp_ptr = GlobalState().GetAccNp();

  bool ok = true;
  switch (GetType())
  {
    case INPAR::STR::pred_constdis:
    case INPAR::STR::pred_constdispres:
    {
      ImplInt().predict_const_dis_consist_vel_acc(*disnp_ptr, *velnp_ptr, *accnp_ptr);
      break;
    }
    case INPAR::STR::pred_constvel:
    {
      ok = ImplInt().predict_const_vel_consist_acc(*disnp_ptr, *velnp_ptr, *accnp_ptr);
      break;
    }
    case INPAR::STR::pred_constacc:
    {
      ok = ImplInt().PredictConstAcc(*disnp_ptr, *velnp_ptr, *accnp_ptr);
      break;
    }
    case INPAR::STR::pred_constdisvelacc:
    case INPAR::STR::pred_constdisvelaccpres:
    {
      disnp_ptr->Update(1.0, *GlobalState().GetDisN(), 0.0);
      velnp_ptr->Update(1.0, *GlobalState().GetVelN(), 0.0);
      accnp_ptr->Update(1.0, *GlobalState().GetAccN(), 0.0);
      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported ConstDisVelAccPress predictor!");
      break;
    }
  }
  ImplInt().ModelEval().Predict(GetType());

  // If the const predictors failed e.g. due to too little history information,
  // we use the tangdis predictor as fallback predictor.
  if (not ok) tangdis_ptr_->Compute(grp);
}

FOUR_C_NAMESPACE_CLOSE
