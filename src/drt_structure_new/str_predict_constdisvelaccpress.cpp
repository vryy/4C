/*-----------------------------------------------------------*/
/*!

\brief implementation of predictor for either constant displacement, velocity or acceleration

\maintainer Anh-Tu Vuong

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_predict_constdisvelaccpress.H"
#include "str_timint_base.H"
#include "str_impl_generic.H"
#include "str_predict_factory.H"
#include "str_model_evaluator.H"

#include <Epetra_Vector.h>


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

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::ConstDisVelAccPress::Compute(NOX::Abstract::Group& grp)
{
  CheckInitSetup();

  Teuchos::RCP<Epetra_Vector>& disnp_ptr = GlobalState().GetMutableDisNp();
  Teuchos::RCP<Epetra_Vector>& velnp_ptr = GlobalState().GetMutableVelNp();
  Teuchos::RCP<Epetra_Vector>& accnp_ptr = GlobalState().GetMutableAccNp();

  bool ok = true;
  switch (GetType())
  {
    case INPAR::STR::pred_constdis:
    case INPAR::STR::pred_constdispres:
    {
      ImplInt().PredictConstDisConsistVelAcc(*disnp_ptr, *velnp_ptr, *accnp_ptr);
      break;
    }
    case INPAR::STR::pred_constvel:
    {
      ok = ImplInt().PredictConstVelConsistAcc(*disnp_ptr, *velnp_ptr, *accnp_ptr);
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
      dserror("Unsupported ConstDisVelAccPress predictor!");
      break;
    }
  }
  ImplInt().ModelEval().Predict(GetType());

  // If the const predictors failed e.g. due to too little history information,
  // we use the tangdis predictor as fallback predictor.
  if (not ok) tangdis_ptr_->Compute(grp);

  return;
}
