/*-----------------------------------------------------------*/
/*!
\file str_predict_generic.cpp

\brief Generic class for all predictors.

\maintainer Michael Hiermeier

\date Sep 1, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_predict_generic.H"
#include "str_dbc.H"
#include "str_timint_base.H"
#include "str_impl_generic.H"
#include "str_model_evaluator.H"

#include "../solver_nonlin_nox/nox_nln_group.H"

#include "../drt_io/io_pstream.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::PREDICT::Generic::Generic()
    : isinit_(false),
      issetup_(false),
      type_(INPAR::STR::pred_vague),
      implint_ptr_(Teuchos::null),
      dbc_ptr_(Teuchos::null),
      noxparams_ptr_(Teuchos::null)
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::Generic::Init(const enum INPAR::STR::PredEnum& type,
    const Teuchos::RCP<STR::IMPLICIT::Generic>& implint_ptr,
    const Teuchos::RCP<STR::Dbc>& dbc_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataIO>& iodata_ptr,
    const Teuchos::RCP<Teuchos::ParameterList>& noxparams_ptr)
{
  issetup_ = false;

  // initialize the predictor type
  type_ = type;
  implint_ptr_ = implint_ptr;
  dbc_ptr_ = dbc_ptr;
  gstate_ptr_ = gstate_ptr;
  iodata_ptr_ = iodata_ptr;
  noxparams_ptr_ = noxparams_ptr;

  isinit_ = true;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::Generic::PrePredict(NOX::Abstract::Group& grp)
{
  CheckInitSetup();
  Print();
  dbc_ptr_->UpdateLocSysManager();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::Generic::Predict(NOX::Abstract::Group& grp)
{
  CheckInitSetup();
  bool& ispredict = gstate_ptr_->GetMutableIsPredict();
  ispredict = true;

  // pre-process the prediction step
  PrePredict(grp);

  // compute the actual prediction step
  Compute(grp);

  // post-process the prediction step
  PostPredict(grp);

  ispredict = false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::Generic::PostPredict(NOX::Abstract::Group& grp)
{
  CheckInitSetup();
  // ToDo The pressure dofs need an extra treatment.

  Dbc().ApplyDirichletBC(GlobalState().GetTimeNp(),
      GlobalState().GetMutableDisNp(),
      GlobalState().GetMutableVelNp(),
      GlobalState().GetMutableAccNp(),
      false);

  // Create the new solution vector
  Teuchos::RCP<NOX::Epetra::Vector> x_vec =
      GlobalState().CreateGlobalVector(DRT::UTILS::vec_init_current_state,
          ImplInt().ModelEvalPtr());
  // resets all isValid flags
  grp.setX(*x_vec);

  NOX::NLN::Group* nlngrp_ptr = dynamic_cast<NOX::NLN::Group*>(&grp);
  if (nlngrp_ptr == NULL)
    dserror("Group cast failed!");
  // evaluate the right hand side and the jacobian
  implint_ptr_->SetIsPredictorState(true);
  nlngrp_ptr->computeFandJacobian();
  implint_ptr_->SetIsPredictorState(false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::string STR::PREDICT::Generic::Name() const
{
  CheckInit();
  return INPAR::STR::PredEnumString(type_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::Generic::CheckInit() const
{
  if (not IsInit())
    dserror("Call Init() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::Generic::CheckInitSetup() const
{
  if (!IsInit() or !IsSetup())
    dserror("Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::IMPLICIT::Generic>& STR::PREDICT::Generic::ImplIntPtr()
{
  CheckInit();
  return implint_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::IMPLICIT::Generic& STR::PREDICT::Generic::ImplInt()
{
  CheckInit();
  return *implint_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::Dbc>& STR::PREDICT::Generic::DbcPtr()
{
  CheckInit();
  return dbc_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::Dbc& STR::PREDICT::Generic::Dbc()
{
  CheckInit();
  return *dbc_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& STR::PREDICT::Generic::
    GlobalStatePtr()
{
  CheckInit();
  return gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataGlobalState& STR::PREDICT::Generic::GlobalState()
{
  CheckInit();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TIMINT::BaseDataIO>& STR::PREDICT::Generic::IODataPtr()
{
  CheckInit();
  return iodata_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataIO& STR::PREDICT::Generic::IOData()
{
  CheckInit();
  return *iodata_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::BaseDataGlobalState& STR::PREDICT::Generic::GlobalState() const
{
  CheckInit();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Teuchos::ParameterList>& STR::PREDICT::Generic::NoxParamsPtr()
{
  CheckInit();
  return noxparams_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::ParameterList& STR::PREDICT::Generic::NoxParams()
{
  CheckInit();
  return *noxparams_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::Generic::Print() const
{
  CheckInitSetup();
  if (gstate_ptr_->GetMyRank() == 0 and
      iodata_ptr_->GetPrint2ScreenEveryNStep() and
      gstate_ptr_->GetStepN()%iodata_ptr_->GetPrint2ScreenEveryNStep()==0)
  {
    IO::cout << "=== Structural predictor: " << Name().c_str() <<
        " ===" << IO::endl;
  }
}
