/*-----------------------------------------------------------*/
/*! \file

\brief Generic class for all predictors.


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_predict_generic.hpp"

#include "4C_io_pstream.hpp"
#include "4C_solver_nonlin_nox_group.hpp"
#include "4C_structure_new_dbc.hpp"
#include "4C_structure_new_impl_generic.hpp"
#include "4C_structure_new_model_evaluator.hpp"
#include "4C_structure_new_timint_base.hpp"

FOUR_C_NAMESPACE_OPEN

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
    const Teuchos::RCP<STR::IMPLICIT::Generic>& implint_ptr, const Teuchos::RCP<STR::Dbc>& dbc_ptr,
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
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::Generic::PrePredict(::NOX::Abstract::Group& grp)
{
  CheckInitSetup();
  Print();
  dbc_ptr_->UpdateLocSysManager();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::Generic::Predict(::NOX::Abstract::Group& grp)
{
  CheckInitSetup();
  bool& ispredict = gstate_ptr_->IsPredict();
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
void STR::PREDICT::Generic::PostPredict(::NOX::Abstract::Group& grp)
{
  CheckInitSetup();

  Dbc().ApplyDirichletBC(GlobalState().GetTimeNp(), GlobalState().GetDisNp(),
      GlobalState().GetVelNp(), GlobalState().GetAccNp(), false);

  // Create the new solution vector
  Teuchos::RCP<::NOX::Epetra::Vector> x_vec = GlobalState().CreateGlobalVector(
      TIMINT::BaseDataGlobalState::VecInitType::init_current_state, ImplInt().ModelEvalPtr());
  // resets all isValid flags
  grp.setX(*x_vec);

  NOX::NLN::Group* nlngrp_ptr = dynamic_cast<NOX::NLN::Group*>(&grp);
  FOUR_C_ASSERT(nlngrp_ptr != nullptr, "Group cast failed!");
  // evaluate the right hand side and the jacobian
  implint_ptr_->SetIsPredictorState(true);
  nlngrp_ptr->computeFandJacobian();
  implint_ptr_->SetIsPredictorState(false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::string STR::PREDICT::Generic::Name() const
{
  CheckInit();
  return INPAR::STR::PredEnumString(type_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::Generic::CheckInit() const { FOUR_C_ASSERT(IsInit(), "Call Init() first!"); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::PREDICT::Generic::CheckInitSetup() const
{
  FOUR_C_ASSERT(IsInit() and IsSetup(), "Call Init() and Setup() first!");
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
Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& STR::PREDICT::Generic::GlobalStatePtr()
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
  if (gstate_ptr_->GetMyRank() == 0 and iodata_ptr_->get_print2_screen_every_n_step() and
      gstate_ptr_->GetStepN() % iodata_ptr_->get_print2_screen_every_n_step() == 0)
  {
    IO::cout << "=== Structural predictor: " << Name().c_str() << " ===" << IO::endl;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::PREDICT::Generic::pre_apply_force_external(Epetra_Vector& fextnp) const
{
  // do nothing
  return false;
}

FOUR_C_NAMESPACE_CLOSE
